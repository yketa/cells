#include <cmath>
#include <numbers>
#include <assert.h>

#include "system.hpp"
#include "tools.hpp"

void VertexModel::integrate(double const& dt,
    double const& delta, double const& epsilon) {

    // get forces

    std::map<long int, std::vector<double>> const forces = getForces();

    // integrate

    long int vertexIndex;
    std::vector<double> uposition;
    for (auto it=forces.begin(); it != forces.end(); ++it) {
        vertexIndex = it->first;
        uposition = vertices.at(vertexIndex).getUPosition();
        for (int dim=0; dim < 2; dim++) {
            uposition[dim] += forces.at(vertexIndex)[dim]*dt;           // Euler integration of position
        }
        vertices[vertexIndex].setUPosition( // unwrapped position
            uposition);
        vertices[vertexIndex].setPosition(  // (wrapped) position
            wrap(vertices[vertexIndex].getUPosition()));
    }
    for (SPVertex* sPVertex : sPVertices.getValues()) {
        sPVertex->settheta(
            sPVertex->gettheta()
                + std::sqrt(2*sPVertex->getDr()*dt)*random.gauss());    // Euler integration of angle
    }

    // move cell centres

    std::vector<double> cellUPosition(0), initialCellPosition(0);
    std::vector<long int> neighbourVerticesIndices(0); int numberNeighbours;
    std::vector<double> disp(0);
    for (Cell* cell : cells.getValues()) {
        vertexIndex = cell->getVertexIndex();

        cellUPosition = vertices.at(vertexIndex).getUPosition();
        initialCellPosition = vertices.at(vertexIndex).getPosition();

        neighbourVerticesIndices = getNeighbourVertices(vertexIndex)[0];
        numberNeighbours = neighbourVerticesIndices.size();
        for (long int neighbourVertexIndex : neighbourVerticesIndices) {
            disp = wrapDiff(
                initialCellPosition,
                vertices.at(neighbourVertexIndex).getPosition());

            for (int dim=0; dim < 2; dim++) {
                cellUPosition[dim] += disp[dim]/numberNeighbours;
            }
        }
        vertices[vertexIndex].setUPosition( // unwrapped position
            cellUPosition);
        vertices[vertexIndex].setPosition(  // (wrapped) position
            wrap(vertices.at(vertexIndex).getUPosition()));
    }

    // perform T1s

    doT1(delta, epsilon);

    // update time

    time += dt;
}

std::map<long int,std::vector<double>> const VertexModel::getForces() {

    // VERTEX SELF-PROPULSION FORCE

    long int vertexIndex;
    std::map<long int, std::vector<double>> sPForces;
    double sp; std::vector<double> meanSP(2, 0);
    long int const nSPVertices = sPVertices.size();
    for (SPVertex* sPVertex : sPVertices.getValues()) { // loop over self-propelled vertices
        vertexIndex = sPVertex->getVertexIndex();
        sPForces[vertexIndex] = {0, 0};

        for (int dim=0; dim < 2; dim++) {

            sp = sPVertex->getv0()
                *std::cos(sPVertex->gettheta() - dim*std::numbers::pi/2);
            sPForces[vertexIndex][dim] = sp;            // add self-propulsion force
            meanSP[dim] += sp/nSPVertices;              // compute average self-propulsion force
        }
    }

    for (SPVertex* sPVertex : sPVertices.getValues()) {
        vertexIndex = sPVertex->getVertexIndex();
        for (int dim=0; dim < 2; dim++) {
            sPForces[vertexIndex][dim] -= meanSP[dim];  // remove average self-propulsion force
        }
    }

    // VERTEX MODEL INTERACTION FORCE

    std::map<long int, std::vector<double>> forces;
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        forces[it->first] = {0, 0};                                 // initialise zero force at each vertex
        if (sPVertices.in(it->first)) {
            for (int dim=0; dim < 2; dim++) {
                forces[it->first][dim] += sPForces[it->first][dim]; // add self-propulsion forces to self-propelled vertices
            }
        }
    }

    long int cellVertexIndex;
    std::vector<long int> neighbourVerticesIndices(0); int numberNeighbours;
    long int neighbourVertexIndex;
    long int previousNeighbourVertexIndex, nextNeighbourVertexIndex;
    std::vector<double> fromTo(0), fromToBis(0);
    std::vector<double> crossFromTo(0), crossFromToBis(0);
    for (Cell* cell : cells.getValues()) {          // loop over cells
        cellVertexIndex = cell->getVertexIndex();

        cell->setArea(      // update area
            getVertexToNeighboursArea(cellVertexIndex));
        cell->setPerimeter( // update perimeter
            getVertexToNeighboursPerimeter(cellVertexIndex));

        neighbourVerticesIndices = getNeighbourVertices(cellVertexIndex)[0];    // indices of neighbours of cell centre (i.e. cell corners) are in anticlockwise order
        numberNeighbours = neighbourVerticesIndices.size();
        for (int i=0; i < numberNeighbours; i++) {  // loop over vertices in the cell
            neighbourVertexIndex = neighbourVerticesIndices[i];
            assert(!cells.in(neighbourVertexIndex));

            previousNeighbourVertexIndex =
                neighbourVerticesIndices[pmod(i - 1, numberNeighbours)];
            nextNeighbourVertexIndex =
                neighbourVerticesIndices[pmod(i + 1, numberNeighbours)];

            // area term

            fromTo =    // vector from cell to previous vertex
                wrapTo(cellVertexIndex, previousNeighbourVertexIndex,
                    false);
            crossFromTo = cross2z(fromTo);
            fromToBis = // vector from cell to next vertex
                wrapTo(cellVertexIndex, nextNeighbourVertexIndex,
                    false);
            crossFromToBis = cross2z(fromToBis);
            for (int dim=0; dim < 2; dim++) {

                forces[neighbourVertexIndex][dim] +=
                    (cell->getkA()/2.)*(cell->getArea() - cell->getA0())*(
                        crossFromTo[dim] - crossFromToBis[dim]);
            }

            // perimeter term

            fromTo =    // unit vector from vertex to previous vertex
                wrapTo(neighbourVertexIndex, previousNeighbourVertexIndex,
                    true);
            fromToBis = // unit vector from vertex to next vertex
                wrapTo(neighbourVertexIndex, nextNeighbourVertexIndex,
                    true);
            for (int dim=0; dim < 2; dim++) {

                forces[neighbourVertexIndex][dim] +=
                    cell->getkP()*(cell->getPerimeter() - cell->getP0())*(
                        fromTo[dim] + fromToBis[dim]);
            }
        }
    }

    return forces;
}

void VertexModel::doT1(double const& delta, double const& epsilon) {

    // identify small junctions

    std::vector<long int> halfEdgeIndices(0);
    long int fromMergeIndex, toMergeIndex;
    for (Junction* junction : junctions.getValues()) {
        if (getBoundary()) {
            fromMergeIndex =            // (first) vertex to be (potentially) merged into neighbour
                halfEdges.at(junction->getHalfEdgeIndex()).getFromIndex();
            toMergeIndex =              // (second) vertex towards which neighbour is (potentially) merged
                halfEdges.at(junction->getHalfEdgeIndex()).getToIndex();
            if (vertices.at(fromMergeIndex).getBoundary()
                || vertices.at(toMergeIndex).getBoundary()) {
                continue;               // boundary vertex cannot be merged
            }
        }
        if (getEdgeLength(junction->getHalfEdgeIndex()) < delta) {
            halfEdgeIndices.push_back(junction->getHalfEdgeIndex());
        }
    }
    random.shuffle(halfEdgeIndices);    // do T1s in a random order

    // perform T1s

    std::vector<std::vector<long int>> neighbours;
    std::vector<long int> neighboursFromMerge, halfEdgesNeighboursFromMerge;
    std::vector<long int> neighboursToMerge, halfEdgesNeighboursToMerge;
    int numberNeighbours;
    std::vector<long int> halfEdgeToCellsIndices;
    long int vertexIndex;
    long int createHalfEdgeIndex0, createHalfEdgeIndex1;
    double angle;
    for (long int mergeHalfEdgeIndex : halfEdgeIndices) {
        nT1++;

        // identify half-edge to split to create new junction

        fromMergeIndex = halfEdges.at(mergeHalfEdgeIndex).getFromIndex();   // (first) vertex to be merged into neighbour
        toMergeIndex = halfEdges.at(mergeHalfEdgeIndex).getToIndex();       // (second) vertex towards which neighbour is merged

        neighbours = getNeighbourVertices(fromMergeIndex);  // neighbours of first vertex and half-edges towards them
        neighboursFromMerge = neighbours[0];
        halfEdgesNeighboursFromMerge = neighbours[1];

        neighbours = getNeighbourVertices(toMergeIndex);    // neighbours of second vertex and half-edges towards them
        neighboursToMerge = neighbours[0];
        halfEdgesNeighboursToMerge = neighbours[1];

        numberNeighbours = neighboursFromMerge.size();
        assert(numberNeighbours == (int) halfEdgesNeighboursFromMerge.size());
        halfEdgeToCellsIndices.clear();
        for (int i=0; i < numberNeighbours; i++) {
            vertexIndex = neighboursFromMerge[i];
            if (cells.in(vertexIndex)                                   // cell centres neighbours of the first vertex...
                || vertices.at(vertexIndex).getBoundary()) {            // ... or boundary neighbours of the first vertex...
                if (!inVec(neighboursToMerge, vertexIndex)) {           // ... which are not neighbours to the second vertex
                    halfEdgeToCellsIndices.push_back(
                        halfEdgesNeighboursFromMerge[i]);
                }
            }
        }
        assert(halfEdgeToCellsIndices.size() > 0);
        createHalfEdgeIndex0 = random.pick(halfEdgeToCellsIndices);    // randomly pick one

        numberNeighbours = neighboursToMerge.size();
        assert(numberNeighbours == (int) halfEdgesNeighboursToMerge.size());
        halfEdgeToCellsIndices.clear();
        for (int i=0; i < numberNeighbours; i++) {
            vertexIndex = neighboursToMerge[i];
            if (cells.in(vertexIndex)                                   // cell centres neighbours of the second vertex...
                || vertices.at(vertexIndex).getBoundary()) {            // ... or boundary neighbours of the second vertex...
                if (!inVec(neighboursFromMerge, vertexIndex)) {         // ... which are not neighbours to the first vertex
                    halfEdgeToCellsIndices.push_back(
                        halfEdgesNeighboursToMerge[i]);
                }
            }
        }
        assert(halfEdgeToCellsIndices.size() > 0);
        createHalfEdgeIndex1 = random.pick(halfEdgeToCellsIndices);    // randomly pick one

        angle = std::numbers::pi/2. // create new junction orthogonal to the deleted junction
            + angle2(getHalfEdgeVector(mergeHalfEdgeIndex));

//         // output T1 to standard error
// 
//         long int cellA0, cellB0, cellA1, cellB1;
//         cellA0 = halfEdges[                                     // first cell for which junction is created
//             halfEdges[mergeHalfEdgeIndex]
//                 .getPreviousIndex()]
//                     .getFromIndex();
//         cellB0 = halfEdges[                                     // second cell for which junction is deleted
//             halfEdges[halfEdges[mergeHalfEdgeIndex].getPairIndex()]
//                     .getPreviousIndex()]
//                         .getFromIndex();
//         cellA1 = halfEdges[createHalfEdgeIndex0].getToIndex();  // first cell for which junction is created
//         cellB1 = halfEdges[createHalfEdgeIndex1].getToIndex();  // second cell for whch junction is created
//         std::cerr << "T1: #" << mergeHalfEdgeIndex
//             << " ("
//             << std::min(cellA0, cellB0) << "|" << std::max(cellA0, cellB0)
//             << " to "
//             << std::min(cellA1, cellB1) << "|" << std::max(cellA1, cellB1)
//             << ")" << std::endl;

        // merge vertices

        mergeVertices(mergeHalfEdgeIndex);

        // create new vertex

        createJunction(createHalfEdgeIndex0, createHalfEdgeIndex1,
            angle, delta + epsilon);
    }
//     if (halfEdgeIndices.size() > 0) { checkMesh(); }
}

long int const VertexModel::mergeVertices(long int const& halfEdgeIndex) {

    // identify vertices to merge

    long int const fromMergeIndex = // (first) vertex to be merge into neighbour
        halfEdges.at(halfEdgeIndex).getFromIndex();
    long int const toMergeIndex =   // (second) vertex towards which neighbour is merged
        halfEdges.at(halfEdgeIndex).getToIndex();

    // relabel half-edges origins and destinations

    long int previousHalfEdgeIndex =            // half-edge pointing to the vertex to be removed in the first triangle to be removed
        halfEdges.at(halfEdgeIndex).getPreviousIndex();
    long int const endPreviousHalfEdgeIndex =   // half-edge pointing to the vertex to be removed in the second triangle to be removed
        halfEdges.at(halfEdgeIndex).getPairIndex();

    long int pairHalfEdgeIndex;
    while (true) {

        pairHalfEdgeIndex = halfEdges.at(previousHalfEdgeIndex).getPairIndex(); // half-edge pointing from the vertex to be removed and to be relabelled

        // half-edge coming from vertex to be removed
        assert(
            fromMergeIndex
                == halfEdges.at(pairHalfEdgeIndex).getFromIndex());
        halfEdges[pairHalfEdgeIndex].setFromIndex(toMergeIndex);

        // half-edge going to vertex to be removed
        assert(
            fromMergeIndex
                == halfEdges.at(previousHalfEdgeIndex).getToIndex());
        halfEdges[previousHalfEdgeIndex].setToIndex(toMergeIndex);

        if (previousHalfEdgeIndex == endPreviousHalfEdgeIndex) {    // all triangles which the vertex to be removed belongs to have been explored and half-edges origins and destinations relabelled
            assert(pairHalfEdgeIndex == halfEdgeIndex);
            break;
        }

        previousHalfEdgeIndex = // half-edge pointing to the vertex to be removed and to be relabelled
            halfEdges.at(pairHalfEdgeIndex).getPreviousIndex();
    }

    // reassign identifying half-edges of points belonging to triangles to be removed

    long int nextHalfEdgeIndex, thirdVertexIndex;

    previousHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPreviousIndex(); // half-edge pointing to the vertex to be removed in the first triangle to be removed
    pairHalfEdgeIndex = halfEdges.at(previousHalfEdgeIndex).getPairIndex(); // half-edge pointing from the vertex to be removed

    assert(toMergeIndex == halfEdges.at(pairHalfEdgeIndex).getFromIndex());
    vertices[toMergeIndex].setHalfEdgeIndex(pairHalfEdgeIndex);

    nextHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getNextIndex();         // half-edge pointing from the merge vertex in the first triangle to be removed
    pairHalfEdgeIndex = halfEdges.at(nextHalfEdgeIndex).getPairIndex();     // half-edge pointing to the merge vertex
    thirdVertexIndex = halfEdges.at(pairHalfEdgeIndex).getFromIndex();      // index of third vertex in the first triangle to be removed

    vertices[thirdVertexIndex].setHalfEdgeIndex(pairHalfEdgeIndex);

    pairHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPairIndex();         // half-edge pointing to the vertex to be removed in the second triangle to be removed
    nextHalfEdgeIndex = halfEdges.at(pairHalfEdgeIndex).getNextIndex();     // half-edge pointing from the vertex to be removed in the second triangle to be removed
    pairHalfEdgeIndex = halfEdges.at(nextHalfEdgeIndex).getPairIndex();     // half-edge pointing to the vertex to be removed
    thirdVertexIndex = halfEdges.at(pairHalfEdgeIndex).getFromIndex();      // index of third vertex in the second triangle to be removed

    vertices[thirdVertexIndex].setHalfEdgeIndex(pairHalfEdgeIndex);

    // relabel half-edge pairs and delete/merge junctions

    long int fromHalfEdgeIndex, toHalfEdgeIndex;

    fromHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPreviousIndex();     // half-edge pointing to the vertex to be removed in the first triangle to be removed
    fromHalfEdgeIndex = halfEdges.at(fromHalfEdgeIndex).getPairIndex();     // half-edge pointing from the vertex to be removed paired to the first triangle to be removed

    toHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getNextIndex();           // half-edge pointing from the merge vertex in the first triangle to be removed
    toHalfEdgeIndex = halfEdges.at(toHalfEdgeIndex).getPairIndex();         // half-edge pointing to the vertex to be removed paired to the first triangle to be removed

    halfEdges[fromHalfEdgeIndex].setPairIndex(toHalfEdgeIndex);
    halfEdges[toHalfEdgeIndex].setPairIndex(fromHalfEdgeIndex);

    fromHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPairIndex();         // half-edge pointing from the merge vertex to the vertex to be removed in the second triangle to be removed
    fromHalfEdgeIndex = halfEdges.at(fromHalfEdgeIndex).getPreviousIndex(); // half-edge pointing to the merge vertex in the second triangle to be removed
    fromHalfEdgeIndex = halfEdges.at(fromHalfEdgeIndex).getPairIndex();     // half-edge pointing from the merge vertex paired to the second triangle to be removed

    toHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPairIndex();           // half-edge pointing from the merge vertex to the vertex to be removed in the second triangle to be removed
    toHalfEdgeIndex = halfEdges.at(toHalfEdgeIndex).getNextIndex();         // half-edge pointing from the vertex to be removed in the second triangle to be removed
    toHalfEdgeIndex = halfEdges.at(toHalfEdgeIndex).getPairIndex();         // half-edge pointing to the vertex to be removed paired to the second triangle to be removed

    halfEdges[fromHalfEdgeIndex].setPairIndex(toHalfEdgeIndex);
    halfEdges[toHalfEdgeIndex].setPairIndex(fromHalfEdgeIndex);

    // move merge vertex

    std::vector<double> const diff =
        wrapTo(toMergeIndex, fromMergeIndex, false);
    std::vector<double> position =
        vertices.at(toMergeIndex).getPosition();
    for (int dim=0; dim < 2; dim++) {
        position[dim] += diff[dim]/2;       // move vertex to half point between two merged vertices
    }
    vertices[toMergeIndex].setPosition(     // wrap position
        wrap(position));
    vertices[toMergeIndex].setUPosition(    // set unwrapped position back to wrapped position
        vertices.at(toMergeIndex).getPosition());

    // delete faces, junctions, half-edges, and vertex

    assert(halfEdgeIndex == halfEdges.at(halfEdgeIndex).getIndex());    // index of half-edge from deleted vertex
    pairHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPairIndex();     // index of half-edge to deleted vertex

    long int previousIndex, nextIndex;

    // first face
    previousIndex = halfEdges.at(halfEdgeIndex).getPreviousIndex();
    nextIndex = halfEdges.at(halfEdgeIndex).getNextIndex();
    faces.erase({halfEdgeIndex, previousIndex, nextIndex});             // delete face (first triangle)
    halfEdges.erase(previousIndex); halfEdges.erase(nextIndex);         // delete all half-edges from first triangle but the one from the erased junction

    // second face
    previousIndex = halfEdges.at(pairHalfEdgeIndex).getPreviousIndex();
    nextIndex = halfEdges.at(pairHalfEdgeIndex).getNextIndex();
    faces.erase({pairHalfEdgeIndex, previousIndex, nextIndex});         // delete face (second triangle)
    halfEdges.erase(previousIndex); halfEdges.erase(nextIndex);         // delete all half-edges from second triangle but the one from the erased junction

    junctions.erase({halfEdgeIndex, pairHalfEdgeIndex});                // delete junction
    halfEdges.erase(halfEdgeIndex); halfEdges.erase(pairHalfEdgeIndex); // delete half-edges in the erased junction

    sPVertices.erase(fromMergeIndex);                                   // delete self-propelled vertex
    vertices.erase(fromMergeIndex);                                     // delete vertex

    return fromMergeIndex;
}

long int const VertexModel::createJunction(
    long int const& halfEdgeIndex0, long int const& halfEdgeIndex1,
    double const& angle, double const& length) {

    // junctions cannot be split

    assert(!junctions.in(halfEdgeIndex0));
    assert(!junctions.in(halfEdgeIndex1));

    // create new vertex

    long int const vertexIndex = halfEdges.at(halfEdgeIndex0).getFromIndex();
    assert(vertexIndex == halfEdges.at(halfEdgeIndex1).getFromIndex()); // check that both half-edges go out of the same vertex
    vertices[vertexIndex].setHalfEdgeIndex(halfEdgeIndex0);             // re-assign identifying half-edge to `halfEdgeIndex0' which will remain attached to `vertexIndex'

    long int const newVertexIndex = maxKey(vertices) + 1;
    std::vector<double> const position =
        vertices.at(vertexIndex).getPosition();
    vertices.emplace(newVertexIndex, Vertex(
        // vertices is std::map so add with std::map::emplace
        newVertexIndex,     // vertexIndex
        position,           // position
        halfEdgeIndex1));   // halfEdgeIndex
    sPVertices[newVertexIndex] = {
        // sPVertices is MultiIntKeyDict<SPVertex> so add with initialiser-list (sent to SPVertex)
        newVertexIndex};    // vertexIndex

    // relabel origins and destinations of half of the half-edges

    long int halfEdgeIndex = halfEdgeIndex1;
    long int previousHalfEdgeIndex;
    while (halfEdgeIndex != halfEdgeIndex0) {

        assert(
            vertexIndex == halfEdges.at(halfEdgeIndex).getFromIndex());
        halfEdges[halfEdgeIndex].setFromIndex(newVertexIndex);

        previousHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPreviousIndex();
        assert(
            vertexIndex == halfEdges.at(previousHalfEdgeIndex).getToIndex());
        halfEdges[previousHalfEdgeIndex].setToIndex(newVertexIndex);

        halfEdgeIndex = halfEdges.at(previousHalfEdgeIndex).getPairIndex();
    }

    // create new half-edges, faces, and junctions

    Counter c(maxKey(halfEdges) + 1);

    long int const newHalfEdgeIndex0 = c();
    long int const previousNewHalfEdgeIndex0 = c();
    long int const nextNewHalfEdgeIndex0 = c();
    long int const newHalfEdgeIndex1 = c();
    long int const previousNewHalfEdgeIndex1 = c();
    long int const nextNewHalfEdgeIndex1 = c();

    // first face
    halfEdges.emplace(newHalfEdgeIndex0, HalfEdge(
        // halfEdges is std::map so add with std::map::emplace
        newHalfEdgeIndex0,                              // halfEdgeIndex
        halfEdges.at(halfEdgeIndex0).getToIndex(),      // fromIndex
        vertexIndex,                                    // toIndex
        previousNewHalfEdgeIndex0,                      // previousIndex
        nextNewHalfEdgeIndex0,                          // nextIndex
        halfEdgeIndex0));                               // pairIndex
    halfEdges.emplace(previousNewHalfEdgeIndex0, HalfEdge(
        previousNewHalfEdgeIndex0,                      // halfEdgeIndex
        newVertexIndex,                                 // fromIndex
        halfEdges.at(halfEdgeIndex0).getToIndex(),      // toIndex
        nextNewHalfEdgeIndex0,                          // previousIndex
        newHalfEdgeIndex0,                              // nextIndex
        halfEdges.at(halfEdgeIndex0).getPairIndex()));  // pairIndex
    halfEdges.emplace(nextNewHalfEdgeIndex0, HalfEdge(
        nextNewHalfEdgeIndex0,                          // halfEdgeIndex
        vertexIndex,                                    // fromIndex
        newVertexIndex,                                 // toIndex
        newHalfEdgeIndex0,                              // previousIndex
        previousNewHalfEdgeIndex0,                      // nextIndex
        nextNewHalfEdgeIndex1));                        // pairIndex
    faces
        [{newHalfEdgeIndex0, previousNewHalfEdgeIndex0, nextNewHalfEdgeIndex0}]
        // faces is MultiIntKeyDict<Face> so add with initialiser-list (sent to Face)
        = {newHalfEdgeIndex0};                          // create new face
    halfEdges[halfEdges.at(halfEdgeIndex0).getPairIndex()].setPairIndex(
        previousNewHalfEdgeIndex0);                     // relabel pair of old pair of halfEdgeIndex0
    halfEdges[halfEdgeIndex0].setPairIndex(
        newHalfEdgeIndex0);                             // relabel pair of halfEdgeIndex0

    // second face
    halfEdges.emplace(newHalfEdgeIndex1, HalfEdge(
        newHalfEdgeIndex1,                              // halfEdgeIndex
        halfEdges.at(halfEdgeIndex1).getToIndex(),      // fromIndex
        newVertexIndex,                                 // toIndex
        previousNewHalfEdgeIndex1,                      // previousIndex
        nextNewHalfEdgeIndex1,                          // nextIndex
        halfEdgeIndex1));                               // pairIndex
    halfEdges.emplace(previousNewHalfEdgeIndex1, HalfEdge(
        previousNewHalfEdgeIndex1,                      // halfEdgeIndex
        vertexIndex,                                    // fromIndex
        halfEdges.at(halfEdgeIndex1).getToIndex(),      // toIndex
        nextNewHalfEdgeIndex1,                          // previousIndex
        newHalfEdgeIndex1,                              // nextIndex
        halfEdges.at(halfEdgeIndex1).getPairIndex()));  // pairIndex
    halfEdges.emplace(nextNewHalfEdgeIndex1, HalfEdge(
        nextNewHalfEdgeIndex1,                          // halfEdgeIndex
        newVertexIndex,                                 // fromIndex
        vertexIndex,                                    // toIndex
        newHalfEdgeIndex1,                              // previousIndex
        previousNewHalfEdgeIndex1,                      // nextIndex
        nextNewHalfEdgeIndex0));                        // pairIndex
    faces
        [{newHalfEdgeIndex1, previousNewHalfEdgeIndex1, nextNewHalfEdgeIndex1}]
        = {newHalfEdgeIndex1};                          // create new face
    halfEdges[halfEdges.at(halfEdgeIndex1).getPairIndex()].setPairIndex(
        previousNewHalfEdgeIndex1);                     // relabel pair of old pair of halfEdgeIndex1
    halfEdges[halfEdgeIndex1].setPairIndex(
        newHalfEdgeIndex1);                             // relabel pair of halfEdgeIndex1

    // create junction

    junctions[{nextNewHalfEdgeIndex0, nextNewHalfEdgeIndex1}] = {
        // faces is MultiIntKeyDict<Junction> so add with initialiser-list (sent to Junction)
        nextNewHalfEdgeIndex0};    // halfEdgeIndex

    // move vertices apart

    std::vector<double> direction = {std::cos(angle), std::sin(angle)};
    if (cross2(getHalfEdgeVector(newHalfEdgeIndex0), direction) < 0) {
        for (int dim=0; dim < 2; dim++) { direction[dim] *= -1; }   // enforces that the newly created face has anticlockwise orientation
    }

    std::vector<double> newPosition(position), newPositionBis(position);
    for (int dim=0; dim < 2; dim++) {
        newPosition[dim] -= direction[dim]*length/2.;
        newPositionBis[dim] += direction[dim]*length/2.;
    }
    vertices[vertexIndex].setPosition(      // wrap position
        wrap(newPosition));
    vertices[vertexIndex].setUPosition(     // set unwrapped position back to wrapped position
        vertices.at(vertexIndex).getPosition());
    vertices[newVertexIndex].setPosition(   // wrap position
        wrap(newPositionBis));
    vertices[newVertexIndex].setUPosition(  // set unwrapped position back to wrapped position
        vertices.at(newVertexIndex).getPosition());

    return newVertexIndex;
}

void VertexModel::initRegularTriangularLattice(
    long int const& size, double const& junctionLength) {

    assert(size%6 == 0);

    reset();                                                // reset vertices, junctions, and cells
    systemSize[0] = size*junctionLength;                    // length of the periodic box in the x-direction
    systemSize[1] = size*junctionLength*std::sqrt(3.)/2.;   // length of the periodic box in the y-direction

    auto getIndex = [&size](long int const& line, long int const& column) {
        int const col = pmod((column + qpmod(line, size)*(size/2)), size);  // the system is periodic along an inclined vertical axis
        int const lin = pmod(line, size);                                   // the system is periodic along the horizontal axis
        return (long int const) lin*size + col;
    };

    // loop over vertices
    long int vertexIndex;
    std::vector<double> position(2, 0);
    long int A, B, C, D, E, F, G;
    for (int line=0; line < size; line++) {
        for (int column=0; column < size; column++) {
            vertexIndex = line*size + column;

            // create vertex
            position[0] = junctionLength*(0.25 + line*0.5 + column);
            position[1] = junctionLength*(0.5 + line)*std::sqrt(3.)/2.;
            vertices.emplace(vertexIndex, Vertex(
                // vertices is std::map so add with std::map::emplace
                vertexIndex,
                wrap(position),     // wrapped position of the vertex on the regular triangular lattice
                6*vertexIndex + 0,  // index of a single half-edge going out of this vertex
                false));            // there is no open boundary

            // create cell or self-propelled vertex
            if ((line - column)%3 == 0) {                           // condition for vertex to be a cell centre
                cells[vertexIndex] = {
                    // cells is MultiIntKeyDict<Cell> so add with initialiser-list (sent to Cell)
                    vertexIndex,                            // vertexIndex
                    (junctionLength*junctionLength)         // A0
                        *(3./2.)/std::tan(std::numbers::pi/6.),
                    p0};                                    // p0
            }
            else {                          // ... or a self-propelled vertex
                sPVertices[vertexIndex] = {
                    // sPVertices is MultiIntKeyDict<SPVertex> so add with initialiser-list (sent to SPVertex)
                    vertexIndex,                            // vertexIndex
                    2*std::numbers::pi*random.random01(),   // theta
                    v0,                                     // v0
                    Dr};                                    // Dr
            }

            // vertex indices for half-edge construction
            A = getIndex(line, column);
            assert(A == vertexIndex);
            B = getIndex(line, column + 1);
            C = getIndex(line + 1, column);
            D = getIndex(line + 1, column - 1);
            E = getIndex(line, column - 1);
            F = getIndex(line - 1, column);
            G = getIndex(line - 1, column + 1);

            // create half-edges
            // (0) from vertex (A) to right (B)
            halfEdges.emplace(6*A + 0, HalfEdge(
                // vertices is std::map so add with std::map::emplace
                6*A + 0,    // halfEdgeIndex
                A,          // fromIndex
                B,          // toIndex
                6*A + 4,    // previousIndex
                6*A + 2,    // nextIndex
                6*A + 1));  // pairIndex
            // (1) from right (B) to vertex (A)
            halfEdges.emplace(6*A + 1, HalfEdge(
                6*A + 1,    // halfEdgeIndex
                B,          // fromIndex
                A,          // toIndex
                6*G + 5,    // previousIndex
                6*F + 3,    // nextIndex
                6*A + 0));  // pairIndex
            // (2) from right (B) to top (C)
            halfEdges.emplace(6*A + 2, HalfEdge(
                6*A + 2,    // halfEdgeIndex
                B,          // fromIndex
                C,          // toIndex
                6*A + 0,    // previousIndex
                6*A + 4,    // nextIndex
                6*A + 3));  // pairIndex
            // (3) from top (C) to right (B)
            halfEdges.emplace(6*A + 3, HalfEdge(
                6*A + 3,    // halfEdgeIndex
                C,          // fromIndex
                B,          // toIndex
                6*C + 1,    // previousIndex
                6*B + 5,    // nextIndex
                6*A + 2));  // pairIndex
            // (4) from top (C) to vertex (A)
            halfEdges.emplace(6*A + 4, HalfEdge(
                6*A + 4,    // halfEdgeIndex
                C,          // fromIndex
                A,          // toIndex
                6*A + 2,    // previousIndex
                6*A + 0,    // nextIndex
                6*A + 5));  // pairIndex
            // (5) from vertex (A) to top (C)
            halfEdges.emplace(6*A + 5, HalfEdge(
                6*A + 5,    // halfEdgeIndex
                A,          // fromIndex
                C,          // toIndex
                6*E + 3,    // previousIndex
                6*D + 1,    // nextIndex
                6*A + 4));  // pairIndex

            // create faces
            // (0) right above vertex (A)
            faces[{6*A + 0, 6*A + 2, 6*A + 4}] = {
                // faces is MultiIntKeyDict<Face> so add with initialiser-list (sent to Face)
                6*A + 0};   // halfEdgeIndex
            faces[{6*A + 1, 6*F + 3, 6*G + 5}] = {
                6*A + 1};   // halfEdgeIndex

            // create junctions between vertices which are not cell centres
            if (((A/size - A%size)%3 != 0) && ((B/size - B%size)%3 != 0)) {
                junctions[{6*A + 0, 6*A + 1}] = {
                    // junctions is MultiIntKeyDict<Junction> so add with initialiser-list (sent to Junction)
                    6*A + 0};   // halfEdgeIndex
            }
            if (((B/size - B%size)%3 != 0) && ((C/size - C%size)%3 != 0)) {
                junctions[{6*A + 2, 6*A + 3}] = {
                    6*A + 2};   // halfEdgeIndex
            }
            if (((C/size - C%size)%3 != 0) && ((A/size - A%size)%3 != 0)) {
                junctions[{6*A + 4, 6*A + 5}] = {
                    6*A + 4};   // halfEdgeIndex
            }
        }
    }

    checkMesh();    // check mesh construction
}

void VertexModel::initOpenRegularTriangularLattice(
    long int const& size, double const& junctionLength) {
    // (TESTING) replace a single cell with a hole

    assert(getBoundary());  // mesh must be set to check for boundaries
    assert(size%6 == 0);

    reset();                                                // reset vertices, junctions, and cells
    systemSize[0] = size*junctionLength;                    // length of the periodic box in the x-direction
    systemSize[1] = size*junctionLength*std::sqrt(3.)/2.;   // length of the periodic box in the y-direction

    auto getIndex = [&size](long int const& line, long int const& column) {
        int const col = pmod((column + qpmod(line, size)*(size/2)), size);  // the system is periodic along an inclined vertical axis
        int const lin = pmod(line, size);                                   // the system is periodic along the horizontal axis
        return (long int const) lin*size + col;
    };

    // loop over vertices
    long int vertexIndex;
    std::vector<double> position(2, 0);
    long int A, B, C, D, E, F, G;
    for (int line=0; line < size; line++) {
        for (int column=0; column < size; column++) {
            vertexIndex = line*size + column;

            // create vertex
            position[0] = junctionLength*(0.25 + line*0.5 + column);
            position[1] = junctionLength*(0.5 + line)*std::sqrt(3.)/2.;
            vertices.emplace(vertexIndex, Vertex(
                // vertices is std::map so add with std::map::emplace
                vertexIndex,
                wrap(position),                                 // wrapped position of the vertex on the regular triangular lattice
                6*vertexIndex + 0,                              // index of a single half-edge going out of this vertex
                (line == size/2 - 1 && column == size/2 - 1))); // (TESTING) create 1 boundary vertex (here a hole)

            // create cell or self-propelled vertex
            if ((line - column)%3 == 0                              // condition for vertex to be a cell centre
                && !(line == size/2 - 1 && column == size/2 - 1)) { // ignore created boundary vertex
                cells[vertexIndex] = {
                    // cells is MultiIntKeyDict<Cell> so add with initialiser-list (sent to Cell)
                    vertexIndex,                                    // vertexIndex
                    (junctionLength*junctionLength)                 // A0
                        *(3./2.)/std::tan(std::numbers::pi/6.),
                    p0};                                            // p0
            }
            else {                          // ... or a self-propelled vertex
                sPVertices[vertexIndex] = {
                    // sPVertices is MultiIntKeyDict<SPVertex> so add with initialiser-list (sent to SPVertex)
                    vertexIndex,                            // vertexIndex
                    2*std::numbers::pi*random.random01(),   // theta
                    v0,                                     // v0
                    Dr};                                    // Dr
            }

            // vertex indices for half-edge construction
            A = getIndex(line, column);
            assert(A == vertexIndex);
            B = getIndex(line, column + 1);
            C = getIndex(line + 1, column);
            D = getIndex(line + 1, column - 1);
            E = getIndex(line, column - 1);
            F = getIndex(line - 1, column);
            G = getIndex(line - 1, column + 1);

            // create half-edges
            // (0) from vertex (A) to right (B)
            halfEdges.emplace(6*A + 0, HalfEdge(
                // vertices is std::map so add with std::map::emplace
                6*A + 0,    // halfEdgeIndex
                A,          // fromIndex
                B,          // toIndex
                6*A + 4,    // previousIndex
                6*A + 2,    // nextIndex
                6*A + 1));  // pairIndex
            // (1) from right (B) to vertex (A)
            halfEdges.emplace(6*A + 1, HalfEdge(
                6*A + 1,    // halfEdgeIndex
                B,          // fromIndex
                A,          // toIndex
                6*G + 5,    // previousIndex
                6*F + 3,    // nextIndex
                6*A + 0));  // pairIndex
            // (2) from right (B) to top (C)
            halfEdges.emplace(6*A + 2, HalfEdge(
                6*A + 2,    // halfEdgeIndex
                B,          // fromIndex
                C,          // toIndex
                6*A + 0,    // previousIndex
                6*A + 4,    // nextIndex
                6*A + 3));  // pairIndex
            // (3) from top (C) to right (B)
            halfEdges.emplace(6*A + 3, HalfEdge(
                6*A + 3,    // halfEdgeIndex
                C,          // fromIndex
                B,          // toIndex
                6*C + 1,    // previousIndex
                6*B + 5,    // nextIndex
                6*A + 2));  // pairIndex
            // (4) from top (C) to vertex (A)
            halfEdges.emplace(6*A + 4, HalfEdge(
                6*A + 4,    // halfEdgeIndex
                C,          // fromIndex
                A,          // toIndex
                6*A + 2,    // previousIndex
                6*A + 0,    // nextIndex
                6*A + 5));  // pairIndex
            // (5) from vertex (A) to top (C)
            halfEdges.emplace(6*A + 5, HalfEdge(
                6*A + 5,    // halfEdgeIndex
                A,          // fromIndex
                C,          // toIndex
                6*E + 3,    // previousIndex
                6*D + 1,    // nextIndex
                6*A + 4));  // pairIndex

            // create faces
            // (0) right above vertex (A)
            faces[{6*A + 0, 6*A + 2, 6*A + 4}] = {
                // faces is MultiIntKeyDict<Face> so add with initialiser-list (sent to Face)
                6*A + 0};   // halfEdgeIndex
            faces[{6*A + 1, 6*F + 3, 6*G + 5}] = {
                6*A + 1};   // halfEdgeIndex

            // create junctions between vertices which are not cell centres
            if (((A/size - A%size)%3 != 0) && ((B/size - B%size)%3 != 0)) {
                junctions[{6*A + 0, 6*A + 1}] = {
                    // junctions is MultiIntKeyDict<Junction> so add with initialiser-list (sent to Junction)
                    6*A + 0};   // halfEdgeIndex
            }
            if (((B/size - B%size)%3 != 0) && ((C/size - C%size)%3 != 0)) {
                junctions[{6*A + 2, 6*A + 3}] = {
                    6*A + 2};   // halfEdgeIndex
            }
            if (((C/size - C%size)%3 != 0) && ((A/size - A%size)%3 != 0)) {
                junctions[{6*A + 4, 6*A + 5}] = {
                    6*A + 4};   // halfEdgeIndex
            }
        }
    }

    checkMesh();    // check mesh construction
}

void VertexModel::initOpenRegularHexagonalLattice(
    long int const& nCells, double const& junctionLength) {

    long int const nnCells = sqrt(nCells);
    assert(nnCells*nnCells == nCells);

    reset();                                    // reset vertices, junctions, and cells
    systemSize[0] = 2*nnCells*junctionLength;   // length of the periodic box in the x-direction
    systemSize[1] = 2*nnCells*junctionLength;   // length of the periodic box in the y-direction

    // first loop over cells to create vertices and self-propelled vertices
    double const x0 = (nnCells/2. + 0.5)*junctionLength;                        // x-position of first cell centre
    double const y0 = (nnCells - (nnCells + 1)*sqrt(3)/4 + 0.5)*junctionLength; // y-position of first cell centre
    auto vertexIndex = [&nnCells]
        (long int const& col, long int const& line, long int const& k)
            { return 7*(col + line*nnCells) + k; };
    std::vector<double> position(2, 0);
    std::map<std::vector<double>, long int> vertexPosMap;
    std::map<long int, long int> vertexIndexMap;
    std::map<long int, long int> exteriorVertexNext;
    for (long int col=0; col < nnCells; col++) {
        for (long int line=0; line < nnCells; line++) {

            // create cell centre vertex
            position[0] = x0 + (col + (col%2)/2.)*junctionLength;
            position[1] = y0 + line*(1. + sqrt(3)/2.)*junctionLength;
            vertices.emplace(vertexIndex(col, line, 0), Vertex(
                // vertices is std::map so add with std::map::emplace
                vertexIndex(col, line, 0),
                wrap(position)));   // wrapped position of the vertex
            cells[vertexIndex(col, line, 0)] = {
                // cells is MultiIntKeyDict<Cell> so add with initialiser-list (sent to Cell)
                vertexIndex(col, line, 0),
                (junctionLength*junctionLength) // A0
                    *(3./2.)/std::tan(std::numbers::pi/6.),
                p0};                            // p0

            // create cell corner vertices
            for (long int k=1; k <= 6; k++) {
                position[0] =
                    vertices[vertexIndex(col, line, 0)].getPosition()[0]
                        + cos(-std::numbers::pi/2. + k*std::numbers::pi/6.);
                position[1] =
                    vertices[vertexIndex(col, line, 0)].getPosition()[1]
                        + sin(-std::numbers::pi/2. + k*std::numbers::pi/6.);

                if (vertexPosMap.find(wrap(position)) != vertexPosMap.end()) {  // vertex already exists
                    // relations between vertices for initialisation
                    vertexIndexMap[vertexIndex(col, line, k)] =
                        vertexPosMap.find(position)->second;
                }
                else {                                                          // vertex does not exist
                    // create vertex
                    vertices.emplace(vertexIndex(col, line, k), Vertex(
                        // vertices is std::map so add with std::map::emplace
                        vertexIndex(col, line, k),
                        wrap(position)));   // wrapped position of the vertex
                    // crate self-propelled vertex
                    sPVertices[vertexIndex(col, line, k)] = {
                        // sPVertices is MultiIntKeyDict<SPVertex> so add with initialiser-list (sent to SPVertex)
                        vertexIndex(col, line, k),              // vertexIndex
                        2*std::numbers::pi*random.random01(),   // theta
                        v0,                                     // v0
                        Dr};                                    // Dr
                    // relations between vertices for initialisation
                    vertexIndexMap[vertexIndex(col, line, k)] =
                        vertexIndex(col, line, k);
                    vertexPosMap[wrap(position)] =
                        vertexIndex(col, line, k);
                    // relations between exterior vertices for initialisation
                    if (exteriorVertexNext.find(vertexIndex(col, line, k))
                        == exteriorVertexNext.end()) {  // next vertex of exterior vertex not computed
                        if (line == 0) {                // vertices on the bottom line
                            if (k == 1) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 0);
                            }
                            else if (k == 0) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 5);
                            }
                        }
                        if (line == nnCells - 1) {      // vertices on the top line
                            if (k == 4) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 3);
                            }
                            else if (k == 3) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 2);
                            }
                        }
                        if (col == 0) {                 // vertices on the left column
                            if (k == 5) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 4);
                            }
                            else if (k == 4 && line%2 == 0) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 3);
                            }
                        }
                        if (col == nnCells - 1) {       // vertices on the right column
                            if (k == 2) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 1);
                            }
                            else if (k == 1 && line%2 == 1) {
                                exteriorVertexNext[vertexIndex(col, line, k)] =
                                    vertexIndex(col, line, 0);
                            }
                        }
                    }
                }
            }
        }
    }

    // boundary vertex
    long int const boundaryIndex = vertexIndex(nnCells - 1, nnCells - 1, 7);
    vertices.emplace(boundaryIndex, Vertex(
        // vertices is std::map so add with std::map::emplace
        boundaryIndex,
        std::vector<double>(2, 0),
        -1,
        true));

    // second loop over cells to create interior half-edges and faces
    long int fromIndex, toIndex;
    Counter c;
    long int triangle[3];
    std::map<std::tuple<long int, long int>, long int> halfEdgeIndexMap;
    for (long int col=0; col < nnCells; col++) {
        for (long int line=0; line < nnCells; line++) {

            // create inner half-edges
            for (long int k=1; k <= 6; k++) {
                triangle[0] = c(); triangle[1] = c(); triangle[2] = c();
                // half-edge from cell centre to first corner
                fromIndex = vertexIndexMap[vertexIndex(col, line, 0)];
                toIndex = vertexIndexMap[vertexIndex(col, line, k)];
                if (
                    halfEdgeIndexMap.find(std::make_tuple(fromIndex, toIndex))
                    == halfEdgeIndexMap.end()) {    // half-edge does not exist
                    halfEdges.emplace(triangle[0], HalfEdge(
                        // halfEdges is std::map so add with std::map::emplace
                        triangle[0],    // halfEdgeIndex
                        fromIndex,      // fromIndex
                        toIndex,        // toIndex
                        triangle[2],    // previousIndex
                        triangle[1]));  // nextIndex
                    halfEdgeIndexMap[std::make_tuple(fromIndex, toIndex)] =
                        triangle[0];
                    vertices[fromIndex].setHalfEdgeIndex(triangle[0]);  // this will be changed multiple times so is redundant
                }
                // half-edge from first corner to second corner
                fromIndex = vertexIndexMap[vertexIndex(col, line, k)];
                toIndex = vertexIndexMap[vertexIndex(col, line, k + 1)];
                if (
                    halfEdgeIndexMap.find(std::make_tuple(fromIndex, toIndex))
                    == halfEdgeIndexMap.end()) {    // half-edge does not exist
                    halfEdges.emplace(triangle[1], HalfEdge(
                        // halfEdges is std::map so add with std::map::emplace
                        triangle[1],    // halfEdgeIndex
                        fromIndex,      // fromIndex
                        toIndex,        // toIndex
                        triangle[0],    // previousIndex
                        triangle[2]));  // nextIndex
                    halfEdgeIndexMap[std::make_tuple(fromIndex, toIndex)] =
                        triangle[1];
                    vertices[fromIndex].setHalfEdgeIndex(triangle[1]);  // this will be changed multiple times so is redundant
                }
                // half-edge from second corner to cell centre
                fromIndex = vertexIndexMap[vertexIndex(col, line, k + 1)];
                toIndex = vertexIndexMap[vertexIndex(col, line, 0)];
                if (
                    halfEdgeIndexMap.find(std::make_tuple(fromIndex, toIndex))
                    == halfEdgeIndexMap.end()) {    // half-edge does not exist
                    halfEdges.emplace(triangle[2], HalfEdge(
                        // halfEdges is std::map so add with std::map::emplace
                        triangle[2],    // halfEdgeIndex
                        fromIndex,      // fromIndex
                        toIndex,        // toIndex
                        triangle[1],    // previousIndex
                        triangle[0]));  // nextIndex
                    halfEdgeIndexMap[std::make_tuple(fromIndex, toIndex)] =
                        triangle[2];
                    vertices[fromIndex].setHalfEdgeIndex(triangle[2]);  // this will be changed multiple times so is redundant
                }
                // face
                faces[{triangle[0], triangle[1], triangle[2]}] = {
                    // faces is MultiIntKeyDict<Face> so add with initialiser-list (sent to Face)
                    triangle[0]};   // halfEdgeIndex
            }
        }
    }

    // create external half-edges
    std::map<long int, long int>::iterator it;
    for (it=exteriorVertexNext.begin(); it != exteriorVertexNext.end(); ++it) {
        triangle[0] = c(); triangle[1] = c(); triangle[2] = c();
        // half-edge from first exterior vertex to second exterior vertex
        fromIndex = vertexIndexMap[it->first];
        toIndex = vertexIndexMap[it->second];
        assert(
            halfEdgeIndexMap.find(std::make_tuple(fromIndex, toIndex))
            == halfEdgeIndexMap.end()); // half-edge should not exist
        halfEdges.emplace(triangle[0], HalfEdge(
            // halfEdges is std::map so add with std::map::emplace
            triangle[0],    // halfEdgeIndex
            fromIndex,      // fromIndex
            toIndex,        // toIndex
            triangle[2],    // previousIndex
            triangle[1]));  // nextIndex
            halfEdgeIndexMap[std::make_tuple(fromIndex, toIndex)] =
                triangle[0];
        vertices[fromIndex].setHalfEdgeIndex(triangle[0]);  // this will be changed multiple times so is redundant
        // half-edge from second exterior vertex to boundary vertex
        fromIndex = vertexIndexMap[it->second];
        toIndex = boundaryIndex;
        assert(
            halfEdgeIndexMap.find(std::make_tuple(fromIndex, toIndex))
            == halfEdgeIndexMap.end()); // half-edge should not exist
        halfEdges.emplace(triangle[1], HalfEdge(
            // halfEdges is std::map so add with std::map::emplace
            triangle[1],    // halfEdgeIndex
            fromIndex,      // fromIndex
            toIndex,        // toIndex
            triangle[0],    // previousIndex
            triangle[2]));  // nextIndex
            halfEdgeIndexMap[std::make_tuple(fromIndex, toIndex)] =
                triangle[1];
        vertices[fromIndex].setHalfEdgeIndex(triangle[1]);  // this will be changed multiple times so is redundant
        // half-edge from boundary vertex to first boundary vertex
        fromIndex = boundaryIndex;
        toIndex = vertexIndexMap[it->first];
        assert(
            halfEdgeIndexMap.find(std::make_tuple(fromIndex, toIndex))
            == halfEdgeIndexMap.end()); // half-edge should not exist
        halfEdges.emplace(triangle[2], HalfEdge(
            // halfEdges is std::map so add with std::map::emplace
            triangle[2],    // halfEdgeIndex
            fromIndex,      // fromIndex
            toIndex,        // toIndex
            triangle[1],    // previousIndex
            triangle[0]));  // nextIndex
            halfEdgeIndexMap[std::make_tuple(fromIndex, toIndex)] =
                triangle[2];
        vertices[fromIndex].setHalfEdgeIndex(triangle[2]);  // this will be changed multiple times so is redundant
    }

    // associate pair of half-edges
    long int halfEdgeIndex, pairHalfEdgeIndex;
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
        fromIndex = (it->second).getFromIndex();
        toIndex = (it->second).getToIndex();
        halfEdgeIndex = (it->second).getIndex();
        assert(halfEdgeIndexMap.find(std::make_tuple(toIndex, fromIndex))
            != halfEdgeIndexMap.end()); // pair half-edge should exist
        pairHalfEdgeIndex =
            halfEdgeIndexMap.find(std::make_tuple(toIndex, fromIndex))->second;
        halfEdges[halfEdgeIndex].setPairIndex(pairHalfEdgeIndex);
        junctions[{halfEdgeIndex, pairHalfEdgeIndex}] = {
            // junctions is MultiIntKeyDict<Junction> so add with initialiser-list (sent to Junction)
            halfEdgeIndex};   // halfEdgeIndex
    }

    checkMesh();    // check mesh construction
}

