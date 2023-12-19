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

