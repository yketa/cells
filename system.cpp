#include <cmath>
#include <numbers>
#include <assert.h>

#include "system.hpp"
#include "tools.hpp"

void VertexModel::integrate(double const& dt,
    double const& delta, double const& epsilon) {

    // get forces

    std::map<long int, std::vector<double>> const forces = getForces();

    long int vertexIndex;

    // integrate

    for (auto it=forces.begin(); it != forces.end(); ++it) {
        vertexIndex = it->first;
        for (int dim=0; dim < 2; dim++) {
            vertices[vertexIndex].getPosition()[dim] +=
                forces.at(vertexIndex)[dim]*dt;     // Euler integration
        }
        wrap(vertices[vertexIndex].getPosition());  // wrap with respect to boundary conditions
    }
    for (SPVertex* sPVertex : sPVertices.getValues()) {
        *(sPVertex->gettheta()) +=
            std::sqrt(2*sPVertex->getDr()*dt)*random.gauss();
    }

    // move cell centres

    double* cellPosition; std::vector<double> initialCellPosition(0);
    std::vector<long int> neighbourVerticesIndices(0); int numberNeighbours;
    std::vector<double> disp(0);
    for (Cell* cell : cells.getValues()) {
        vertexIndex = cell->getVertexIndex();

        cellPosition = vertices[vertexIndex].getPosition();
        initialCellPosition = {cellPosition[0], cellPosition[1]};

        neighbourVerticesIndices = getNeighbourVertices(vertexIndex)[0];
        numberNeighbours = neighbourVerticesIndices.size();
        for (long int neighbourVertexIndex : neighbourVerticesIndices) {
            disp = wrapDiff(
                &(initialCellPosition[0]),
                vertices[neighbourVertexIndex].getPosition());

            for (int dim=0; dim < 2; dim++) {
                cellPosition[dim] += disp[dim]/numberNeighbours;
            }
        }
    }

    // perform T1s

    doT1(delta, epsilon);

    // update time

    time += dt;
}

std::map<long int,std::vector<double>> const VertexModel::getForces() {

    std::map<long int,std::vector<double>> forces;
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        forces[it->first] = {0, 0}; // initialise forces at each vertex
    }

    long int cellVertexIndex;
    std::vector<long int> neighbourVerticesIndices(0); int numberNeighbours;
    long int neighbourVertexIndex;
    long int previousNeighbourVertexIndex, nextNeighbourVertexIndex;
    std::vector<double> fromTo(0), fromToBis(0);
    double areaTerm;
    for (Cell* cell : cells.getValues()) {          // loop over cells
        cellVertexIndex = cell->getVertexIndex();

        *cell->getArea() =      // update area
            getVertexToNeighboursArea(cellVertexIndex);
        *cell->getPerimeter() = // update perimeter
            getVertexToNeighboursPerimeter(cellVertexIndex);

        neighbourVerticesIndices = getNeighbourVertices(cellVertexIndex)[0];
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
            fromToBis = // vector from cell to next vertex
                wrapTo(cellVertexIndex, nextNeighbourVertexIndex,
                    false);
            for (int dim=0; dim < 2; dim++) {

                areaTerm =
                    -cell->getkA()*(*cell->getArea() - cell->getA0())*(1./2.)*(
                        cross2z(fromTo)[dim] - cross2z(fromToBis)[dim]);

                forces[neighbourVertexIndex][dim] += areaTerm;
                forces[cellVertexIndex][dim] -= areaTerm;   // enforce Newton's third law w/ cell centrs
            }

            // perimeter term

            fromTo =    // unit vector from next vertex to vertex
                wrapTo(nextNeighbourVertexIndex, neighbourVertexIndex,
                    true);
            fromToBis = // unit vector from previous vertex to vertex
                wrapTo(previousNeighbourVertexIndex, neighbourVertexIndex,
                    true);
            for (int dim=0; dim < 2; dim++) {

                forces[neighbourVertexIndex][dim] +=
                    -cell->getkP()*(*cell->getPerimeter() - cell->getP0())*(
                        fromTo[dim] + fromToBis[dim]);
            }
        }
    }

    // self-propelled vertices

    long int vertexIndex;
    for (SPVertex* sPVertex : sPVertices.getValues()) {
        vertexIndex = sPVertex->getVertexIndex();

        for (int dim=0; dim < 2; dim++) {

            forces[vertexIndex][dim] +=
                sPVertex->getv0()
                    *std::cos(*sPVertex->gettheta() - dim*std::numbers::pi/2);
        }
    }

    return forces;
}

void VertexModel::doT1(double const& delta, double const& epsilon) {

    // identify small junctions

    std::vector<long int> halfEdgeIndices(0);
    for (Junction* junction : junctions.getValues()) {
        if (getEdgeLength(junction->getHalfEdgeIndex()) < delta) {
            halfEdgeIndices.push_back(junction->getHalfEdgeIndex());
        }
    }
    random.shuffle(halfEdgeIndices);    // do T1s in a random order

    // perform T1s

    long int fromMergeIndex, toMergeIndex;
    std::vector<std::vector<long int>> neighbours;
    std::vector<long int> neighboursFromMerge, halfEdgesNeighboursFromMerge;
    std::vector<long int> neighboursToMerge, halfEdgesNeighboursToMerge;
    int numberNeighbours;
    std::vector<long int> halfEdgeToCellsIndices;
    long int vertexIndex;
    long int createHalfEdgeIndex0, createHalfEdgeIndex1;
    double angle;
    for (long int mergeHalfEdgeIndex : halfEdgeIndices) {
        std::cout << "merge: " << mergeHalfEdgeIndex << std::endl;

        // identify half-edge to split to create new junction

        fromMergeIndex = *halfEdges[mergeHalfEdgeIndex].getFromIndex(); // (first) vertex to be merge into neighbour
        toMergeIndex = *halfEdges[mergeHalfEdgeIndex].getToIndex();     // (second) vertex towards which neighbour is merged

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
            if (cells.in(vertexIndex)) {                                // cell centres neighbours of the first vertex...
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
            if (cells.in(vertexIndex)) {                                // cell centres neighbours of the second vertex...
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

        // merge vertices

        mergeVertices(mergeHalfEdgeIndex);

        // create new vertex

        createJunction(createHalfEdgeIndex0, createHalfEdgeIndex1,
            angle, delta + epsilon);
    }

    if (halfEdgeIndices.size() > 0) { checkMesh(); }
}

long int const VertexModel::mergeVertices(long int const& halfEdgeIndex) {

    // identify vertices to merge

    long int const fromMergeIndex = // (first) vertex to be merge into neighbour
        *halfEdges[halfEdgeIndex].getFromIndex();
    long int const toMergeIndex =   // (second) vertex towards which neighbour is merged
        *halfEdges[halfEdgeIndex].getToIndex();

    // relabel half-edges origins and destinations

    long int previousHalfEdgeIndex =            // half-edge pointing to the vertex to be removed in the first triangle to be removed
        *halfEdges[halfEdgeIndex].getPreviousIndex();
    long int const endPreviousHalfEdgeIndex =   // half-edge pointing to the vertex to be removed in the second triangle to be removed
        *halfEdges[halfEdgeIndex].getPairIndex();

    long int pairHalfEdgeIndex;
    while (true) {

        pairHalfEdgeIndex = *halfEdges[previousHalfEdgeIndex].getPairIndex();   // half-edge pointing from the vertex to be removed and to be relabelled

        // half-edge coming from vertex to be removed
        assert
            (fromMergeIndex == *halfEdges[pairHalfEdgeIndex].getFromIndex());
        *halfEdges[pairHalfEdgeIndex].getFromIndex() = toMergeIndex;

        // half-edge going to vertex to be removed
        assert
            (fromMergeIndex == *halfEdges[previousHalfEdgeIndex].getToIndex());
        *halfEdges[previousHalfEdgeIndex].getToIndex() = toMergeIndex;

        if (previousHalfEdgeIndex == endPreviousHalfEdgeIndex) {    // all triangles which the vertex to be removed belongs to have been explored and half-edges origins and destinations relabelled
            assert(pairHalfEdgeIndex == halfEdgeIndex);
            break;
        }

        previousHalfEdgeIndex = // half-edge pointing to the vertex to be removed and to be relabelled
            *halfEdges[pairHalfEdgeIndex].getPreviousIndex();
    }

    // reassign identifying half-edges of points belonging to triangles to be removed

    long int nextHalfEdgeIndex, thirdVertexIndex;

    previousHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPreviousIndex();   // half-edge pointing to the vertex to be removed in the first triangle to be removed
    pairHalfEdgeIndex = *halfEdges[previousHalfEdgeIndex].getPairIndex();   // half-edge pointing from the vertex to be removed

    assert(toMergeIndex == *halfEdges[pairHalfEdgeIndex].getFromIndex());
    *vertices[toMergeIndex].getHalfEdgeIndex() = pairHalfEdgeIndex;

    nextHalfEdgeIndex = *halfEdges[halfEdgeIndex].getNextIndex();       // half-edge pointing from the merge vertex in the first triangle to be removed
    pairHalfEdgeIndex = *halfEdges[nextHalfEdgeIndex].getPairIndex();   // half-edge pointing to the merge vertex
    thirdVertexIndex = *halfEdges[pairHalfEdgeIndex].getFromIndex();    // index of third vertex in the first triangle to be removed

    *vertices[thirdVertexIndex].getHalfEdgeIndex() = pairHalfEdgeIndex;

    pairHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPairIndex();       // half-edge pointing to the vertex to be removed in the second triangle to be removed
    nextHalfEdgeIndex = *halfEdges[pairHalfEdgeIndex].getNextIndex();   // half-edge pointing from the vertex to be removed in the second triangle to be removed
    pairHalfEdgeIndex = *halfEdges[nextHalfEdgeIndex].getPairIndex();   // half-edge pointing to the vertex to be removed
    thirdVertexIndex = *halfEdges[pairHalfEdgeIndex].getFromIndex();    // index of third vertex in the second triangle to be removed

    *vertices[thirdVertexIndex].getHalfEdgeIndex() = pairHalfEdgeIndex;

    // relabel half-edge pairs and delete/merge junctions

    long int fromHalfEdgeIndex, toHalfEdgeIndex;

    fromHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPreviousIndex();       // half-edge pointing to the vertex to be removed in the first triangle to be removed
    fromHalfEdgeIndex = *halfEdges[fromHalfEdgeIndex].getPairIndex();       // half-edge pointing from the vertex to be removed paired to the first triangle to be removed

    toHalfEdgeIndex = *halfEdges[halfEdgeIndex].getNextIndex();             // half-edge pointing from the merge vertex in the first triangle to be removed
    toHalfEdgeIndex = *halfEdges[toHalfEdgeIndex].getPairIndex();           // half-edge pointing to the vertex to be removed paired to the first triangle to be removed

    *halfEdges[fromHalfEdgeIndex].getPairIndex() = toHalfEdgeIndex;
    *halfEdges[toHalfEdgeIndex].getPairIndex() = fromHalfEdgeIndex;

    fromHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPairIndex();           // half-edge pointing from the merge vertex to the vertex to be removed in the second triangle to be removed
    fromHalfEdgeIndex = *halfEdges[fromHalfEdgeIndex].getPreviousIndex();   // half-edge pointing to the merge vertex in the second triangle to be removed
    fromHalfEdgeIndex = *halfEdges[fromHalfEdgeIndex].getPairIndex();       // half-edge pointing from the merge vertex paired to the second triangle to be removed

    toHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPairIndex();             // half-edge pointing from the merge vertex to the vertex to be removed in the second triangle to be removed
    toHalfEdgeIndex = *halfEdges[toHalfEdgeIndex].getNextIndex();           // half-edge pointing from the vertex to be removed in the second triangle to be removed
    toHalfEdgeIndex = *halfEdges[toHalfEdgeIndex].getPairIndex();           // half-edge pointing to the vertex to be removed paired to the second triangle to be removed

    *halfEdges[fromHalfEdgeIndex].getPairIndex() = toHalfEdgeIndex;
    *halfEdges[toHalfEdgeIndex].getPairIndex() = fromHalfEdgeIndex;

    // move merge vertex

    std::vector<double> const diff =
        wrapTo(toMergeIndex, fromMergeIndex, false);
    for (int dim=0; dim < 2; dim++) {
        vertices[toMergeIndex].getPosition()[dim] += diff[dim]/2;   // move vertex to half point between two merged vertices
    }

    // delete faces, junctions, half-edges, and vertex

    assert(halfEdgeIndex == halfEdges[halfEdgeIndex].getIndex());   // index of half-edge from deleted vertex
    pairHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPairIndex();   // index of half-edge to deleted vertex

    long int previousIndex, nextIndex;

    // first face
    previousIndex = *halfEdges[halfEdgeIndex].getPreviousIndex();
    nextIndex = *halfEdges[halfEdgeIndex].getNextIndex();
    faces.erase({halfEdgeIndex, previousIndex, nextIndex});             // delete face (first triangle)
    halfEdges.erase(previousIndex); halfEdges.erase(nextIndex);         // delete all half-edges from first triangle but the one from the erased junction

    // second face
    previousIndex = *halfEdges[pairHalfEdgeIndex].getPreviousIndex();
    nextIndex = *halfEdges[pairHalfEdgeIndex].getNextIndex();
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

    long int const vertexIndex = *halfEdges[halfEdgeIndex0].getFromIndex();
    assert(vertexIndex == *halfEdges[halfEdgeIndex1].getFromIndex());   // check that both half-edges go out of the same vertex
    *vertices[vertexIndex].getHalfEdgeIndex() = halfEdgeIndex0;         // re-assign identifying half-edge to `halfEdgeIndex0' which will remain attached to `vertexIndex'

    long int const newVertexIndex = maxKey(vertices);
    vertices.emplace(newVertexIndex, Vertex(
        // vertices is std::map so add with std::map::emplace
        newVertexIndex,                         // vertexIndex
        vertices[vertexIndex].getPosition(),    // position
        halfEdgeIndex1));                       // halfEdgeIndex
    sPVertices[newVertexIndex] = {
        // sPVertices is MultiIntKeyDict<SPVertex> so add with initialiser-list (sent to SPVertex)
        newVertexIndex};                        // vertexIndex

    // relabel origins and destinations of half of the half-edges

    long int halfEdgeIndex = halfEdgeIndex1;
    long int previousHalfEdgeIndex;
    while (halfEdgeIndex != halfEdgeIndex0) {

        assert(vertexIndex == *halfEdges[halfEdgeIndex].getFromIndex());
        *halfEdges[halfEdgeIndex].getFromIndex() = newVertexIndex;

        previousHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPreviousIndex();
        assert(vertexIndex == *halfEdges[previousHalfEdgeIndex].getToIndex());
        *halfEdges[previousHalfEdgeIndex].getToIndex() = newVertexIndex;

        halfEdgeIndex = *halfEdges[previousHalfEdgeIndex].getPairIndex();
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
        *halfEdges[halfEdgeIndex0].getToIndex(),        // fromIndex
        vertexIndex,                                    // toIndex
        previousNewHalfEdgeIndex0,                      // previousIndex
        nextNewHalfEdgeIndex0,                          // nextIndex
        halfEdgeIndex0));                               // pairIndex
    halfEdges.emplace(previousNewHalfEdgeIndex0, HalfEdge(
        previousNewHalfEdgeIndex0,                      // halfEdgeIndex
        newVertexIndex,                                 // fromIndex
        *halfEdges[halfEdgeIndex0].getToIndex(),        // toIndex
        nextNewHalfEdgeIndex0,                          // previousIndex
        newHalfEdgeIndex0,                              // nextIndex
        *halfEdges[halfEdgeIndex0].getPairIndex()));    // pairIndex
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
    *halfEdges[*halfEdges[halfEdgeIndex0].getPairIndex()].getPairIndex() =
        previousNewHalfEdgeIndex0;                      // relabel pair of old pair of halfEdgeIndex0
    *halfEdges[halfEdgeIndex0].getPairIndex() =
        newHalfEdgeIndex0;                              // relabel pair of halfEdgeIndex0

    // second face
    halfEdges.emplace(newHalfEdgeIndex1, HalfEdge(
        newHalfEdgeIndex1,                              // halfEdgeIndex
        *halfEdges[halfEdgeIndex1].getToIndex(),        // fromIndex
        newVertexIndex,                                 // toIndex
        previousNewHalfEdgeIndex1,                      // previousIndex
        nextNewHalfEdgeIndex1,                          // nextIndex
        halfEdgeIndex1));                               // pairIndex
    halfEdges.emplace(previousNewHalfEdgeIndex1, HalfEdge(
        previousNewHalfEdgeIndex1,                      // halfEdgeIndex
        vertexIndex,                                    // fromIndex
        *halfEdges[halfEdgeIndex1].getToIndex(),        // toIndex
        nextNewHalfEdgeIndex1,                          // previousIndex
        newHalfEdgeIndex1,                              // nextIndex
        *halfEdges[halfEdgeIndex1].getPairIndex()));    // pairIndex
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
    *halfEdges[*halfEdges[halfEdgeIndex1].getPairIndex()].getPairIndex() =
        previousNewHalfEdgeIndex1;                      // relabel pair of old pair of halfEdgeIndex1
    *halfEdges[halfEdgeIndex1].getPairIndex() =
        newHalfEdgeIndex1;                              // relabel pair of halfEdgeIndex1

    // create junction

    junctions[{nextNewHalfEdgeIndex0, nextNewHalfEdgeIndex1}] = {
        // faces is MultiIntKeyDict<Junction> so add with initialiser-list (sent to Junction)
        nextNewHalfEdgeIndex0};    // halfEdgeIndex

    // move vertices apart

    std::vector<double> direction = {std::cos(angle), std::sin(angle)};
    if (cross2(getHalfEdgeVector(newHalfEdgeIndex0), direction) < 0) {
        for (int dim=0; dim < 2; dim++) { direction[dim] *= -1; }   // enforces that the newly created face has anticlockwise orientation
    }

    for (int dim=0; dim < 2; dim++) {
        vertices[vertexIndex].getPosition()[dim] -=
            direction[dim]*length/2.;
        vertices[newVertexIndex].getPosition()[dim] +=
            direction[dim]*length/2.;
    }

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
    double position[2];
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
                position,                               // position of the vertex on the regular triangular lattice
                6*vertexIndex + 0));                    // index of a single half-edge going out of this vertex
            wrap(vertices[vertexIndex].getPosition());  // wrap coordinate around periodic boundary condition

            // create cell or self-propelled vertex
            if ((line - column)%3 == 0) {   // condition for vertex to be a cell centre...
                cells[vertexIndex] = {
                    // cells is MultiIntKeyDict<Cell> so add with initialiser-list (sent to Cell)
                    vertexIndex,                    // vertexIndex
                    (junctionLength*junctionLength) //A0
                        *(3./2.)/std::tan(std::numbers::pi/6)};
            }
            else {                          // ... or a self-propelled vertex
                sPVertices[vertexIndex] = {
                    // sPVertices is MultiIntKeyDict<SPVertex> so add with initialiser-list (sent to SPVertex)
                    vertexIndex};
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

