#include <cmath>
#include <numbers>
#include <assert.h>

#include "system.hpp"
#include "tools.hpp"

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
    std::vector<long int> proxyCellPos(2, 0);   // these two proxies are there to always exactly compare integer positions (and not approximately compare floats)
    std::vector<long int> proxyVertexPos(2, 0);
    std::vector<std::vector<long int>> proxyIncrements =
        {{0, -2}, {1, -1}, {1, 1}, {0, 2}, {-1, 1}, {-1, -1}};
    std::map<std::vector<long int>, long int> vertexPosMap;
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
            proxyCellPos[0] = 2*col + col%2;
            proxyCellPos[1] = 3*line;
            for (long int k=1; k <= 6; k++) {
                for (int dim=0; dim < 2; dim++) {
                    proxyVertexPos[dim] =
                        proxyCellPos[dim] + proxyIncrements[k][dim];
                }
                position[0] =
                    vertices[vertexIndex(col, line, 0)].getPosition()[0]
                        + cos(-std::numbers::pi/2. + k*std::numbers::pi/6.);
                position[1] =
                    vertices[vertexIndex(col, line, 0)].getPosition()[1]
                        + sin(-std::numbers::pi/2. + k*std::numbers::pi/6.);

                if (vertexPosMap.find(proxyVertexPos) != vertexPosMap.end()) {  // vertex already exists
                    // relations between vertices for initialisation
                    vertexIndexMap[vertexIndex(col, line, k)] =
                        vertexPosMap.find(proxyVertexPos)->second;
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
                    vertexPosMap[proxyVertexPos] =
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

