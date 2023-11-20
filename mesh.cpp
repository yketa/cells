#include <iostream>
#include <cmath>
#include <assert.h>

#include "tools.hpp"
#include "mesh.hpp"

void Mesh::wrap(double* position) {

    for (int dim=0; dim < 2; dim++) {
        position[dim] = pmod(position[dim], systemSize[dim]);
    }
}

std::vector<double> const Mesh::wrapDiff(
    double* const& fromPos, double* const& toPos) {

    std::vector<double> disp(2, 0);
    for (int dim=0; dim < 2; dim++) {
        disp[dim] = std::remainder(toPos[dim] - fromPos[dim], systemSize[dim]);
    }
    return disp;
}

std::vector<double> const Mesh::wrapTo(
    long int const& fromVertexIndex, long int const& toVertexIndex,
    bool const& unit) {

    std::vector<double> fromTo = wrapDiff(
        vertices[fromVertexIndex].getPosition(),
        vertices[toVertexIndex].getPosition());

    if (unit) {
        double norm = std::sqrt(fromTo[0]*fromTo[0] + fromTo[1]*fromTo[1]);
        for (int dim=0; dim < 2; dim++) {
            fromTo[dim] /= norm;
        }
    }

    return fromTo;
}

std::vector<double> const Mesh::getHalfEdgeVector(
    long int const& halfEdgeIndex, bool const& unit) {

    return wrapTo(
        *halfEdges[halfEdgeIndex].getFromIndex(),
        *halfEdges[halfEdgeIndex].getToIndex(),
        unit);
}

double const Mesh::getEdgeLength(long int const& halfEdgeIndex) {

    std::vector<double> const halfEdgeVector =
        getHalfEdgeVector(halfEdgeIndex, false);
    return std::sqrt(
        halfEdgeVector[0]*halfEdgeVector[0]
        + halfEdgeVector[1]*halfEdgeVector[1]);
}

std::vector<std::vector<long int>> const Mesh::getNeighbourVertices(
    long int const& vertexIndex) {

    std::vector<long int> neighbourVerticesIndices(0);
    std::vector<long int> halfEdgesToNeighboursIndices(0);

    // find destination vertex in half-edge construction
    long int halfEdgeIndex =
        *vertices[vertexIndex].getHalfEdgeIndex();
    assert(halfEdgeIndex >= 0);                                         // check that the half-edge exists
    assert(*halfEdges[halfEdgeIndex].getFromIndex() == vertexIndex);    // check that the half-edge goes out of this vertex
    long int const firstNeighbourVertexIndex =
        *halfEdges[halfEdgeIndex].getToIndex();
    assert(firstNeighbourVertexIndex >= 0);                             // check that the first destination vertex exists

    // loop around neighbours
    long int toVertexIndex = -1;
    long int previousHalfEdgeIndex;
    while (toVertexIndex != firstNeighbourVertexIndex) {

        // get next half-edge in anticlockwise order (previous - pair)
        previousHalfEdgeIndex = *halfEdges[halfEdgeIndex].getPreviousIndex();
        halfEdgeIndex = *halfEdges[previousHalfEdgeIndex].getPairIndex();

        // find destination vertex in half-edge construction
        assert(halfEdgeIndex >= 0);                                         // check that the half-edge exists
        assert(*halfEdges[halfEdgeIndex].getFromIndex() == vertexIndex);    // check that the half-edge goes out of this vertex
        toVertexIndex = *halfEdges[halfEdgeIndex].getToIndex();
        assert(toVertexIndex >= 0);                                         // check that the destination vertex exists

        neighbourVerticesIndices.push_back(toVertexIndex);
        halfEdgesToNeighboursIndices.push_back(halfEdgeIndex);
    }

    return {neighbourVerticesIndices, halfEdgesToNeighboursIndices};
}

double const Mesh::getVertexToNeighboursArea(
    long int const& vertexIndex) {

    std::vector<long int> const halfEdgeIndices =   // neighbours (which should be in anticlockwise order)
        getNeighbourVertices(vertexIndex)[1];
    long int const numberNeighbours = halfEdgeIndices.size();
    assert(numberNeighbours >= 3);

    // compute area
    double area = 0;
    for (long int i=0; i < numberNeighbours; i++) {
        area += cross2( // area of triangle
            getHalfEdgeVector(
                halfEdgeIndices[pmod(i + 0, numberNeighbours)]),
            getHalfEdgeVector(
                halfEdgeIndices[pmod(i + 1, numberNeighbours)]))
            /2.;
    }

    return area;
}

double const Mesh::getVertexToNeighboursPerimeter(
    long int const& vertexIndex) {

    std::vector<long int> const vertexIndices =     // neighbours (which should be in anticlockwise order)
        getNeighbourVertices(vertexIndex)[0];
    long int const numberNeighbours = vertexIndices.size();
    assert(numberNeighbours >= 3);

    // compute perimeter
    double perimeter = 0;
    std::vector<double> diff;
    for (long int i=0; i < numberNeighbours; i++) {
        diff = wrapTo(
            vertexIndices[pmod(i + 0, numberNeighbours)],
            vertexIndices[pmod(i + 1, numberNeighbours)],
            false);
        perimeter += std::sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return perimeter;
}

void Mesh::checkMesh() {

    std::vector<long int> vertexIndices(0);     // vector of vertex indices
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        vertexIndices.push_back(it->first);
    }
    std::vector<long int> halfEdgeIndices(0);   // vector of half-edge indices
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
        halfEdgeIndices.push_back(it->first);
    }

    long int halfEdgeIndex, halfEdgeIndexBis, pairHalfEdgeIndex;
    long int fromVertex, toVertex;
    std::vector<long int> triangle(3, 0);
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {  // loop over all half-edges
        halfEdgeIndex = it->first;

        if (!inVec<long int>(halfEdgeIndices, halfEdgeIndex)) { continue; }

        triangle = {    // three half-edges forming a triangle
            halfEdgeIndex,
            *halfEdges[halfEdgeIndex].getNextIndex(),
            *halfEdges[halfEdgeIndex].getPreviousIndex()};

        for (int i=0; i < 3; i++) {

            fromVertex = *halfEdges[triangle[i]].getFromIndex();
            toVertex = *halfEdges[triangle[i]].getToIndex();

            if (inVec<long int>(vertexIndices, fromVertex)) {       // looping over all half-edges should enconter all vertices

                halfEdgeIndexBis = *vertices[fromVertex].getHalfEdgeIndex();
                assert(             // check consistency with identifying half-edge
                    *halfEdges[halfEdgeIndexBis].getFromIndex() == fromVertex);

                eraseInVec<long int>(vertexIndices, fromVertex);
            }

            assert(cross2(          // check that the triangle has anticlockwise orientation
                getHalfEdgeVector(triangle[pmod(i + 0, 3)]),
                getHalfEdgeVector(triangle[pmod(i + 1, 3)])) > 0);

            pairHalfEdgeIndex = *halfEdges[triangle[i]].getPairIndex();
            assert(triangle[i] ==   // check the indexing of pair
                *halfEdges[pairHalfEdgeIndex].getPairIndex());
            assert(toVertex ==
                *halfEdges[pairHalfEdgeIndex].getFromIndex());
            assert(fromVertex ==
                *halfEdges[pairHalfEdgeIndex].getToIndex());

            assert(triangle[i] ==   // check the indexing of previous
                *halfEdges[triangle[pmod(i + 1, 3)]].getPreviousIndex());
            assert(toVertex ==
                *halfEdges[triangle[pmod(i + 1, 3)]].getFromIndex());

            assert(triangle[i] ==   // check the indexing of next
                *halfEdges[triangle[pmod(i - 1, 3)]].getNextIndex());
            assert(fromVertex ==
                *halfEdges[triangle[pmod(i - 1, 3)]].getToIndex());

            eraseInVec<long int>(halfEdgeIndices, triangle[i]);
        }
    }

    assert(vertexIndices.size() == 0);
    assert(halfEdgeIndices.size() == 0);
    std::cerr << "OK mesh construction" << std::endl;
}

