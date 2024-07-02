#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>

#include "assert.hpp"
#include "mesh.hpp"
#include "tools.hpp"

std::vector<double> Mesh::wrap(
    std::vector<double> const& position)
    const {

    std::vector<double> wposition(2, 0);
    for (int dim=0; dim < 2; dim++) {
        wposition[dim] = pmod(position[dim], systemSize[dim]);
    }
    return wposition;
}

std::vector<double> Mesh::wrapDiff(
    std::vector<double> const& fromPos,
    std::vector<double> const& toPos,
    bool const& unit)
    const {

    // compute difference bector
    std::vector<double> disp(2, 0);
    for (int dim=0; dim < 2; dim++) {
        disp[dim] = std::remainder(toPos[dim] - fromPos[dim], systemSize[dim]);
    }

    // normalise
    if (unit) {
        double const norm = std::sqrt(disp[0]*disp[0] + disp[1]*disp[1]);
        for (int dim=0; dim < 2; dim++) {
            disp[dim] /= norm;
        }
    }

    return disp;
}

std::vector<double> Mesh::wrapTo(
    long int const& fromVertexIndex, long int const& toVertexIndex,
    bool const& unit)
    const {

    return wrapDiff(
        vertices.at(fromVertexIndex).getPosition(),
        vertices.at(toVertexIndex).getPosition(),
        unit);
}

long int Mesh::getHalfEdgeIndex(
    long int const& fromVertexIndex, long int const& toVertexIndex)
    const {

    std::vector<std::vector<long int>> neighbours =     // neighbours and half-edge to neighbours
        getNeighbourVertices(fromVertexIndex);
    long int const nNeighbours = neighbours[0].size();  // number of neighbours

    for (long int i=0; i < nNeighbours; i++) {  // loop over neighbours
        if (neighbours[0][i] == toVertexIndex) {
            return neighbours[1][i];
        }
    }
    throw std::runtime_error("There is not half-edge from vertex index "
        + std::to_string(fromVertexIndex) + " to vertex index "
        + std::to_string(toVertexIndex) + ".");
}

std::vector<double> Mesh::getHalfEdgeVector(
    long int const& halfEdgeIndex, bool const& unit)
    const {

    return wrapTo(
        halfEdges.at(halfEdgeIndex).getFromIndex(),
        halfEdges.at(halfEdgeIndex).getToIndex(),
        unit);
}

double Mesh::getEdgeLength(long int const& halfEdgeIndex) const {

    std::vector<double> const halfEdgeVector =
        getHalfEdgeVector(halfEdgeIndex, false);
    return std::sqrt(
        halfEdgeVector[0]*halfEdgeVector[0]
        + halfEdgeVector[1]*halfEdgeVector[1]);
}

std::vector<std::vector<long int>> Mesh::getNeighbourVertices(
    long int const& vertexIndex)
    const {

    std::vector<long int> neighbourVerticesIndices(0);
    std::vector<long int> halfEdgesToNeighboursIndices(0);

    // find destination vertex in half-edge construction
    long int halfEdgeIndex =
        vertices.at(vertexIndex).getHalfEdgeIndex();
    assert(halfEdgeIndex >= 0);                                             // check that the half-edge exists
    assert(halfEdges.at(halfEdgeIndex).getFromIndex() == vertexIndex);      // check that the half-edge goes out of this vertex
    long int const firstNeighbourVertexIndex =
        halfEdges.at(halfEdgeIndex).getToIndex();
    assert(firstNeighbourVertexIndex >= 0);                                 // check that the first destination vertex exists

    // loop around neighbours
    long int toVertexIndex;
    long int previousHalfEdgeIndex;
    while (true) {

        // get next half-edge in anticlockwise order (previous - pair)
        previousHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPreviousIndex();
        halfEdgeIndex = halfEdges.at(previousHalfEdgeIndex).getPairIndex();

        // find destination vertex in half-edge construction
        assert(halfEdgeIndex >= 0);                                         // check that the half-edge exists
        assert(halfEdges.at(halfEdgeIndex).getFromIndex() == vertexIndex);  // check that the half-edge goes out of this vertex
        toVertexIndex = halfEdges.at(halfEdgeIndex).getToIndex();
        assert(toVertexIndex >= 0);                                         // check that the destination vertex exists

        neighbourVerticesIndices.push_back(toVertexIndex);
        halfEdgesToNeighboursIndices.push_back(halfEdgeIndex);

        if (toVertexIndex == firstNeighbourVertexIndex) {   // all neighbours have been found
            break;
        }
    }

    return {neighbourVerticesIndices, halfEdgesToNeighboursIndices};
}

double Mesh::getVertexToNeighboursArea(
    long int const& vertexIndex)
    const {

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

double Mesh::getVertexToNeighboursPerimeter(
    long int const& vertexIndex)
    const {

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

void Mesh::moveToNeigboursBarycentre(
    long int const& vertexIndex) {

    // compute centre of mass
    std::vector<double> posCM = {0, 0};             // position of centre of mass (= barycentre)
    std::vector<long int> const neighbourIndices =  // indices of neighbouring vertices
        getNeighbourVertices(vertexIndex)[0];
    long int N = neighbourIndices.size();           // number of vertices
    for (long int neighbourIndex : neighbourIndices) {
        std::vector<double> const pos =
            (vertices.at(neighbourIndex)).getPosition();
        for (int dim=0; dim < 2; dim++)
            { posCM[dim] += pos.at(dim)/N; }
    }

    // set position
    vertices[vertexIndex].setPosition(  // wrap position
        wrap(posCM));
    vertices[vertexIndex].setUPosition( // set unwrapped position back to wrapped position
        vertices.at(vertexIndex).getPosition());
}

TopoChangeEdgeInfoType Mesh::deleteEdge(
    long int const& halfEdgeIndex) {

    std::vector<long int> deletedHalfEdgeIndices(0);

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

    // relabel half-edge pairs

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

    // delete half-edges and vertex

    assert(halfEdgeIndex == halfEdges.at(halfEdgeIndex).getIndex());    // index of half-edge from deleted vertex
    pairHalfEdgeIndex = halfEdges.at(halfEdgeIndex).getPairIndex();     // index of half-edge to deleted vertex

    long int previousIndex, nextIndex;

    // first face
    previousIndex = halfEdges.at(halfEdgeIndex).getPreviousIndex();
    nextIndex = halfEdges.at(halfEdgeIndex).getNextIndex();
    // delete all half-edges from first triangle but the one from the deleted edge
    deletedHalfEdgeIndices.push_back(previousIndex);
    deletedHalfEdgeIndices.push_back(nextIndex);

    // second face
    previousIndex = halfEdges.at(pairHalfEdgeIndex).getPreviousIndex();
    nextIndex = halfEdges.at(pairHalfEdgeIndex).getNextIndex();
    // delete all half-edges from second triangle but the one from the deleted edge
    deletedHalfEdgeIndices.push_back(previousIndex);
    deletedHalfEdgeIndices.push_back(nextIndex);

    // delete half-edges in the deleted edge
    deletedHalfEdgeIndices.push_back(halfEdgeIndex);
    deletedHalfEdgeIndices.push_back(pairHalfEdgeIndex);

    vertices.erase(fromMergeIndex);                                         // delete vertex
    for (long int index : deletedHalfEdgeIndices) halfEdges.erase(index);   // delete half-edges

    return std::make_tuple(fromMergeIndex, deletedHalfEdgeIndices);
}

TopoChangeEdgeInfoType Mesh::createEdge(
    long int const& halfEdgeIndex0, long int const& halfEdgeIndex1,
    double const& angle, double const& length,
    std::string type0, std::string type1) {

    std::vector<long int> createdHalfEdgeIndices(0);

    // create new vertex

    long int const vertexIndex = halfEdges.at(halfEdgeIndex0).getFromIndex();
    assert(vertexIndex == halfEdges.at(halfEdgeIndex1).getFromIndex()); // check that both half-edges go out of the same vertex
    vertices[vertexIndex].setHalfEdgeIndex(halfEdgeIndex0);             // re-assign identifying half-edge to `halfEdgeIndex0' which will remain attached to `vertexIndex'

    long int const newVertexIndex = maxKey(vertices) + 1;
    std::vector<double> const position =
        vertices.at(vertexIndex).getPosition();
    vertices.emplace(newVertexIndex, Vertex(
        // vertices is std::map so add with std::map::emplace
        newVertexIndex,                     // vertexIndex
        position,                           // position
        halfEdgeIndex1,                     // halfEdgeIndex
        false,                              // boundary
        vertices[vertexIndex].getType()));  // type (inherit from parent vertex)

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

    // create new half-edges

    Counter c(maxKey(halfEdges) + 1);

    long int const newHalfEdgeIndex0 = c();
    createdHalfEdgeIndices.push_back(newHalfEdgeIndex0);
    long int const previousNewHalfEdgeIndex0 = c();
    createdHalfEdgeIndices.push_back(previousNewHalfEdgeIndex0);
    long int const nextNewHalfEdgeIndex0 = c();
    createdHalfEdgeIndices.push_back(nextNewHalfEdgeIndex0);
    long int const newHalfEdgeIndex1 = c();
    createdHalfEdgeIndices.push_back(newHalfEdgeIndex1);
    long int const previousNewHalfEdgeIndex1 = c();
    createdHalfEdgeIndices.push_back(previousNewHalfEdgeIndex1);
    long int const nextNewHalfEdgeIndex1 = c();
    createdHalfEdgeIndices.push_back(nextNewHalfEdgeIndex1);

    // first triangle
    halfEdges.emplace(newHalfEdgeIndex0, HalfEdge(
        // halfEdges is std::map so add with std::map::emplace
        newHalfEdgeIndex0,                                                  // halfEdgeIndex
        halfEdges.at(halfEdgeIndex0).getToIndex(),                          // fromIndex
        vertexIndex,                                                        // toIndex
        previousNewHalfEdgeIndex0,                                          // previousIndex
        nextNewHalfEdgeIndex0,                                              // nextIndex
        halfEdgeIndex0,                                                     // pairIndex
        halfEdges[halfEdgeIndex0].getType()));                              // type
    halfEdges.emplace(previousNewHalfEdgeIndex0, HalfEdge(
        previousNewHalfEdgeIndex0,                                          // halfEdgeIndex
        newVertexIndex,                                                     // fromIndex
        halfEdges.at(halfEdgeIndex0).getToIndex(),                          // toIndex
        nextNewHalfEdgeIndex0,                                              // previousIndex
        newHalfEdgeIndex0,                                                  // nextIndex
        halfEdges.at(halfEdgeIndex0).getPairIndex(),                        // pairIndex
        halfEdges[halfEdges.at(halfEdgeIndex0).getPairIndex()].getType())); // type
    halfEdges.emplace(nextNewHalfEdgeIndex0, HalfEdge(
        nextNewHalfEdgeIndex0,                                              // halfEdgeIndex
        vertexIndex,                                                        // fromIndex
        newVertexIndex,                                                     // toIndex
        newHalfEdgeIndex0,                                                  // previousIndex
        previousNewHalfEdgeIndex0,                                          // nextIndex
        nextNewHalfEdgeIndex1,                                              // pairIndex
        type0));                                                            // type
    halfEdges[halfEdges.at(halfEdgeIndex0).getPairIndex()].setPairIndex(
        previousNewHalfEdgeIndex0);                                         // relabel pair of old pair of halfEdgeIndex0
    halfEdges[halfEdgeIndex0].setPairIndex(
        newHalfEdgeIndex0);                                                 // relabel pair of halfEdgeIndex0

    // second triangle
    halfEdges.emplace(newHalfEdgeIndex1, HalfEdge(
        newHalfEdgeIndex1,                                                  // halfEdgeIndex
        halfEdges.at(halfEdgeIndex1).getToIndex(),                          // fromIndex
        newVertexIndex,                                                     // toIndex
        previousNewHalfEdgeIndex1,                                          // previousIndex
        nextNewHalfEdgeIndex1,                                              // nextIndex
        halfEdgeIndex1,                                                     // pairIndex
        halfEdges[halfEdgeIndex1].getType()));                              // type
    halfEdges.emplace(previousNewHalfEdgeIndex1, HalfEdge(
        previousNewHalfEdgeIndex1,                                          // halfEdgeIndex
        vertexIndex,                                                        // fromIndex
        halfEdges.at(halfEdgeIndex1).getToIndex(),                          // toIndex
        nextNewHalfEdgeIndex1,                                              // previousIndex
        newHalfEdgeIndex1,                                                  // nextIndex
        halfEdges.at(halfEdgeIndex1).getPairIndex(),                        // pairIndex
        halfEdges[halfEdges.at(halfEdgeIndex1).getPairIndex()].getType())); // type
    halfEdges.emplace(nextNewHalfEdgeIndex1, HalfEdge(
        nextNewHalfEdgeIndex1,                                              // halfEdgeIndex
        newVertexIndex,                                                     // fromIndex
        vertexIndex,                                                        // toIndex
        newHalfEdgeIndex1,                                                  // previousIndex
        previousNewHalfEdgeIndex1,                                          // nextIndex
        nextNewHalfEdgeIndex0,                                              // pairIndex
        type1));                                                            // type
    halfEdges[halfEdges.at(halfEdgeIndex1).getPairIndex()].setPairIndex(
        previousNewHalfEdgeIndex1);                                         // relabel pair of old pair of halfEdgeIndex1
    halfEdges[halfEdgeIndex1].setPairIndex(
        newHalfEdgeIndex1);                                                 // relabel pair of halfEdgeIndex1

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

    return std::make_tuple(newVertexIndex, createdHalfEdgeIndices);
}

void Mesh::scale(double const& scalingFactor) {

    auto scaleVector = [&scalingFactor](std::vector<double>& vector)
        { for (int dim=0; dim < 2; dim++) vector[dim] *= scalingFactor; };

    // scale box
    scaleVector(systemSize);

    // scale vertices
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        std::vector<double> position = (it->second).getPosition();
        scaleVector(position);
        (it->second).setPosition(position);
        std::vector<double> uposition = (it->second).getUPosition();
        scaleVector(uposition);
        (it->second).setUPosition(uposition);
    }
}

void Mesh::checkMesh() const {

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
    bool boundary = false;
    long int thirdVertex;
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {  // loop over all half-edges
        halfEdgeIndex = it->first;

        if (!inVec<long int>(halfEdgeIndices, halfEdgeIndex)) { continue; }

        triangle = {    // three half-edges forming a triangle
            halfEdgeIndex,
            halfEdges.at(halfEdgeIndex).getNextIndex(),
            halfEdges.at(halfEdgeIndex).getPreviousIndex()};

        fromVertex =
            halfEdges.at(halfEdgeIndex)
                .getFromIndex();
        toVertex =
            halfEdges.at(halfEdgeIndex)
                .getToIndex();
        thirdVertex =
            halfEdges.at(halfEdges.at(halfEdgeIndex).getNextIndex())
                .getToIndex();
        boundary = (    // does the triangle has an outer boundary corner
            vertices.at(fromVertex).getBoundary()
            || vertices.at(toVertex).getBoundary()
            || vertices.at(thirdVertex).getBoundary());

        for (int i=0; i < 3; i++) {

            fromVertex = halfEdges.at(triangle[i]).getFromIndex();
            toVertex = halfEdges.at(triangle[i]).getToIndex();

            if (inVec<long int>(vertexIndices, fromVertex)) {       // looping over all half-edges should enconter all vertices

                halfEdgeIndexBis = vertices.at(fromVertex).getHalfEdgeIndex();
                assert(             // check consistency with identifying half-edge
                    halfEdges.at(halfEdgeIndexBis).getFromIndex()
                        == fromVertex);

                eraseInVec<long int>(vertexIndices, fromVertex);
            }

            if (!boundary) {        // only for inner triangles
                assert(cross2(      // check that the triangle has anticlockwise orientation
                    getHalfEdgeVector(triangle[pmod(i + 0, 3)]),
                    getHalfEdgeVector(triangle[pmod(i + 1, 3)])) > 0);
            }

            pairHalfEdgeIndex = halfEdges.at(triangle[i]).getPairIndex();
            assert(triangle[i] ==   // check the indexing of pair
                halfEdges.at(pairHalfEdgeIndex).getPairIndex());
            assert(toVertex ==
                halfEdges.at(pairHalfEdgeIndex).getFromIndex());
            assert(fromVertex ==
                halfEdges.at(pairHalfEdgeIndex).getToIndex());

            assert(triangle[i] ==   // check the indexing of next
                halfEdges.at(triangle[pmod(i + 1, 3)]).getPreviousIndex());
            assert(toVertex ==
                halfEdges.at(triangle[pmod(i + 1, 3)]).getFromIndex());

            assert(triangle[i] ==   // check the indexing of previous
                halfEdges.at(triangle[pmod(i - 1, 3)]).getNextIndex());
            assert(fromVertex ==
                halfEdges.at(triangle[pmod(i - 1, 3)]).getToIndex());

            eraseInVec<long int>(halfEdgeIndices, triangle[i]);
        }
    }

    assert(vertexIndices.size() == 0);
    assert(halfEdgeIndices.size() == 0);
    std::cerr << "OK mesh construction" << std::endl;
}

void cerrTopoChangeEdgeInfo(TopoChangeEdgeInfoType const& info) {
    std::cerr << "vertex:     " << std::get<0>(info) << std::endl;
    std::cerr << "half-edges: ";
    for (long int i : std::get<1>(info)) { std::cerr << i << " "; }
    std::cerr << std::endl;
}

