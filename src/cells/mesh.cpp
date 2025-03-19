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

std::vector<double> Mesh::getHalfEdgeCentre(
    long int const& halfEdgeIndex)
    const {

    std::vector<double> const fromVertexPos =
        vertices.at(halfEdges.at(halfEdgeIndex).getFromIndex()).getPosition();
    std::vector<double> const fromToDiff =
        getHalfEdgeVector(halfEdgeIndex);

    return wrap({
        fromVertexPos[0] + fromToDiff[0]/2,
        fromVertexPos[1] + fromToDiff[1]/2});
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

void Mesh::moveToNeighboursBarycentre(
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

void Mesh::moveToNeighboursCentroid(
    long int const& vertexIndex) {

    double const area = getVertexToNeighboursArea(vertexIndex);

    std::vector<long int> const neighbourVerticesIndices =
        getNeighbourVertices(vertexIndex)[0];
    long int const numberNeighbours = neighbourVerticesIndices.size();

    std::vector<double> const uposr = vertices.at(vertexIndex).getUPosition();  // reference unwrapped position
    std::vector<double> diff = {0, 0};                                          // difference from vertex position to reference position

    for (long int i=0; i < numberNeighbours; i++) {

        long int const index0 =             // index of first neighbour
            neighbourVerticesIndices[i];
        std::vector<double> const diff0 =   // wrapped difference from unwrapped reference position to unwrapped unwrapped position of first neighbour
            wrapDiff(uposr, vertices.at(index0).getUPosition());

        long int const index1 =             // index of second (next) neighbour
            neighbourVerticesIndices[pmod(i + 1, numberNeighbours)];
        std::vector<double> const diff1 =   // wrapped difference from unwrapped reference position to unwrapped unwrapped position of second neighbour
            wrapDiff(uposr, vertices.at(index1).getUPosition());

        for (int dim=0; dim < 2; dim++) {
            diff[dim] +=                    // position of centroid in the space where reference position is the centre
                (diff0[dim] + diff1[dim])*(
                    diff0[0]*diff1[1] - diff1[0]*diff0[1]
                )/(6*area);
        }
    }

    vertices[vertexIndex].setUPosition( // unwrapped position
        {uposr.at(0) + diff.at(0), uposr.at(1) + diff.at(1)});
    vertices[vertexIndex].setPosition(  // (wrapped) position
        wrap(vertices.at(vertexIndex).getUPosition()));
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
        newVertexIndex,                         // vertexIndex
        position,                               // position
        halfEdgeIndex1,                         // halfEdgeIndex
        false,                                  // boundary
        vertices.at(vertexIndex).getType()));   // type (inherit from parent vertex)

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
        newHalfEdgeIndex0,                                                      // halfEdgeIndex
        halfEdges.at(halfEdgeIndex0).getToIndex(),                              // fromIndex
        vertexIndex,                                                            // toIndex
        previousNewHalfEdgeIndex0,                                              // previousIndex
        nextNewHalfEdgeIndex0,                                                  // nextIndex
        halfEdgeIndex0,                                                         // pairIndex
        halfEdges.at(halfEdgeIndex0).getType()));                               // type
    halfEdges.emplace(previousNewHalfEdgeIndex0, HalfEdge(
        previousNewHalfEdgeIndex0,                                              // halfEdgeIndex
        newVertexIndex,                                                         // fromIndex
        halfEdges.at(halfEdgeIndex0).getToIndex(),                              // toIndex
        nextNewHalfEdgeIndex0,                                                  // previousIndex
        newHalfEdgeIndex0,                                                      // nextIndex
        halfEdges.at(halfEdgeIndex0).getPairIndex(),                            // pairIndex
        halfEdges.at(halfEdges.at(halfEdgeIndex0).getPairIndex()).getType()));  // type
    halfEdges.emplace(nextNewHalfEdgeIndex0, HalfEdge(
        nextNewHalfEdgeIndex0,                                                  // halfEdgeIndex
        vertexIndex,                                                            // fromIndex
        newVertexIndex,                                                         // toIndex
        newHalfEdgeIndex0,                                                      // previousIndex
        previousNewHalfEdgeIndex0,                                              // nextIndex
        nextNewHalfEdgeIndex1,                                                  // pairIndex
        type0));                                                                // type
    halfEdges[halfEdges.at(halfEdgeIndex0).getPairIndex()].setPairIndex(
        previousNewHalfEdgeIndex0);                                             // relabel pair of old pair of halfEdgeIndex0
    halfEdges[halfEdgeIndex0].setPairIndex(
        newHalfEdgeIndex0);                                                     // relabel pair of halfEdgeIndex0

    // second triangle
    halfEdges.emplace(newHalfEdgeIndex1, HalfEdge(
        newHalfEdgeIndex1,                                                      // halfEdgeIndex
        halfEdges.at(halfEdgeIndex1).getToIndex(),                              // fromIndex
        newVertexIndex,                                                         // toIndex
        previousNewHalfEdgeIndex1,                                              // previousIndex
        nextNewHalfEdgeIndex1,                                                  // nextIndex
        halfEdgeIndex1,                                                         // pairIndex
        halfEdges.at(halfEdgeIndex1).getType()));                               // type
    halfEdges.emplace(previousNewHalfEdgeIndex1, HalfEdge(
        previousNewHalfEdgeIndex1,                                              // halfEdgeIndex
        vertexIndex,                                                            // fromIndex
        halfEdges.at(halfEdgeIndex1).getToIndex(),                              // toIndex
        nextNewHalfEdgeIndex1,                                                  // previousIndex
        newHalfEdgeIndex1,                                                      // nextIndex
        halfEdges.at(halfEdgeIndex1).getPairIndex(),                            // pairIndex
        halfEdges.at(halfEdges.at(halfEdgeIndex1).getPairIndex()).getType()));  // type
    halfEdges.emplace(nextNewHalfEdgeIndex1, HalfEdge(
        nextNewHalfEdgeIndex1,                                                  // halfEdgeIndex
        newVertexIndex,                                                         // fromIndex
        vertexIndex,                                                            // toIndex
        newHalfEdgeIndex1,                                                      // previousIndex
        previousNewHalfEdgeIndex1,                                              // nextIndex
        nextNewHalfEdgeIndex0,                                                  // pairIndex
        type1));                                                                // type
    halfEdges[halfEdges.at(halfEdgeIndex1).getPairIndex()].setPairIndex(
        previousNewHalfEdgeIndex1);                                             // relabel pair of old pair of halfEdgeIndex1
    halfEdges[halfEdgeIndex1].setPairIndex(
        newHalfEdgeIndex1);                                                     // relabel pair of halfEdgeIndex1

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

TopoChangeEdgeInfoType Mesh::mergeVertices(
    long int const& fromVertexIndex, long int const& toVertexIndex) {

    std::vector<long int> deletedHalfEdgeIndices(0);

    // identify connecting edge

    long int index = -1;                                        // index of half-edge belonging to edge separating vertices
    std::vector<long int> const halfEdgesToNeighboursIndices =  // half-edge to neighbouring vertices
        getNeighbourVertices(fromVertexIndex)[1];
    for (long int halfEdgeIndex : halfEdgesToNeighboursIndices) {
        long int const nextHalfEdgeIndex =
            halfEdges.at(halfEdgeIndex).getNextIndex();
        long int const vertexIndex =
            halfEdges.at(
                halfEdges.at(
                    halfEdges.at(
                        nextHalfEdgeIndex
                        ).getPairIndex()
                    ).getNextIndex()
                ).getToIndex();
        if (vertexIndex == toVertexIndex) {
            index = nextHalfEdgeIndex;
            break;
        }
    }

    assert(index != -1);                                        // check there is actually an edge separating the vertices
    long int const sepHalfEdgeIndex = index;

    assert(                                                     // check vertices have the same bondary indicator
        vertices.at(fromVertexIndex).getBoundary()
            == vertices.at(toVertexIndex).getBoundary());
    assert(                                                     // check vertices are the same type
        vertices.at(fromVertexIndex).getType()
            == vertices.at(toVertexIndex).getType());

    // flag interior half-edges for deletion

    deletedHalfEdgeIndices.push_back(
        sepHalfEdgeIndex);
    deletedHalfEdgeIndices.push_back(
        halfEdges.at(sepHalfEdgeIndex).getPreviousIndex());
    deletedHalfEdgeIndices.push_back(
        halfEdges.at(sepHalfEdgeIndex).getNextIndex());

    deletedHalfEdgeIndices.push_back(
        halfEdges.at(sepHalfEdgeIndex).getPairIndex());
    deletedHalfEdgeIndices.push_back(
        halfEdges.at(
            halfEdges.at(sepHalfEdgeIndex).getPairIndex()
            ).getPreviousIndex());
    deletedHalfEdgeIndices.push_back(
        halfEdges.at(
            halfEdges.at(sepHalfEdgeIndex).getPairIndex()
            ).getNextIndex());

    // reassign half-edges origin and destination

    for (long int halfEdgeIndex : halfEdgesToNeighboursIndices) {
        halfEdges[
            halfEdgeIndex
            ].setFromIndex(toVertexIndex);
        halfEdges[
            halfEdges.at(halfEdgeIndex).getPreviousIndex()
            ].setToIndex(toVertexIndex);
    }

    // reassign half-edge pairs

    long int const inFromHalfEdgeIndex =
        halfEdges.at(
            halfEdges.at(
                sepHalfEdgeIndex
                ).getPreviousIndex()
            ).getPairIndex();
    long int const outFromHalfEdgeIndex =
        halfEdges.at(
            halfEdges.at(
                sepHalfEdgeIndex
                ).getNextIndex()
            ).getPairIndex();

    long int const inToHalfEdgeIndex =
        halfEdges.at(
            halfEdges.at(
                halfEdges.at(
                    sepHalfEdgeIndex
                    ).getPairIndex()
                ).getPreviousIndex()
            ).getPairIndex();
    long int const outToHalfEdgeIndex =
        halfEdges.at(
            halfEdges.at(
                halfEdges.at(
                    sepHalfEdgeIndex
                    ).getPairIndex()
                ).getNextIndex()
            ).getPairIndex();

    halfEdges[inFromHalfEdgeIndex].setPairIndex(outToHalfEdgeIndex);
    halfEdges[outToHalfEdgeIndex].setPairIndex(inFromHalfEdgeIndex);

    halfEdges[outFromHalfEdgeIndex].setPairIndex(inToHalfEdgeIndex);
    halfEdges[inToHalfEdgeIndex].setPairIndex(outFromHalfEdgeIndex);

    // set half-edge going out of vertex

    vertices[toVertexIndex].setHalfEdgeIndex(
        outToHalfEdgeIndex);
    vertices[halfEdges.at(outFromHalfEdgeIndex).getToIndex()].setHalfEdgeIndex(
        inToHalfEdgeIndex);
    vertices[halfEdges.at(outToHalfEdgeIndex).getToIndex()].setHalfEdgeIndex(
        inFromHalfEdgeIndex);

    // move vertex

    std::vector<double> const fromVertexUPosition =
        vertices.at(fromVertexIndex).getUPosition();
    std::vector<double> const toVertexUPosition =
        vertices.at(toVertexIndex).getUPosition();

    std::vector<double> const midUPosition = {
        (fromVertexUPosition[0] + toVertexUPosition[0])/2,
        (fromVertexUPosition[1] + toVertexUPosition[1])/2};

    vertices[toVertexIndex].setPosition(wrap(midUPosition));
    vertices[toVertexIndex].setUPosition(midUPosition);

    // delete vertex half-edges

    vertices.erase(fromVertexIndex);                                        // delete vertex
    for (long int index : deletedHalfEdgeIndices) halfEdges.erase(index);   // delete half-edges

    return std::make_tuple(fromVertexIndex, deletedHalfEdgeIndices);
}

TopoChangeEdgeInfoType Mesh::splitVertices(
    long int const& halfEdgeIndex) {

    std::vector<long int> createdHalfEdgeIndices(0);

    // create new vertex

    long int const fromVertexIndex =            // origin vertex of the half-edge index
        halfEdges.at(halfEdgeIndex).getFromIndex();
    long int const toVertexIndex =              // destination vertex of the half-edge index
        halfEdges.at(halfEdgeIndex).getToIndex();

    long int const newVertexIndex = maxKey(vertices) + 1;
    vertices.emplace(newVertexIndex, Vertex(
        // vertices is std::map so add with std::map::emplace
        newVertexIndex,                             // vertexIndex
        getHalfEdgeCentre(halfEdgeIndex),           // position
        halfEdgeIndex,                              // halfEdgeIndex
        false,                                      // boundary
        vertices.at(fromVertexIndex).getType()));   // type (inherit from origin vertex

    // create new half-edges

    long int const previousHalfEdgeIndex =
        halfEdges.at(halfEdgeIndex).getPreviousIndex();
    long int const nextHalfEdgeIndex =
        halfEdges.at(halfEdgeIndex).getNextIndex();
    long int const pairHalfEdgeIndex =
        halfEdges.at(halfEdgeIndex).getPairIndex();
    long int const previousPairHalfEdgeIndex =
        halfEdges.at(pairHalfEdgeIndex).getPreviousIndex();
    long int const nextPairHalfEdgeIndex =
        halfEdges.at(pairHalfEdgeIndex).getNextIndex();

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
        newHalfEdgeIndex0,                                                      // halfEdgeIndex
        fromVertexIndex,                                                        // fromIndex
        newVertexIndex,                                                         // toIndex
        previousHalfEdgeIndex,                                                  // previousIndex
        nextNewHalfEdgeIndex0,                                                  // nextIndex
        pairHalfEdgeIndex,                                                      // pairIndex
        halfEdges.at(halfEdgeIndex).getType()));                                // type
    halfEdges.emplace(previousNewHalfEdgeIndex0, HalfEdge(
        previousNewHalfEdgeIndex0,                                              // halfEdgeIndex
        halfEdges.at(nextHalfEdgeIndex).getToIndex(),                           // fromIndex
        newVertexIndex,                                                         // toIndex
        nextHalfEdgeIndex,                                                      // previousIndex
        halfEdgeIndex,                                                          // nextIndex
        nextNewHalfEdgeIndex0,                                                  // pairIndex
        halfEdges.at(nextHalfEdgeIndex).getType()));                            // type
    halfEdges.emplace(nextNewHalfEdgeIndex0, HalfEdge(
        nextNewHalfEdgeIndex0,                                                  // halfEdgeIndex
        newVertexIndex,                                                         // fromIndex
        halfEdges.at(nextHalfEdgeIndex).getToIndex(),                           // toIndex
        newHalfEdgeIndex0,                                                      // previousIndex
        previousHalfEdgeIndex,                                                  // nextIndex
        previousNewHalfEdgeIndex0,                                              // pairIndex
        halfEdges.at(nextHalfEdgeIndex).getType()));                            // type

    // second triangle
    halfEdges.emplace(newHalfEdgeIndex1, HalfEdge(
        newHalfEdgeIndex1,                                                      // halfEdgeIndex
        toVertexIndex,                                                          // fromIndex
        newVertexIndex,                                                         // toIndex
        previousPairHalfEdgeIndex,                                              // previousIndex
        nextNewHalfEdgeIndex1,                                                  // nextIndex
        halfEdgeIndex,                                                          // pairIndex
        halfEdges.at(pairHalfEdgeIndex).getType()));                            // type
    halfEdges.emplace(previousNewHalfEdgeIndex1, HalfEdge(
        previousNewHalfEdgeIndex1,                                              // halfEdgeIndex
        halfEdges.at(nextPairHalfEdgeIndex).getToIndex(),                       // fromIndex
        newVertexIndex,                                                         // toIndex
        nextPairHalfEdgeIndex,                                                  // previousIndex
        pairHalfEdgeIndex,                                                      // nextIndex
        nextNewHalfEdgeIndex1,                                                  // pairIndex
        halfEdges.at(nextPairHalfEdgeIndex).getType()));                        // type
    halfEdges.emplace(nextNewHalfEdgeIndex1, HalfEdge(
        nextNewHalfEdgeIndex1,                                                  // halfEdgeIndex
        newVertexIndex,                                                         // fromIndex
        halfEdges.at(nextPairHalfEdgeIndex).getToIndex(),                       // toIndex
        newHalfEdgeIndex1,                                                      // previousIndex
        previousPairHalfEdgeIndex,                                              // nextIndex
        previousNewHalfEdgeIndex1,                                              // pairIndex
        halfEdges.at(nextPairHalfEdgeIndex).getType()));                        // type

    // relabel original half-edges

    vertices[fromVertexIndex].setHalfEdgeIndex(newHalfEdgeIndex0);

    halfEdges[halfEdgeIndex].setFromIndex(newVertexIndex);
    halfEdges[halfEdgeIndex].setPreviousIndex(previousNewHalfEdgeIndex0);
    halfEdges[halfEdgeIndex].setPairIndex(newHalfEdgeIndex1);

    halfEdges[previousHalfEdgeIndex].setPreviousIndex(
        nextNewHalfEdgeIndex0);
    halfEdges[previousHalfEdgeIndex].setNextIndex(
        newHalfEdgeIndex0);

    halfEdges[nextHalfEdgeIndex].setNextIndex(previousNewHalfEdgeIndex0);

    vertices[toVertexIndex].setHalfEdgeIndex(newHalfEdgeIndex1);

    halfEdges[pairHalfEdgeIndex].setFromIndex(newVertexIndex);
    halfEdges[pairHalfEdgeIndex].setPreviousIndex(previousNewHalfEdgeIndex1);
    halfEdges[pairHalfEdgeIndex].setPairIndex(newHalfEdgeIndex0);

    halfEdges[previousPairHalfEdgeIndex].setPreviousIndex(
        nextNewHalfEdgeIndex1);
    halfEdges[previousPairHalfEdgeIndex].setNextIndex(
        newHalfEdgeIndex1);

    halfEdges[nextPairHalfEdgeIndex].setNextIndex(previousNewHalfEdgeIndex1);

    return std::make_tuple(newVertexIndex, createdHalfEdgeIndices);
}

std::tuple<long int, long int> Mesh::swapEdge(
    long int const& halfEdgeIndex,
    std::string const& type0, std::string const& type1) {

    long int const& previousHalfEdgeIndex =
        halfEdges.at(halfEdgeIndex).getPreviousIndex();
    long int const& nextHalfEdgeIndex =
        halfEdges.at(halfEdgeIndex).getNextIndex();
    long int const& pairHalfEdgeIndex =
        halfEdges.at(halfEdgeIndex).getPairIndex();
    long int const& previousPairHalfEdgeIndex =
        halfEdges.at(pairHalfEdgeIndex).getPreviousIndex();
    long int const& nextPairHalfEdgeIndex =
        halfEdges.at(pairHalfEdgeIndex).getNextIndex();

    // delete half-edges

    halfEdges.erase(halfEdgeIndex);
    halfEdges.erase(pairHalfEdgeIndex);

    // create half-edges

    halfEdges.emplace(halfEdgeIndex, HalfEdge(
        // halfEdges is std::map so add with std::map::emplace
        halfEdgeIndex,                                      // halfEdgeIndex
        halfEdges.at(nextPairHalfEdgeIndex).getToIndex(),   // fromIndex
        halfEdges.at(nextHalfEdgeIndex).getToIndex(),       // toIndex
        nextPairHalfEdgeIndex,                              // previousIndex
        previousHalfEdgeIndex,                              // nextIndex
        pairHalfEdgeIndex,                                  // pairIndex
        type0));                                            // type
    halfEdges.emplace(pairHalfEdgeIndex, HalfEdge(
        pairHalfEdgeIndex,                                  // halfEdgeIndex
        halfEdges.at(nextHalfEdgeIndex).getToIndex(),       // fromIndex
        halfEdges.at(nextPairHalfEdgeIndex).getToIndex(),   // toIndex
        nextHalfEdgeIndex,                                  // previousIndex
        previousPairHalfEdgeIndex,                          // nextIndex
        halfEdgeIndex,                                      // pairIndex
        type1));                                            // type

    // relabel half-edges

    vertices[halfEdges.at(nextHalfEdgeIndex).getFromIndex()]
        .setHalfEdgeIndex(nextHalfEdgeIndex);

    halfEdges[previousHalfEdgeIndex].setPreviousIndex(halfEdgeIndex);
    halfEdges[previousHalfEdgeIndex].setNextIndex(nextPairHalfEdgeIndex);

    halfEdges[nextHalfEdgeIndex].setPreviousIndex(previousPairHalfEdgeIndex);
    halfEdges[nextHalfEdgeIndex].setNextIndex(pairHalfEdgeIndex);

    vertices[halfEdges.at(nextPairHalfEdgeIndex).getFromIndex()]
        .setHalfEdgeIndex(nextPairHalfEdgeIndex);

    halfEdges[previousPairHalfEdgeIndex].setPreviousIndex(pairHalfEdgeIndex);
    halfEdges[previousPairHalfEdgeIndex].setNextIndex(nextHalfEdgeIndex);

    halfEdges[nextPairHalfEdgeIndex].setPreviousIndex(previousHalfEdgeIndex);
    halfEdges[nextPairHalfEdgeIndex].setNextIndex(halfEdgeIndex);

    return std::make_tuple(halfEdgeIndex, pairHalfEdgeIndex);
}

long int Mesh::changeToBoundary(
    long int const& vertexIndex) {

    // copy attributes

    std::vector<double> const position =
        vertices.at(vertexIndex).getPosition();
    std::vector<double> const uposition =
        vertices.at(vertexIndex).getUPosition();
    long int const halfEdgeIndex =
        vertices.at(vertexIndex).getHalfEdgeIndex();

    // delete vertex

    vertices.erase(vertexIndex);

    // create boundary vertex

    vertices.emplace(vertexIndex, Vertex(
        // vertices is std::map so add with std::map::emplace
        vertexIndex,
        position,
        halfEdgeIndex,
        true,   // boundary
        ""));   // type
    vertices[vertexIndex].setUPosition(uposition);

    return vertexIndex;
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

void Mesh::setSystemSize(std::vector<double> const& systemSize_) {

    std::vector<double> const disp =    // displacement of system centre
        {(systemSize_[0] - systemSize[0])/2,
        (systemSize_[1] - systemSize[1])/2};

    // move vertices
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        // wrapped position
        std::vector<double> const oldpos =
            (it->second).getPosition();
        std::vector<double> const newpos =
            {oldpos[0] + disp[0], oldpos[1] + disp[1]};
        if (!(it->second).getBoundary()) {
            assert(newpos[0] < systemSize_[0] && newpos[1] < systemSize_[1]);
        }
        (it->second).setPosition(newpos);
        // unwrapped position
        std::vector<double> const oldupos =
            (it->second).getUPosition();
        std::vector<double> const newupos =
            {oldupos[0] + disp[0], oldupos[1] + disp[1]};
        (it->second).setUPosition(newupos);
    }

    // change system size
    systemSize = systemSize_;
}

void Mesh::checkMesh(bool const& checkOrientations) const {

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

            if (!boundary && checkOrientations) {    // only for inner triangles
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

