#include <cmath>
#include <numbers>

#include "assert.hpp"
#include "system.hpp"
#include "tools.hpp"

void VertexModel::integrate(double const& dt,
    double const& delta, double const& epsilon) {

    // integrate positions

    long int vertexIndex;
    std::vector<double> uposition;
    if (dt > 0) {
        for (auto it=velocities.begin(); it != velocities.end(); ++it) {
            vertexIndex = it->first;
            uposition = vertices.at(vertexIndex).getUPosition();
            for (int dim=0; dim < 2; dim++) {
                uposition[dim] += velocities[vertexIndex][dim]*dt;  // Euler integration of position
            }
            vertices[vertexIndex].setUPosition(                     // unwrapped position
                uposition);
            vertices[vertexIndex].setPosition(                      // (wrapped) position
                wrap(vertices[vertexIndex].getUPosition()));
        }
    }

    // perform T1s and update internal degrees of freedom

    if (dt > 0) { doT1(delta, epsilon); }

    // move cell centres

    #if true
    std::vector<long int> const cellCentreVertexIndices =
        getVertexIndicesByType("centre");
    for (long int vertexIndex : cellCentreVertexIndices)    // only consider cell centres
        { moveToNeighboursCentroid(vertexIndex); }
    #endif

    // integrate internal degrees of freedom

    for (auto it=halfEdgeForces.begin(); it != halfEdgeForces.end(); ++it)
        { (it->second)->integrate(dt); }
    for (auto it=vertexForces.begin(); it != vertexForces.end(); ++it)
        { (it->second)->integrate(dt); }

    // compute forces and integrate velocities

    integrateVelocities(dt);

    // update time

    time += dt;
}

void VertexModel::integrateVelocities(double const& dt) {

    // clear forces and velocities
    forces.clear();
    for (auto it=vertices.begin(); it != vertices.end(); ++it)
        { if (!(it->second).getBoundary()) { forces[it->first] = {0, 0}; } }

    // compute forces
    for (auto it=halfEdgeForces.begin(); it != halfEdgeForces.end(); ++it)
        { (it->second)->addAllForces(); }
    for (auto it=vertexForces.begin(); it != vertexForces.end(); ++it)
        { (it->second)->addAllForces(); }

    // remove centre of mass force
    #if true
    long int const numberVertices = forces.size();
    double avForce[2] = {0, 0};
    for (int dim=0; dim < 2; dim++) {
        for (auto it=forces.begin(); it != forces.end(); ++it) {
            avForce[dim] += forces[it->first][dim]/numberVertices;
        }
        for (auto it=forces.begin(); it != forces.end(); ++it) {
            forces[it->first][dim] -= avForce[dim];
        }
    }
    #endif

    // get velocities
    integrator->integrate(dt);
}

void VertexModel::doT1(double const& delta, double const& epsilon) {

    // identify small junctions

    std::vector<long int> halfEdgeIndices(0);
    long int fromMergeIndex, toMergeIndex;
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
        if ((it->second).getType() != "junction") { continue; } // loop over junctions
        fromMergeIndex =                                        // (first) vertex to be (potentially) merged into neighbour
            halfEdges.at((it->second).getIndex()).getFromIndex();
        toMergeIndex =                                          // (second) vertex towards which neighbour is (potentially) merged
            halfEdges.at((it->second).getIndex()).getToIndex();
        if (vertices.at(fromMergeIndex).getBoundary()
            || vertices.at(toMergeIndex).getBoundary()) {
            continue;               // boundary vertex cannot be merged
        }
        if (getEdgeLength((it->second).getIndex()) < delta) {
            halfEdgeIndices.push_back((it->second).getIndex());
        }
    }
    if (halfEdgeIndices.size() == 0) return;
    random.shuffle(halfEdgeIndices);    // do T1s in a random order

    // perform T1s

    std::vector<std::vector<long int>> neighbours;
    std::vector<long int> neighboursFromMerge, halfEdgesNeighboursFromMerge;
    std::vector<long int> neighboursToMerge, halfEdgesNeighboursToMerge;
    int numberNeighbours;
    std::vector<long int> halfEdgeToCellsIndices;
    std::vector<long int> halfEdgeToBoundariesIndices;
    long int vertexIndex;
    long int createHalfEdgeIndex0, createHalfEdgeIndex1;            // saving both half-edge index...
    long int createHalfEdgePairIndex0, createHalfEdgePairIndex1;    // ... and pair half-edge index in case either one is deleted in the process
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
        halfEdgeToBoundariesIndices.clear();
        for (int i=0; i < numberNeighbours; i++) {
            vertexIndex = neighboursFromMerge[i];
            if ((vertices.at(vertexIndex)).getType() == "centre" && // centre neighbours of the first vertex...
                !inVec(neighboursToMerge, vertexIndex)) {           // ... which are not neighbours of the second vertex
                halfEdgeToCellsIndices.push_back(
                    halfEdgesNeighboursFromMerge[i]);
            }
            if (vertices.at(vertexIndex).getBoundary()) {           // boundary neighbours of the first vertex

                halfEdgeToBoundariesIndices.push_back(
                    halfEdgesNeighboursFromMerge[i]);
            }
        }
        assert(halfEdgeToCellsIndices.size() > 0
            || halfEdgeToBoundariesIndices.size() > 0);
        if (halfEdgeToCellsIndices.size() > 0) {                    // randomly pick one preferably towards a cell
            createHalfEdgeIndex0 = random.pick(halfEdgeToCellsIndices);
        }
        else {
            createHalfEdgeIndex0 = random.pick(halfEdgeToBoundariesIndices);
        }
        createHalfEdgePairIndex0 =
            (halfEdges.at(createHalfEdgeIndex0)).getPairIndex();

        numberNeighbours = neighboursToMerge.size();
        assert(numberNeighbours == (int) halfEdgesNeighboursToMerge.size());
        halfEdgeToCellsIndices.clear();
        halfEdgeToBoundariesIndices.clear();
        for (int i=0; i < numberNeighbours; i++) {
            vertexIndex = neighboursToMerge[i];
            if ((vertices.at(vertexIndex)).getType() == "centre" && // centre neighbours of the second vertex...
                !inVec(neighboursFromMerge, vertexIndex)) {         // ... which are not neighbours of the first vertex
                halfEdgeToCellsIndices.push_back(
                    halfEdgesNeighboursToMerge[i]);
            }
            if (vertices.at(vertexIndex).getBoundary()) {           // boundary neighbours of the second vertex

                halfEdgeToBoundariesIndices.push_back(
                    halfEdgesNeighboursToMerge[i]);
            }
        }
        assert(halfEdgeToCellsIndices.size() > 0
            || halfEdgeToBoundariesIndices.size() > 0);
        if (halfEdgeToCellsIndices.size() > 0) {                    // randomly pick one preferably towards a cell
            createHalfEdgeIndex1 = random.pick(halfEdgeToCellsIndices);
        }
        else {
            createHalfEdgeIndex1 = random.pick(halfEdgeToBoundariesIndices);
        }
        createHalfEdgePairIndex1 =
            (halfEdges.at(createHalfEdgeIndex1)).getPairIndex();

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

        // delete edge

        TopoChangeEdgeInfoType const del =
            deleteEdge(mergeHalfEdgeIndex);

        for (auto it=halfEdgeForces.begin(); it != halfEdgeForces.end(); ++it)
            { (it->second)->deleteEdge(del); }
        for (auto it=vertexForces.begin(); it != vertexForces.end(); ++it)
            { (it->second)->deleteEdge(del); }

        // create edge

        TopoChangeEdgeInfoType const cre =
            createEdge(
                (!inVec(std::get<1>(del), createHalfEdgeIndex0)) ?  // was this half-edge NOT deleted?
                    createHalfEdgeIndex0 :
                    halfEdges[createHalfEdgePairIndex0].getPairIndex(),
                (!inVec(std::get<1>(del), createHalfEdgeIndex1)) ?  // was this half-edge NOT deleted?
                    createHalfEdgeIndex1 :
                    halfEdges[createHalfEdgePairIndex1].getPairIndex(),
                angle, delta + epsilon,
                "junction", "");

        for (auto it=halfEdgeForces.begin(); it != halfEdgeForces.end(); ++it)
            { (it->second)->createEdge(cre); }
        for (auto it=vertexForces.begin(); it != vertexForces.end(); ++it)
            { (it->second)->createEdge(cre); }

    }
//     checkMesh({"junction"});
}

long int VertexModel::splitCell(
    long int const& halfEdgeIndex0, long int const& halfEdgeIndex1) {

    long int const cellVertexIndex =    // index of cell centre vertex
        halfEdges.at(halfEdges.at(halfEdgeIndex0).getNextIndex()).getToIndex();
    assert(cellVertexIndex ==           // check that two half-edges belong to the boundary of the same cell
        halfEdges.at(halfEdges.at(halfEdgeIndex1).getNextIndex()).getToIndex()
        );
    assert(vertices.at(cellVertexIndex).getType() == "centre");

    TopoChangeEdgeInfoType info;

    // create new vertices

    std::vector<long int> halfEdgeToSplitIndices;
    std::vector<long int> const halfEdgeIndices =
        {halfEdgeIndex0, halfEdgeIndex1};
    for (long int halfEdgeIndex : halfEdgeIndices) {            // loop over half-edges in the middle of which to create vertices
        long int const newVertexIndex =
            std::get<0>(splitVertices(halfEdgeIndex));          // add new vertex in the middle of half-edge
        halfEdgeToSplitIndices.push_back(
            getHalfEdgeIndex(cellVertexIndex, newVertexIndex)); // half-edge from cell centre to the created vertex
    }

    // create new cell centre

    long int const newCellVertexIndex = std::get<0>(
        createEdge(                         // split the half-edges from the initial centre to the newly created vertices in order to create the new cell centre
            halfEdgeToSplitIndices.at(0), halfEdgeToSplitIndices.at(1),
            0, 0, "", ""));
    long int const halfEdgeToSwapIndex =    // half-edge from initial cell centre to the created cell centre
        getHalfEdgeIndex(cellVertexIndex, newCellVertexIndex);

    // flip edge to create cell junction

    swapEdge(halfEdgeToSwapIndex, "junction");      // swap edge and create junction
    moveToNeighboursCentroid(cellVertexIndex);      // move initial cell vertex to its barycentre
    moveToNeighboursCentroid(newCellVertexIndex);   // move new cell vertex to its barycentre

//     checkMesh({"junction"});
    return newCellVertexIndex;
}

long int VertexModel::splitCellAtMax(
    long int const& cellVertexIndex) {

    double maxDist = 0;
    std::vector<long int> const centreToBoundaryHalfEdges = // half-edges from cell centre to cell corners
        getNeighbourVertices(cellVertexIndex)[1];
    std::vector<long int> maxHalfEdgeIndices;               // pair of half-edges which maximises distance between their centres

    for (long int halfEdgeIndex0 : centreToBoundaryHalfEdges) {
        long int const boundaryHalfEdgeIndex0 =             // first given boundary half-edge
            halfEdges.at(halfEdgeIndex0).getNextIndex();

        for (long int halfEdgeIndex1 : centreToBoundaryHalfEdges) {
            long int const boundaryHalfEdgeIndex1 =         // second given boundary half-edge
                halfEdges.at(halfEdgeIndex1).getNextIndex();

            double const dist = norm2(wrapDiff(             // distance between the centres of first and second given boundary half-edges
                getHalfEdgeCentre(boundaryHalfEdgeIndex0),
                getHalfEdgeCentre(boundaryHalfEdgeIndex1)));
            if (dist > maxDist) {
                maxDist = dist;
                maxHalfEdgeIndices.clear();
                maxHalfEdgeIndices.push_back(boundaryHalfEdgeIndex0);
                maxHalfEdgeIndices.push_back(boundaryHalfEdgeIndex1);
            }
        }
    }

    assert((int) maxHalfEdgeIndices.size() == 2);
    return splitCell(maxHalfEdgeIndices.at(0), maxHalfEdgeIndices.at(1));
}

long int VertexModel::mergeCell(
    long int const& halfEdgeIndex) {

    long int const fromCellVertexIndex =    // centre vertex index of the cell which will be merged from
        halfEdges.at(halfEdges.at(halfEdgeIndex).getNextIndex())
            .getToIndex();
    long int const pairHalfEdgeIndex =      // pair half-edge index belonging to cell which will be merged to
        halfEdges.at(halfEdgeIndex).getPairIndex();
    long int const toCellVertexIndex =      // centre vertex index of the cell which will be merged to
        halfEdges.at(halfEdges.at(pairHalfEdgeIndex).getNextIndex())
            .getToIndex();

    assert(vertices.at(fromCellVertexIndex).getType() == "centre");
    assert(vertices.at(toCellVertexIndex).getType() == "centre");

    // flip cell junction

    swapEdge(halfEdgeIndex);

    // identify edges to remove

    long int const fromToHalfEdgeIndex =    // half-edge connecting two cell centres
        getHalfEdgeIndex(fromCellVertexIndex, toCellVertexIndex);

    long int const halfEdgeToDeleteIndex0 = // index of half-edge connecting the third vertex in the triangle which contains fromToHalfEdgeIndex to a neighbouring vertex to which it will be merged
        halfEdges.at(
            halfEdges.at(
                halfEdges.at(
                    halfEdgeIndex)
                    .getPreviousIndex())
                .getPairIndex())
            .getNextIndex();
    long int const neighbourCellVertexIndex0 =
        halfEdges.at(
            halfEdges.at(
                halfEdges.at(
                    halfEdgeToDeleteIndex0)
                    .getPairIndex())
                .getNextIndex())
            .getToIndex();

    long int const halfEdgeToDeleteIndex1 = // index of half-edge connecting the third vertex in the triangle which contains the pair of fromToHalfEdgeIndex to a neighbouring vertex to which it will be merged
        halfEdges.at(
            halfEdges.at(
                halfEdges.at(
                    halfEdges.at(halfEdgeIndex).getPairIndex())
                    .getPreviousIndex())
                .getPairIndex())
            .getNextIndex();
    long int const neighbourCellVertexIndex1 =
        halfEdges.at(
            halfEdges.at(
                halfEdges.at(
                    halfEdgeToDeleteIndex1)
                    .getPairIndex())
                .getNextIndex())
            .getToIndex();

    // connect both cell centre

    deleteEdge(fromToHalfEdgeIndex);

    // merge vertices

    mergeVertices(halfEdgeToDeleteIndex0);
    mergeVertices(halfEdgeToDeleteIndex1);

    // move cells with merged vertices and destination cell centres
    if (vertices.at(neighbourCellVertexIndex0).getType() == "centre")
        moveToNeighboursCentroid(neighbourCellVertexIndex0);    // first cell with merged vertices centre
    if (vertices.at(neighbourCellVertexIndex1).getType() == "centre")
        moveToNeighboursCentroid(neighbourCellVertexIndex1);    // second cell with merged vertices centre
    moveToNeighboursCentroid(toCellVertexIndex);                // destination cell centre vertex


//     checkMesh({"junction"});
    return toCellVertexIndex;
}

void VertexModel::checkMesh(
    std::vector<std::string> const& halfEdgeTypes,
    bool const& checkOrientations) const {

    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
        if (inVec(halfEdgeTypes, (it->second).getType())) {                 // type of half-edge is in helfEdgeTypes
            if ((it->second).getType()
                == halfEdges.at((it->second).getPairIndex()).getType()) {   // half-edge and pair have identical types
                throw std::runtime_error("Pair half-edges have identical type "
                    + (it->second).getType() + ".");
            }
        }
    }

    Mesh::checkMesh(checkOrientations);
}

