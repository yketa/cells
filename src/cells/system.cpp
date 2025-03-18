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
    std::vector<double> upos(0);
    std::vector<long int> neighbourVerticesIndices(0);
    long int numberNeighbours;
    std::vector<double> disp(0);
    std::vector<long int> const cellCentreVertexIndices =
        getVertexIndicesByType("centre");
    for (long int vertexIndex : cellCentreVertexIndices) {  // only consider cell centres

        double const cellArea = getVertexToNeighboursArea(vertexIndex);

        upos = {0, 0};
        neighbourVerticesIndices = getNeighbourVertices(vertexIndex)[0];
        numberNeighbours = neighbourVerticesIndices.size();
        std::vector<double> const uposr =           // use cell centre as reference rather than vertex
            vertices.at(vertexIndex).getUPosition();
        for (long int i=0; i < numberNeighbours; i++) {
            long int const index0 =
                neighbourVerticesIndices[i];
            std::vector<double> const upos0 =       // from reference position
                wrapDiff(uposr, vertices.at(index0).getUPosition());
            long int const index1 =
                neighbourVerticesIndices[pmod(i + 1, numberNeighbours)];
            std::vector<double> const upos1 =       // from reference position
                wrapDiff(uposr, vertices.at(index1).getUPosition());
            for (int dim=0; dim < 2; dim++) {
                upos[dim] +=
                    (upos0[dim] + upos1[dim])*(
                        upos0[0]*upos1[1] - upos1[0]*upos0[1]
                    )/(6*cellArea);
            }
        }
        for (int dim=0; dim < 2; dim++)             // add reference position
            { upos[dim] += uposr[dim]; }
        vertices[vertexIndex].setUPosition( // unwrapped position
            upos);
        vertices[vertexIndex].setPosition(  // (wrapped) position
            wrap(vertices.at(vertexIndex).getUPosition()));
    }
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
    moveToNeigboursBarycentre(cellVertexIndex);     // move initial cell vertex to its barycentre
    moveToNeigboursBarycentre(newCellVertexIndex);  // move new cell vertex to its barycentre

    checkMesh({"junction"});
    return newCellVertexIndex;
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

