#include <cmath>
#include <numbers>
#include <assert.h>

#include "system.hpp"
#include "tools.hpp"

void VertexModel::integrate(double const& dt,
    double const& delta, double const& epsilon) {

    // get forces

    computeForces();

    // integrate positions

    long int vertexIndex;
    std::vector<double> uposition;
    for (auto it=forces.begin(); it != forces.end(); ++it) {
        vertexIndex = it->first;
        uposition = vertices.at(vertexIndex).getUPosition();
        for (int dim=0; dim < 2; dim++) {
            uposition[dim] += forces[vertexIndex][dim]*dt;  // Euler integration of position
        }
        vertices[vertexIndex].setUPosition(                 // unwrapped position
            uposition);
        vertices[vertexIndex].setPosition(                  // (wrapped) position
            wrap(vertices[vertexIndex].getUPosition()));
    }

    // move cell centres

    std::vector<double> cellUPosition(0), initialCellPosition(0);
    std::vector<long int> neighbourVerticesIndices(0); int numberNeighbours;
    std::vector<double> disp(0);
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        if ((it->second).getType() != "centre") { continue; }   // only consider cell centres
        vertexIndex = (it->second).getIndex();

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

    // integrate internal degrees of freedom

    for (auto it=halfEdgeForces.begin(); it != halfEdgeForces.end(); ++it)
        { (it->second)->integrate(dt); }
    for (auto it=vertexForces.begin(); it != vertexForces.end(); ++it)
        { (it->second)->integrate(dt); }

    // update time

    time += dt;
}

void VertexModel::computeForces() {

    // clear forces
    forces.clear();
    for (auto it=vertices.begin(); it != vertices.end(); ++it)
        { forces[it->first] = {0, 0}; }

    // compute forces
    for (auto it=halfEdgeForces.begin(); it != halfEdgeForces.end(); ++it)
        { (it->second)->addAllForces(); }
    for (auto it=vertexForces.begin(); it != vertexForces.end(); ++it)
        { (it->second)->addAllForces(); }

    // remove centre of mass force
    if (true) {
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
    }
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
            if (vertices[vertexIndex].getType() == "centre"             // cell centres neighbours of the first vertex...
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
            if (vertices[vertexIndex].getType() == "centre"             // cell centres neighbours of the first vertex...
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

        deleteEdge(mergeHalfEdgeIndex);

        // create new vertex

        createEdge(createHalfEdgeIndex0, createHalfEdgeIndex1,
            angle, delta + epsilon,
            "junction", "");
    }
//     if (halfEdgeIndices.size() > 0) { checkMesh(); }
}

void VertexModel::checkMesh(std::vector<std::string> halfEdgeTypes) const {

    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
        if (inVec(halfEdgeTypes, (it->second).getType())) {                 // type of half-edge is in helfEdgeTypes
            if ((it->second).getType()
                == halfEdges.at((it->second).getPairIndex()).getType()) {   // half-edge and pair have identical types
                throw std::runtime_error("Pair half-edges have identical type "
                    + (it->second).getType() + ".");
            }
        }
    }

    Mesh::checkMesh();
}

