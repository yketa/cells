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
        fromMergeIndex =            // (first) vertex to be (potentially) merged into neighbour
            halfEdges.at(junction->getHalfEdgeIndex()).getFromIndex();
        toMergeIndex =              // (second) vertex towards which neighbour is (potentially) merged
            halfEdges.at(junction->getHalfEdgeIndex()).getToIndex();
        if (vertices.at(fromMergeIndex).getBoundary()
            || vertices.at(toMergeIndex).getBoundary()) {
            continue;               // boundary vertex cannot be merged
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
            "junction", "junctionPair");
    }
//     if (halfEdgeIndices.size() > 0) { checkMesh(); }
}

