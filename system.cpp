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
        assert(numberNeighbours == halfEdgesNeighboursFromMerge.size());
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
        assert(numberNeighbours == halfEdgesNeighboursToMerge.size());
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
}

