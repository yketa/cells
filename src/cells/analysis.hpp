/*
Functions for vertex model analysis.
*/

#include <complex>
#include <map>
#include <vector>

using namespace std::complex_literals;  // using comlex numbers

#include "system.hpp"
#include "tools.hpp"

std::map<long int, std::vector<double>>
getCentreVelocities(VertexModel const& vm) {
/*
Compute velocities of the centroid of each cell.
*/

    std::map<long int, std::vector<double>> centreVelocities;

    std::map<long int, std::vector<double>> const velocities =
        vm.getVelocities();
    std::map<long int, Vertex> const vertices =
        vm.getVertices();

    std::vector<long int> const vertexIndices =
        vm.getVertexIndicesByType("centre");
    for (long int vertexIndex : vertexIndices) {

        std::vector<double> const upos =                            // uwrapped position of cell
            vertices.at(vertexIndex).getUPosition();
        double const cellArea =                                     // cell area
            vm.getVertexToNeighboursArea(vertexIndex);

        std::vector<long int> const neighbours =
            vm.getNeighbourVertices(vertexIndex)[0];
        long int const nNeighbours = neighbours.size(); // number of neighbours
        centreVelocities[vertexIndex] = {0, 0};
        for (long int i=0; i < nNeighbours; i++) {      // loop over neighbours

            long int const index0 =
                neighbours[i];
            std::vector<double> const upos0 = vm.wrapDiff(upos, // unwrapped position of neighbour
                vertices.at(index0).getUPosition());
            long int const index1 =
                neighbours[pmod(i + 1, nNeighbours)];
            std::vector<double> const upos1 = vm.wrapDiff(upos, // unwrapped position of next neighbour
                vertices.at(index1).getUPosition());

            for (int dim=0; dim < 2; dim++) {
                centreVelocities[vertexIndex][dim] +=
                    ((velocities.at(index0).at(dim)
                        + velocities.at(index1).at(dim))*(
                        upos0[0]*upos1[1] - upos1[0]*upos0[1])
                    + (upos0[dim] + upos1[dim])*(
                        velocities.at(index0).at(0)*upos1[1]
                        + upos0[0]*velocities.at(index1).at(1)
                        - velocities.at(index1).at(0)*upos0[1]
                        - upos1[0]*velocities.at(index0).at(1))
                    )/(6*cellArea);
            }
        }
    }

    return centreVelocities;
}

std::map<long int, std::vector<std::vector<double>>>
getVectorsToNeighbouringCells(VertexModel const& vm) {
/*
Compute vectors to neighbouring cell centres.
*/

    std::map<long int, std::vector<std::vector<double>>>
        vectorsToNeighbours;

    std::map<long int, Vertex> const vertices =
        vm.getVertices();
    std::map<long int, HalfEdge> const halfEdges =
        vm.getHalfEdges();

    std::vector<long int> const vertexIndices =
        vm.getVertexIndicesByType("centre");
    for (long int vertexIndex : vertexIndices) {
        vectorsToNeighbours.emplace(vertexIndex,
            std::vector<std::vector<double>>(0));

        std::vector<long int> const halfEdgesToNeighbours =
            vm.getNeighbourVertices(vertexIndex)[1];
        for (long int index : halfEdgesToNeighbours) {  // loop over half-edges to neighbour vertices

            long int const neighbourVertexIndex =       // neighbour cell index
                halfEdges.at(
                    halfEdges.at(
                        halfEdges.at(
                            halfEdges.at(index).getNextIndex()
                        ).getPairIndex()
                    ).getNextIndex()
                ).getToIndex();

            if (vertices.at(neighbourVertexIndex).getBoundary())    // ignore boundary vertices
                { continue; }
            assert(                                                 // the vertex should be a cell centre
                vertices.at(neighbourVertexIndex).getType()
                    == "centre");
            vectorsToNeighbours[vertexIndex].push_back(
                vm.wrapTo(vertexIndex, neighbourVertexIndex,        // unnormalised vector to neighbouring cell
                    false));
        }
    }

    return vectorsToNeighbours;
}

std::map<long int, std::complex<double>>
getPAticOrderParameters(VertexModel const& vm, int const& p) {
/*
Compute p-atic order parameters of cell centres.
*/

    std::map<long int, std::complex<double>> psip;                              // p-atic order parameter
    std::map<long int, std::vector<std::vector<double>>> vectorsToNeighbours =  // vectors to neighbouring cell centres
        getVectorsToNeighbouringCells(vm);

    for (auto it=vectorsToNeighbours.begin(); it != vectorsToNeighbours.end();
        ++it) {                                                     // loop over cells

        psip.emplace(it->first, 0);
        std::complex<double> const nNeigh = (it->second).size();
        for (std::vector<double> vectorToNeighbour : it->second) {  // loop over neighbours
            psip[it->first] += exp(1i*(p*angle2(vectorToNeighbour)))/nNeigh;
        }
    }

    return psip;
}

std::map<long int, double> getMaxLengthCells(VertexModel const& vm) {
/*
Return maximum length between two cell corners in each cell, considering
periodic boundary conditions.
*/

    std::vector<long int> const centres =
        vm.getVertexIndicesByType("centre");
    std::map<long int, Vertex> const vertices =
        vm.getVertices();

    std::map<long int, double> maxLength;
    for (long int index : centres) {
        std::vector<long int> const neighbours =
            vm.getNeighbourVertices(index)[0];
        long int const nNeighbours = neighbours.size();
        double max = 0;
        for (long int mu=0; mu < nNeighbours; mu++) {
            for (long int nu=mu + 1; nu < nNeighbours; nu++) {
                std::vector<double> const diff =    // considers periodic boundary conditions
                    vm.wrapTo(neighbours.at(mu), neighbours.at(nu));
                max = std::max(max,
                    sqrt(diff[0]*diff[0] + diff[1]*diff[1]));
            }
        }
        maxLength.emplace(index, max);
    }

    return maxLength;
}

std::map<long int, double> getMaxLengthBoundaries(VertexModel const& vm) {
/*
Return maximum length between two cell corners around each boundary vertex,
ignoring periodic boundary conditions.
*/

    std::map<long int, double> maxLength;
    std::map<long int, Vertex> const vertices = vm.getVertices();
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {
        if (!(it->second).getBoundary()) continue;
        std::vector<long int> const neighbours =
            vm.getNeighbourVertices(it->first)[0];
        long int const nNeighbours = neighbours.size();
        double max = 0;
        for (long int mu=0; mu < nNeighbours; mu++) {
            for (long int nu=mu + 1; nu < nNeighbours; nu++) {
                std::vector<double> const diff = {  // ignores periodic boundary conditions
                    vertices.at(neighbours.at(nu))
                        .getPosition()[0]
                        - vertices.at(neighbours.at(mu))
                            .getPosition()[0],
                    vertices.at(neighbours.at(nu))
                        .getPosition()[1]
                        - vertices.at(neighbours.at(mu))
                            .getPosition()[1]};
                max = std::max(max,
                    sqrt(diff[0]*diff[0] + diff[1]*diff[1]));
            }
        }
        maxLength.emplace(it->first, max);
    }

    return maxLength;
}

std::map<long int, double> getPercentageKeptNeighbours(
    VertexModel const& vm0, VertexModel const& vm1, double const& a) {
/*
Return percentage of kept cell neighbours between two states of a vertex model.
*/

    std::map<long int, long int> nNeigh;    // total number of neighbours
    std::map<long int, long int> nKept;     // number of kept neighbours

    std::vector<long int> const centres0 =
        vm0.getVertexIndicesByType("centre");
    for (long int index0 : centres0)
        { nNeigh.emplace(index0, 0); nKept.emplace(index0, 0); }

    std::map<long int, HalfEdge> const halfEdges0 =
        vm0.getHalfEdges();
    std::map<long int, Vertex> const vertices0 =
        vm0.getVertices();

    std::map<long int, HalfEdge> const halfEdges1 =
        vm1.getHalfEdges();
    std::map<long int, Vertex> const vertices1 =
        vm1.getVertices();

    std::vector<long int> const junctions0 =
        vm0.getHalfEdgeIndicesByType("junction");
    for (long int index0 : junctions0) {

        long int const cellA0 =                             // first cell of pair
            halfEdges0.at(
                halfEdges0.at(index0).getNextIndex()
            ).getToIndex();
        long int const cellB0 =                             // second cell of pair
            halfEdges0.at(
                halfEdges0.at(
                    halfEdges0.at(index0).getPairIndex()
                ).getNextIndex()
            ).getToIndex();
        std::vector<double> const initSep =         // initial separation
            vm0.wrapTo(cellA0, cellB0);

        if (norm2(initSep) > a) { continue; }       // too far for neighbours
        nNeigh[cellA0]++; nNeigh[cellB0]++;

        std::vector<long int> const neighbours1 =           // half-edges to neighbours
            vm1.getNeighbourVertices(cellA0)[1];
        for (long int index1 : neighbours1) {
            long int const cellB1 =                         // tentative second cell of pair
                halfEdges1.at(
                    halfEdges1.at(
                        halfEdges1.at(
                            halfEdges1.at(index1).getNextIndex()
                        ).getPairIndex()
                    ).getNextIndex()
                ).getToIndex();
            if (cellB1 == cellB0) {                 // these are the same cells...

                std::vector<double> const uposA0 =  // initial position of first cell
                    vertices0.at(cellA0).getUPosition();
                std::vector<double> const uposA1 =  // final position of first cell
                    vertices1.at(cellA0).getUPosition();

                std::vector<double> const uposB0 =  // initial position of second cell
                    vertices0.at(cellB0).getUPosition();
                std::vector<double> const uposB1 =  // final position of second cell
                    vertices1.at(cellB0).getUPosition();

                std::vector<double> const finSep =  // final separation
                    {initSep[0]
                        + (uposB1[0] - uposB0[0])
                        - (uposA1[0] - uposA0[0]),
                    initSep[1]
                        + (uposB1[1] - uposB0[1])
                        - (uposA1[1] - uposA0[1])};

                if (!(norm2(finSep) > a))           // ... which are still neighbours
                    { nKept[cellA0]++; nKept[cellB0]++; break; }
            }
        }
    }

    std::map<long int, double> pct;
    for (long int index0 : centres0) {
        pct.emplace(index0,
            nNeigh[index0] == 0 ? 0 :
            (double) nKept[index0]/ (double) nNeigh[index0]);
    }
    return pct;
}

