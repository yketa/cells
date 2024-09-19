/*
Functions for vertex model analysis.
*/

#include <complex>
#include <map>
#include <numbers>
#include <vector>

using namespace std::complex_literals;  // using comlex numbers

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "system.hpp"
#include "tools.hpp"

/*
 *  Vertex model
 *
 */

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

/*
 *  General particles
 *
 */

std::vector<double> _get2DL(pybind11::array_t<double> const& L) {
/*
Return (2,) box size dimensions from array `L' which may be a double or an
array of double.
*/

    auto l = L.unchecked<>();
    if (l.ndim() > 1) {                             // system size is 2D or higher-D: impossible
        throw std::invalid_argument("L must be a 0- or 1-dimensional.");
    }
    else if (l.size() == 0) {                       // no system size: impossible
        throw std::invalid_argument("System size cannot be empty.");
    }
    else if (l.size() == 1) {                       // system size is the same for both dimensions
        if (l.ndim() == 0) {
            return std::vector<double>({l(), l()});
        }
        else {  // l.ndim() == 1
            return std::vector<double>({l(0), l(0)});
        }
    }
    else {                                          // system size is different for the two dimensions
        return std::vector<double>({l(0), l(1)});
    }
}

pybind11::array_t<double> getAllWaveVectors2D(
    pybind11::array_t<double> const& L,
    double const& qmin, double const& qmax) {
/*
Return wave vectors associated to rectangular box of size `L' such that their
norms belong to [`qmin', `qmax'].
Only a single vector of each pair of opposite wave vectors is returned.
*/

    // check system size
    std::vector<double> const systemSize = _get2DL(L);
    double const qx = 2*std::numbers::pi/systemSize[0];
    double const qy = 2*std::numbers::pi/systemSize[1];

    // loop in wave vector space
    std::vector<std::vector<double>> waveVectors(0);
    long long int const xmin = 0;
    long long int const xmax = floor(qmax/qx);
    for (long long int x=xmin; x <= xmax; x++) {
        long long int const ymin =
            std::max(
                0.,
                ceil(sqrt(pow(qmin, 2) - pow(qx*x, 2))/qy));
        long long int const ymax =
            std::min(
                floor(qmax/qy),
                floor(sqrt(pow(qmax, 2) - pow(qx*x, 2))/qy));
        for (long long int y=ymin; y <= ymax; y++) {
            double const qq = sqrt(pow(qx*x, 2) + pow(qy*y, 2));
            if (qq < qmin || qq > qmax) { continue; }   // wave vector norm not within interval
            waveVectors.push_back({qx*x, qy*y});
            if (x != 0 && y != 0) {
                // if x == 0 then (qx x, qy y) and (-qx x, qy y) are identical
                // if y == 0 then (qx x, qy y) and (-qx x, qy y) are opposite
                waveVectors.push_back({-qx*x, qy*y});
            }
        }
    }

    // create and return array
    if (waveVectors.size() == 0) { return pybind11::array_t<double>(); }
    pybind11::array_t<double>
        arr(std::vector<ptrdiff_t>{(long long int) waveVectors.size(), 2});
    auto a = arr.mutable_unchecked<2>();
    for (long long int l=0; l < (long long int) waveVectors.size(); l++) {
        for (int dim=0; dim < 2; dim++) {
            a(l, dim) = waveVectors[l][dim];
        }
    }
    return arr;
}

pybind11::array_t<std::complex<double>> getAllFT2D(
    pybind11::array_t<double> const& positions,
    pybind11::array_t<double> const& L,
    pybind11::array_t<std::complex<double>> const& values,
    double const& qmin, double const& qmax) {
/*
Return 2D Fourier transform of `values' associated to 2D `positions' at 2D wave
vectorswhose norms belong to [`qmin', `qmax'].
*/

    // wave vector norms
    pybind11::array_t<double> qARR = getAllWaveVectors2D(L, qmin, qmax);
    if (qARR.request().size == 0)
        { return pybind11::array_t<std::complex<double>>(0); }
    auto _q = qARR.unchecked<2>();      // direct access to wave vectors
    long int const n = _q.shape(0);

    // check positions and values arrays
    auto r = positions.unchecked<2>();  // direct access to positions
    assert(r.ndim() == 2);
    assert(r.shape(1) == 2);
    long int const N = r.shape(0);
    auto v = values.unchecked<>();      // direct access to first values
    if (v.shape(0) != N) {
        throw std::invalid_argument("Positions and values must have = sizes.");
    }
    if (v.ndim() > 1) {
        throw std::invalid_argument("Values must be 1-dimensional.");
    }

    // compute Fourier transform
    pybind11::array_t<std::complex<double>> ft({n});
    auto FT = ft.mutable_unchecked<1>();
    for (long int l=0; l < n; l++) {        // loop over wave vectors
        FT(l) = 0;                                                              // initialise
        for (long int i=0; i < N; i++) {    // loop over particles
            FT(l) += v(i)*std::exp(-1i*(_q(l, 0)*r(i, 0) + _q(l, 1)*r(i, 1)));  // Fourier term
        }
    }

    return ft;
}

