/*
Functions using VertexModel to speed up plotting.
*/

#ifndef PLOT_HPP
#define PLOT_HPP

#include <cmath>
#include <string>
#include <vector>

#include "system.hpp"

/*
*/

std::vector<std::vector<std::vector<double>>> getLinesHalfEdge(
    VertexModel const& vm, std::vector<long int> const& indices) {
/*
Return vector [[x0, x0'], [y0, y0'], ..., [xN-1, xN-1'], [yN-1, yN-1']] where
the line (xi, yi) -- (xi', yi') corresponds to the half-edge in `vm' whose
index is the i-th element of `indices'.
*/

    std::vector<std::vector<std::vector<double>>> lines(0);

    VerticesType const vertices = vm.getVertices();
    HalfEdgesType const halfEdges = vm.getHalfEdges();

    std::vector<double> fromPos, disp;
    for (long int halfEdgeIndex : indices) {                    // loop over all half-edges
        fromPos =                                               // position of origin vertex
            (vertices.at(
                (halfEdges.at(halfEdgeIndex)).getFromIndex())
            ).getPosition();
        disp = vm.getHalfEdgeVector(halfEdgeIndex, false);      // displacement to destination vertex
        lines.push_back({
            {fromPos[0], fromPos[1]},
            {fromPos[0] + disp[0], fromPos[1] + disp[1]}});
    }

    return lines;
}

std::vector<std::vector<std::vector<double>>> getLinesJunction(
    VertexModel const& vm) {
/*
Return vector [[[x0, y0], [x0', y0']], ..., [[xN-1, yN-1], [xN-1', yN-1']]]
where the line (xi, yi) -- (xi', yi') corresponds to the i-th junction in `vm'.
*/

    std::vector<std::vector<std::vector<double>>> lines(0);

    VerticesType const vertices = vm.getVertices();
    HalfEdgesType const halfEdges = vm.getHalfEdges();

    std::vector<double> fromPos, disp;
    std::vector<long int> const halfEdgeIndices =
        vm.getHalfEdgeIndicesByType("junction");
    for (long int halfEdgeIndex : halfEdgeIndices) {                            // loop over all junctions

        for (long int index :
            {halfEdgeIndex, (halfEdges.at(halfEdgeIndex)).getPairIndex()}) {    // junction half-edge and its pair half-edge

            fromPos = (vertices.at((halfEdges.at(index)).getFromIndex()))
                .getPosition();                                     // position of origin vertex
            disp = vm.getHalfEdgeVector(index, false);              // displacement to destination vertex
            lines.push_back({
                {fromPos[0], fromPos[1]},
                {fromPos[0] + disp[0], fromPos[1] + disp[1]}});
        }
    }

    return lines;
}

std::vector<std::vector<std::vector<double>>> getPolygonsCell(
    VertexModel const& vm) {
/*
Return vector [..., [[xi^0, yi^0], ..., [xi^Ni-1, yi^Ni-1]], ...] where the
point (xi^j, yi^j) is the j-th corner of the i-th cell.
*/

    std::vector<std::vector<std::vector<double>>> polygons(0);

    VerticesType const vertices = vm.getVertices();

    long int n;
    std::vector<long int> halfEdgesToNeighboursIndices;
    std::vector<double> cellPos, disp;
    std::vector<long int> const vertexIndices =
        vm.getVertexIndicesByType("centre");
    for (long int vertexIndex : vertexIndices) {    // loop over all cells
        n = polygons.size();
        polygons.push_back(std::vector<std::vector<double>>(0));

        cellPos = (vertices.at(vertexIndex)).getPosition();             // position of cell centre
        halfEdgesToNeighboursIndices = vm.getNeighbourVertices(vertexIndex)[1];
        for (long int halfEdgeIndex : halfEdgesToNeighboursIndices) {   // loop over cell vertices
            disp = vm.getHalfEdgeVector(halfEdgeIndex, false);          // displacement to vertex
            polygons[n].push_back(                                      // vertex position
                {cellPos[0] + disp[0], cellPos[1] + disp[1]});
        }
    }

    return polygons;
}

std::vector<std::vector<double>> getPrincipalAxesCell(
    VertexModel const& vm) {
/*
Return vector [[x0, y0], ..., [xN-1, yN-1]] where [xi, yi] is the unit vector
corresponding to the principal axis (eigenvector of the inertia tensor
corresponding to the largest eigenvalue) of the i-th cell.
*/

    std::vector<std::vector<double>> principalAxes(0);

    std::vector<long int> const cellIndices =
        vm.getVertexIndicesByType("centre");
    for (long int cellIndex : cellIndices) {    // loop over all cells

        std::vector<std::vector<double>> I = {{0, 0}, {0, 0}};  // inertia tensor

        std::vector<long int> const neighbourIndices =          // cell corners = neighbours of centre
            vm.getNeighbourVertices(cellIndex)[0];
        for (long int neighbourIndex : neighbourIndices) {      // loop over cell corners
            std::vector<double> const r =                       // vector from centre to corner
                vm.wrapTo(cellIndex, neighbourIndex, false);
            for (int alpha=0; alpha < 2; alpha++) {
                for (int beta=0; beta < 2; beta++) {
                    I[alpha][beta] +=
                        ((alpha == beta) ? (r[0]*r[0] + r[1]*r[1]) : 0)
                        - r[alpha]*r[beta];
                }
            }
        }

        std::vector<double> const evec = {                              // (manual) eigenvector with largest eigenvalue
            I[0][0] - I[1][1] - sqrt(
                I[0][0]*I[0][0] + 4*I[0][1]*I[0][1]
                - 2*I[0][0]*I[1][1] + I[1][1]*I[1][1]),
            2*I[0][1]};
        double const evecn = sqrt(evec[0]*evec[0] + evec[1]*evec[1]);   // eigenvector norm
        principalAxes.push_back({evec[0]/evecn, evec[1]/evecn});        // unitary eigenvector
    }

    return principalAxes;
}

#endif

