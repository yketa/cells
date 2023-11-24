/*
Functions using VertexModel to speed up plotting.
*/

#include <vector>

#include "system.hpp"

/*
*/

std::vector<std::vector<double>> getLinesHalfEdge(
    VertexModel& vm) {
/*
Return vector [[x0, x0'], [y0, y0'], ..., [xN-1, xN-1'], [yN-1, yN-1']] where
the line (xi, yi) -- (xi', yi') corresponds to the i-th half-edge in `vm'.
*/

    std::vector<std::vector<double>> lines(0);

    std::map<long int, Vertex> vertices = vm.getVertices();
    std::map<long int, HalfEdge> halfEdges = vm.getHalfEdges();

    std::vector<double> position;
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {    // loop over all vertices
        position = (it->second).getPosition();
        (it->second).setPosition(vm.wrap(position));    // wrap position with respect to periodic boundary conditions
    }

    std::vector<double> fromPos, disp;
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {  // loop over all half-edges
        fromPos = vertices[(it->second).getFromIndex()].getPosition();  // position of origin vertex
        disp = vm.getHalfEdgeVector(it->first, false);                  // displacement to destination vertex
        lines.push_back({fromPos[0], fromPos[0] + disp[0]});            // x-coordinates of line
        lines.push_back({fromPos[1], fromPos[1] + disp[1]});            // y-coordinates of line
    }

    return lines;
}

std::vector<std::vector<double>> getLinesJunction(
    VertexModel& vm) {
/*
Return vector [[x0, x0'], [y0, y0'], ..., [xN-1, xN-1'], [yN-1, yN-1']] where
the line (xi, yi) -- (xi', yi') corresponds to the i-th junction in `vm'.
*/

    std::vector<std::vector<double>> lines(0);

    std::map<long int, Vertex> vertices = vm.getVertices();
    std::map<long int, HalfEdge> halfEdges = vm.getHalfEdges();
    MultiIntKeyDict<Junction> junctions = vm.getJunctions();

    std::vector<double> position;
    for (auto it=vertices.begin(); it != vertices.end(); ++it) {    // loop over all vertices
        position = (it->second).getPosition();
        (it->second).setPosition(vm.wrap(position));    // wrap position with respect to periodic boundary conditions
    }

    std::vector<double> fromPos, disp;
    std::vector<long int> halfEdgeIndices;
    for (auto it=junctions.begin(); it != junctions.end(); ++it) {  // loop over all junctions
        halfEdgeIndices.clear();
        halfEdgeIndices.push_back(
            *it);
        halfEdgeIndices.push_back(
            halfEdges[halfEdgeIndices[0]].getPairIndex());
        for (long int index : halfEdgeIndices) {
            fromPos = vertices[halfEdges[index].getFromIndex()].getPosition();  // position of origin vertex
            disp = vm.getHalfEdgeVector(index, false);                          // displacement to destination vertex
            lines.push_back({fromPos[0], fromPos[0] + disp[0]});                // x-coordinates of line
            lines.push_back({fromPos[1], fromPos[1] + disp[1]});                // y-coordinates of line
        }
    }

    return lines;
}

std::vector<std::vector<std::vector<double>>> getPolygonsCell(
    VertexModel& vm) {
/*
Return vector [..., [[xi^0, yi^0], ..., [xi^Ni-1, yi^Ni-1]], ...] where the
point (xi^j, yi^j) is the j-th corner of the i-th cell.
*/

    std::vector<std::vector<std::vector<double>>> polygons(0);

    std::map<long int, Vertex> vertices = vm.getVertices();
    MultiIntKeyDict<Cell> cells = vm.getCells();

    long int n;
    long int vertexIndex;
    std::vector<long int> halfEdgesToNeighboursIndices;
    std::vector<double> cellPos, disp;
    for (const Cell* cell : cells.getValues()) {    // loop over all cells
        n = polygons.size();
        polygons.push_back(std::vector<std::vector<double>>(0));

        vertexIndex = cell->getVertexIndex();
        cellPos = vertices[vertexIndex].getPosition();          // position of cell centre
        halfEdgesToNeighboursIndices = vm.getNeighbourVertices(vertexIndex)[1];
        for (long int halfEdgeIndex : halfEdgesToNeighboursIndices) {   // loop over cell vertices
            disp = vm.getHalfEdgeVector(halfEdgeIndex, false);  // displacement to vertex
            polygons[n].push_back(                              // vertex position
                {cellPos[0] + disp[0], cellPos[1] + disp[1]});
        }
    }

    return polygons;
}

