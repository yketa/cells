/*
Bind C++ objects to python using pybind11.
https://pybind11.readthedocs.io/en/stable/index.html
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
#include <pybind11/stl.h>   // https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html

#include <string>
#include <vector>
#include <map>

#include "mesh.hpp"
#include "system.hpp"
#include "tools.hpp"

std::vector<std::vector<double>> getLinesHalfEdge(VertexModel& vm) {
/*
Return vector [[x0, x0'], [y0, y0'], ..., [xN-1, xN-1'], [yN-1, yN-1']] where
the line (xi, yi) -- (xi', yi') corresponds to the i-th half-edge in `vm'.
This is aimed as speeding plotting.
*/

    std::vector<std::vector<double>> lines(0);

    std::map<long int, Vertex> vertices = vm.getVertices();
    std::map<long int, HalfEdge> halfEdges = vm.getHalfEdges();

    for (auto it=vertices.begin(); it != vertices.end(); ++it) {    // loop over all vertices
        vm.wrap((it->second).getPosition());    // wrap position with respect to periodic boundary conditions
    }

    double* fromPos;
    std::vector<double> disp;
    for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {  // loop over all half-edges
        fromPos = vertices[*(it->second).getFromIndex()].getPosition(); // position of origin vertex
        disp = vm.getHalfEdgeVector(it->first, false);                  // displacement to destination vertex
        lines.push_back({fromPos[0], fromPos[0] + disp[0]});            // x-coordinates of line
        lines.push_back({fromPos[1], fromPos[1] + disp[1]});            // y-coordinates of line
    }

    return lines;
}

std::vector<std::vector<double>> getLinesJunction(VertexModel& vm) {
/*
Return vector [[x0, x0'], [y0, y0'], ..., [xN-1, xN-1'], [yN-1, yN-1']] where
the line (xi, yi) -- (xi', yi') corresponds to the i-th junction in `vm'.
This is aimed as speeding plotting.
*/

    std::vector<std::vector<double>> lines(0);

    std::map<long int, Vertex> vertices = vm.getVertices();
    std::map<long int, HalfEdge> halfEdges = vm.getHalfEdges();
    MultiIntKeyDict<Junction> junctions = vm.getJunctions();

    for (auto it=vertices.begin(); it != vertices.end(); ++it) {    // loop over all vertices
        vm.wrap((it->second).getPosition());    // wrap position with respect to periodic boundary conditions
    }

    double* fromPos;
    std::vector<double> disp;
    std::vector<long int> halfEdgeIndices;
    for (auto it=junctions.begin(); it != junctions.end(); ++it) {  // loop over all junctions
        halfEdgeIndices.clear();
        halfEdgeIndices.push_back(
            *it);
        halfEdgeIndices.push_back(
            *halfEdges[halfEdgeIndices[0]].getPairIndex());
        for (long int index : halfEdgeIndices) {
            fromPos = vertices[*halfEdges[index].getFromIndex()].getPosition(); // position of origin vertex
            disp = vm.getHalfEdgeVector(index, false);                          // displacement to destination vertex
            lines.push_back({fromPos[0], fromPos[0] + disp[0]});                // x-coordinates of line
            lines.push_back({fromPos[1], fromPos[1] + disp[1]});                // y-coordinates of line
        }
    }

    return lines;
}

/*
declare class template
https://stackoverflow.com/questions/66227840
https://stackoverflow.com/questions/47487888
iterating through object
https://stackoverflow.com/questions/49257947
*/

template<class T> void declare_MultiIntKeyDict_class(
    pybind11::module& m, std::string const& nameT) {

    std::string const type = "MultiIntKeyDict<" + nameT + ">";
    pybind11::class_<MultiIntKeyDict<T>>(m, type.c_str())
        // methods
        .def("__contains__",
            &MultiIntKeyDict<T>::in,
            pybind11::is_operator())
        .def("__getitem__",
            [](MultiIntKeyDict<T>& self, long int const& i) {
                return (T) self[i];
            },
            pybind11::is_operator())
        .def("__iter__",
            [](MultiIntKeyDict<T>& self) {
                return pybind11::make_iterator(self.begin(), self.end());
            },
            pybind11::keep_alive<0, 1>())
        .def("__len__",
            &MultiIntKeyDict<T>::size);
}

PYBIND11_MODULE(bind, m) {

    m.def("getLinesHalfEdge", &getLinesHalfEdge);
    m.def("getLinesJunction", &getLinesJunction);

    /*
     *  [tools.hpp]
     *
     */

    declare_MultiIntKeyDict_class<SPVertex>(m, "SPVertex");
    declare_MultiIntKeyDict_class<Cell>(m, "Cell");
    declare_MultiIntKeyDict_class<Face>(m, "Face");
    declare_MultiIntKeyDict_class<Junction>(m, "Junction");

    /*
     *  [mesh.hpp]
     *
     */

    pybind11::class_<Vertex>(m, "Vertex")
        // attributes
        .def_property_readonly("index",
            &Vertex::getIndex)
        .def_property_readonly("position",
            [](Vertex& self) {
                return pybind11::array_t<double>({2}, self.getPosition());
            })
        .def_property_readonly("halfEdgeIndex",
            [](Vertex& self) { return *self.getHalfEdgeIndex(); });

    pybind11::class_<HalfEdge>(m, "HalfEdge")
        // attributes
        .def_property_readonly("index",
            &HalfEdge::getIndex)
        .def_property_readonly("fromIndex",
            [](HalfEdge& self) { return *self.getFromIndex(); })
        .def_property_readonly("toIndex",
            [](HalfEdge& self) { return *self.getToIndex(); })
        .def_property_readonly("previousIndex",
            [](HalfEdge& self) { return *self.getPreviousIndex(); })
        .def_property_readonly("nextIndex",
            [](HalfEdge& self) { return *self.getNextIndex(); })
        .def_property_readonly("pairIndex",
            [](HalfEdge& self) { return *self.getPairIndex(); });

    pybind11::class_<Mesh>(m, "Mesh");

    /*
     *  [system.hpp]
     *
     *  Also defines some wrappers for functions inherited from Mesh in
     *  mesh.hpp.
     *
     */

    pybind11::class_<SPVertex>(m, "SPVertex")
        // attributes
        .def_property_readonly("vertexIndex",
            &SPVertex::getVertexIndex)
        .def_property_readonly("v0",
            &SPVertex::getv0)
        .def_property_readonly("Dr",
            &SPVertex::getDr)
        .def_property_readonly("theta",
            [](SPVertex& self) { return *self.gettheta(); });

    pybind11::class_<Cell>(m, "Cell")
        // attributes
        .def_property_readonly("vertexIndex",
            &Cell::getVertexIndex)
        .def_property_readonly("area",
            [](Cell& cell) { return *cell.getArea(); })
        .def_property_readonly("kA",
            &Cell::getkA)
        .def_property_readonly("A0",
            &Cell::getA0)
        .def_property_readonly("perimeter",
            [](Cell& cell) { return *cell.getPerimeter(); })
        .def_property_readonly("kP",
            &Cell::getkP)
        .def_property_readonly("p0",
            &Cell::getp0)
        .def_property_readonly("P0",
            &Cell::getP0);

    pybind11::class_<Face>(m, "Face")
        // attributes
        .def_property_readonly("halfEdgeIndex",
            &Face::getHalfEdgeIndex);

    pybind11::class_<Junction>(m, "Junction")
        // attributes
        .def_property_readonly("halfEdgeIndex",
            &Junction::getHalfEdgeIndex);

    pybind11::class_<VertexModel, Mesh>(m, "VertexModel")
        // constructor
        .def(pybind11::init
            <long int const&, double const&, double const&, double const&>(),
            pybind11::arg("seed")=0,
            pybind11::arg("v0")=1e-1,
            pybind11::arg("Dr")=1e-1,
            pybind11::arg("p0")=3.81)
        // attributes
        .def_property_readonly("systemSize",
            [](VertexModel& self) {
                return pybind11::array_t<double>({2}, self.getSystemSize());
            })
        .def_property_readonly("vertices",
            &VertexModel::getVertices)
        .def_property_readonly("halfEdges",
            &VertexModel::getHalfEdges)
        .def_property_readonly("sPVertices",
            &VertexModel::getSPVertices)
        .def_property_readonly("cells",
            &VertexModel::getCells)
        .def_property_readonly("faces",
            &VertexModel::getFaces)
        .def_property_readonly("junctions",
            &VertexModel::getJunctions)
        .def_property_readonly("time",
            &VertexModel::getTime)
        // methods
        .def("wrap",
            [](VertexModel& self, pybind11::array_t<double> pos) {
                pybind11::array_t<double>
                    wrapPos({2}, (double*) pos.request().ptr);
                self.wrap((double*) wrapPos.request().ptr);
                return wrapPos;
            })
        .def("wrapDiff",
            [](VertexModel& self,
                pybind11::array_t<double> fromPos,
                pybind11::array_t<double> toPos) {
                std::vector<double> const disp =
                    self.wrapDiff(
                        (double*) fromPos.request().ptr,
                        (double*) toPos.request().ptr);
                return pybind11::array_t<double>({2}, &(disp[0]));
            })
        .def("reset",
            &VertexModel::reset)
        .def("integrate",
            [](VertexModel& self,
                long int const& niter, double const& dt,
                double const& delta, double const& epsilon) {
                for (long int i=0; i < niter; i++) {
                    self.integrate(dt, delta, epsilon);
                }
            },
            pybind11::arg("niter"),
            pybind11::arg("dt")=0,
            pybind11::arg("delta")=0.1,
            pybind11::arg("epsilon")=0.1)
        .def("getForces",
            &VertexModel::getForces)
//         .def("mergeVertices",
//             &VertexModel::mergeVertices,
//             pybind11::arg("halfEdgeIndex"))
//         .def("createJunction",
//             &VertexModel::createJunction,
//             pybind11::arg("halfEdgeIndex0"),
//             pybind11::arg("halfEdgeIndex1"),
//             pybind11::arg("angle"),
//             pybind11::arg("length")=1)
        .def("initRegularTriangularLattice",
            &VertexModel::initRegularTriangularLattice,
            pybind11::arg("size")=6,
            pybind11::arg("junctionLength")=1)
        .def("checkMesh",
            &VertexModel::checkMesh);
}

