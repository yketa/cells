/*
Bind C++ objects to python using pybind11.
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "mesh.hpp"
#include "system.hpp"
#include "tools.hpp"

PYBIND11_MODULE(bind, m) {

    /*
     * [tools.hpp]
     *
     */

    pybind11::class_<MultiIntKeyDict<SPVertex>>(m, "MultiIntKeyDict<SPVertex>")
        // methods
        .def("__contains__",
            &MultiIntKeyDict<SPVertex>::in,
            pybind11::is_operator())
        .def("__getitem__",
            [](MultiIntKeyDict<SPVertex>& self, long int const& i) {
                return (SPVertex) self[i];
            },
            pybind11::is_operator());

    pybind11::class_<MultiIntKeyDict<Cell>>(m, "MultiIntKeyDict<Cell>")
        // methods
        .def("__contains__",
            &MultiIntKeyDict<Cell>::in,
            pybind11::is_operator())
        .def("__getitem__",
            [](MultiIntKeyDict<Cell>& self, long int const& i) {
                return (Cell) self[i];
            },
            pybind11::is_operator());

    pybind11::class_<MultiIntKeyDict<Face>>(m, "MultiIntKeyDict<Face>")
        // methods
        .def("__contains__",
            &MultiIntKeyDict<Face>::in,
            pybind11::is_operator())
        .def("__getitem__",
            [](MultiIntKeyDict<Face>& self, long int const& i) {
                return (Face) self[i];
            },
            pybind11::is_operator());

    pybind11::class_<MultiIntKeyDict<Junction>>(m, "MultiIntKeyDict<Junction>")
        // methods
        .def("__contains__",
            &MultiIntKeyDict<Junction>::in,
            pybind11::is_operator())
        .def("__getitem__",
            [](MultiIntKeyDict<Junction>& self, long int const& i) {
                return (Junction) self[i];
            },
            pybind11::is_operator());

    /*
     * [mesh.hpp]
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
     * [system.hpp]
     *
     * Also defines some wrappers for functions inherited from Mesh in
     * mesh.hpp.
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
        .def(pybind11::init<long int const&>(),
            pybind11::arg("seed")=0)
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
        .def("initRegularTriangularLattice",
            &VertexModel::initRegularTriangularLattice,
            pybind11::arg("size")=6,
            pybind11::arg("junctionLength")=1)
        .def("checkMesh",
            &VertexModel::checkMesh);
}

