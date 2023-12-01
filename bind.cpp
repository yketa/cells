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
#include "pickle.hpp"
#include "plot.hpp"

/*
declare class template
https://stackoverflow.com/questions/66227840
https://stackoverflow.com/questions/47487888
iterating through object
https://stackoverflow.com/questions/49257947
*/

template<class T> void declare_MultiIntKeyDict_class(
    pybind11::module& m, std::string const& nameT, std::string const& typeT) {

    std::string const type =
        "MultiIntKeyDict" + nameT;
    std::string const doc =
        "Python wrapper around C++ MultiIntKeyDict<" + typeT + ">.";
    pybind11::class_<MultiIntKeyDict<T>>(m, type.c_str(), doc.c_str())
        //.doc("Python wrapper around C++ MultiIntKeyDict<" + nameT + ">.")
//         // constructor
//         .def(pybind11::init<>())
        // attributes
        .def_property_readonly("keys",
            &MultiIntKeyDict<T>::getKeys)
        .def_property_readonly("data",
            &MultiIntKeyDict<T>::getData)
        // methods
        .def("__contains__",
            &MultiIntKeyDict<T>::in,
            "Is key in dictionary?",
            pybind11::is_operator())
        .def("__getitem__",
            [](MultiIntKeyDict<T>& self, long int const& i) {
                return (T) self[i];
            },
            "Get value associated to key.",
            pybind11::is_operator())
//         .def("__setitem__",
//             [](MultiIntKeyDict<T>& self, long int const& i, T const& v) {
//                 self[i] = v;
//             },
//             pybind11::is_operator())
//         .def("__delitem__",
//             [](MultiIntKeyDict<T>& self, long int const& i) {
//                 self.erase(i);
//             },
//             pybind11::is_operator())
        .def("__iter__",
            [](MultiIntKeyDict<T> const& self) {
                return pybind11::make_iterator(self.begin(), self.end());
            },
            "Iterator over keys.",
            pybind11::keep_alive<0, 1>())
        .def("__len__",
            &MultiIntKeyDict<T>::size,
            "Number of values in the dictionary.")
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate_MultiIntKeyDict<T>,
            &pybind11_setstate_MultiIntKeyDict<T>));
}

PYBIND11_MODULE(bind, m) {
    m.doc() =
        "Module bind wraps C++ objects and functions to integrate and plot\n"
        "vertex model.\n";

    /*
     *  [plot.hpp]
     *
     */

    m.def("getLinesHalfEdge", &getLinesHalfEdge);
    m.def("getLinesJunction", &getLinesJunction);
    m.def("getPolygonsCell", &getPolygonsCell);

    /*
     *  [tools.hpp]
     *
     */

    declare_MultiIntKeyDict_class<long int>(m, "", "long int");
    declare_MultiIntKeyDict_class<SPVertex>(m, "SPVertex", "SPVertex");
    declare_MultiIntKeyDict_class<Cell>(m, "Cell", "Cell");
    declare_MultiIntKeyDict_class<Face>(m, "Face", "Face");
    declare_MultiIntKeyDict_class<Junction>(m, "Junction", "Junction");

    /*
     *  [mesh.hpp]
     *
     */

    pybind11::class_<Vertex>(m, "Vertex")
        // attributes
        .def_property_readonly("index",
            &Vertex::getIndex)
        .def_property_readonly("boundary",
            &Vertex::getBoundary)
        .def_property_readonly("position",
            [](Vertex const& self) {
                std::vector<double> const position = self.getPosition();
                return pybind11::array_t<double>({2}, &(position[0]));
            })
        .def_property_readonly("uposition",
            [](Vertex const& self) {
                std::vector<double> const uposition = self.getUPosition();
                return pybind11::array_t<double>({2}, &(uposition[0]));
            })
        .def_property_readonly("halfEdgeIndex",
            &Vertex::getHalfEdgeIndex)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Vertex>,
            &pybind11_setstate<Vertex>));

    pybind11::class_<HalfEdge>(m, "HalfEdge")
        // attributes
        .def_property_readonly("index",
            &HalfEdge::getIndex)
        .def_property_readonly("fromIndex",
            &HalfEdge::getFromIndex)
        .def_property_readonly("toIndex",
            &HalfEdge::getToIndex)
        .def_property_readonly("previousIndex",
            &HalfEdge::getPreviousIndex)
        .def_property_readonly("nextIndex",
            &HalfEdge::getNextIndex)
        .def_property_readonly("pairIndex",
            &HalfEdge::getPairIndex)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<HalfEdge>,
            &pybind11_setstate<HalfEdge>));

    pybind11::class_<Mesh>(m, "Mesh")
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Mesh>,
            &pybind11_setstate<Mesh>));

    /*
     *  [system.hpp]
     *
     *  Also defines some wrappers for functions inherited from Mesh in
     *  mesh.hpp.
     *
     */

    pybind11::class_<SPVertex>(m, "SPVertex", "Self-propelled vertex.")
        // attributes
        .def_property_readonly("vertexIndex",
            &SPVertex::getVertexIndex)
        .def_property_readonly("v0",
            &SPVertex::getv0)
        .def_property_readonly("Dr",
            &SPVertex::getDr)
        .def_property_readonly("theta",
            &SPVertex::gettheta)
        .def(pybind11::pickle(
            &pybind11_getstate<SPVertex>,
            &pybind11_setstate<SPVertex>))
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<SPVertex>,
            &pybind11_setstate<SPVertex>));

    pybind11::class_<Cell>(m, "Cell", "Physical cell.")
        // attributes
        .def_property_readonly("vertexIndex",
            &Cell::getVertexIndex)
        .def_property_readonly("area",
            &Cell::getArea)
        .def_property_readonly("kA",
            &Cell::getkA)
        .def_property_readonly("A0",
            &Cell::getA0)
        .def_property_readonly("perimeter",
            &Cell::getPerimeter)
        .def_property_readonly("kP",
            &Cell::getkP)
        .def_property_readonly("p0",
            &Cell::getp0)
        .def_property_readonly("P0",
            &Cell::getP0)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Cell>,
            &pybind11_setstate<Cell>));

    pybind11::class_<Face>(m, "Face", "Physical face.")
        // attributes
        .def_property_readonly("halfEdgeIndex",
            &Face::getHalfEdgeIndex)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Face>,
            &pybind11_setstate<Face>));

    pybind11::class_<Junction>(m, "Junction", "Physical junction.")
        // attributes
        .def_property_readonly("halfEdgeIndex",
            &Junction::getHalfEdgeIndex)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Junction>,
            &pybind11_setstate<Junction>));

    pybind11::class_<VertexModel, Mesh>(m, "VertexModel",
        "Python wrapper around C++ vertex model simulation object.")
        // constructor
        .def(pybind11::init
            <long int const&, double const&, double const&, double const&>(),
            "Parameters\n"
            "----------\n"
            "seed : int\n"
            "    Random number generator seed. (default: 0)\n"
            "v0 : float\n"
            "    Vertex self-propulsion velocity. (default: 1e-1)\n"
            "Dr : float\n"
            "    Vertex propulsion rotational diffusion constant.\n"
            "    (default: 1e-1)\n"
            "p0 : float\n"
            "    Dimensionless target perimeter of cell. (default: 3.81)",
            pybind11::arg("seed")=0,
            pybind11::arg("v0")=1e-1,
            pybind11::arg("Dr")=1e-1,
            pybind11::arg("p0")=3.81)
        // attributes
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
        .def_property_readonly("v0",
            &VertexModel::getv0)
        .def_property_readonly("Dr",
            &VertexModel::getDr)
        .def_property_readonly("p0",
            &VertexModel::getp0)
        .def_property_readonly("systemSize",
            [](VertexModel const& self) {
                std::vector<double> const systemSize = self.getSystemSize();
                return pybind11::array_t<double>({2}, &(systemSize[0]));
            })
        .def_property_readonly("seed",
            &VertexModel::getSeed)
        .def_property_readonly("time",
            &VertexModel::getTime)
        .def_property_readonly("nT1",
            &VertexModel::getnT1)
        // methods
        .def("wrap",
            [](VertexModel const& self, pybind11::array_t<double> const& pos) {
                double const* posPTR =
                    (double const*) pos.request().ptr;
                std::vector<double> const wpos =
                    self.wrap(std::vector<double>(posPTR, posPTR + 2));
                return pybind11::array_t<double>({2}, &(wpos[0]));
            },
            "Wrap position to positive values with respect to periodic\n"
            "boundary conditions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "position : (2,) float array-like\n"
            "    Position to wrap.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "wposition : (2,) float numpy array\n"
            "    Wrapped position.",
            pybind11::arg("positions"))
        .def("wrapDiff",
            [](VertexModel const& self,
                pybind11::array_t<double> const& fromPos,
                pybind11::array_t<double> const& toPos) {
                double const* fromPosPTR =
                    (double const*) fromPos.request().ptr;
                double const* toPosPTR =
                    (double const*) toPos.request().ptr;
                std::vector<double> const disp =
                    self.wrapDiff(
                        std::vector<double>(fromPosPTR, fromPosPTR + 2),
                        std::vector<double>(toPosPTR, toPosPTR + 2));
                return pybind11::array_t<double>({2}, &(disp[0]));
            },
            "Wrap difference vector with respect to periodic boundary\n"
            "conditions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "fromPos : (2,) float array-like\n"
            "    Pointer to initial point position.\n"
            "toPos : (2,) float array-like\n"
            "    Pointer to final point position.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "disp : (2,) float numpy array\n"
            "    Difference vector.",
            pybind11::arg("fromPos"),
            pybind11::arg("toPos"))
//         .def("reset",
//             &VertexModel::reset)
        .def("nintegrate",
            [](VertexModel& self,
                long int const& niter, double const& dt,
                double const& delta, double const& epsilon) {
                for (long int i=0; i < niter; i++) {
                    self.integrate(dt, delta, epsilon);
                }
            },
            "Integrate self-propelled vertex model and check for and perform\n"
            "T1s (at each step).\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "niter : int\n"
            "    Number of integration steps to perform.\n"
            "dt : float\n"
            "    Integration time step. (default: 0)\n"
            "delta : float\n"
            "    Distance between vertices below which these should be\n"
            "    merged. (default: 0.1)\n"
            "epsilon : float\n"
            "    Create two vertices at distance `delta' + `epsilon' after\n"
            "    T1. (default: 0.1)",
            pybind11::arg("niter"),
            pybind11::arg("dt")=0,
            pybind11::arg("delta")=0.1,
            pybind11::arg("epsilon")=0.1)
        .def("getForces",
            &VertexModel::getForces,
            "Compute forces on vertices from vertex model."
            "\n"
            "Returns\n"
            "-------\n"
            "forces : {int: list} dict\n"
            "    Dictionary which associates vertex indices to force applied\n"
            "     on the vertex.")
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
            "Initialises a regular triangular lattice.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "size : int\n"
            "    Number of vertices in both horizontal and vertical\n"
            "    directions. (default: 6)\n"
            "    NOTE: This must be a multiple of 6.\n"
            "junctionLength : float\n"
            "    Length of nearest neighbour junctions. (default: 1)",
            pybind11::arg("size")=6,
            pybind11::arg("junctionLength")=1)
        .def("checkMesh",
            &VertexModel::checkMesh,
            "Check that the vertices and half-edges define a planar mesh,\n"
            "with anticlockwise triangles.")
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<VertexModel>,
            &pybind11_setstate<VertexModel>));
}

