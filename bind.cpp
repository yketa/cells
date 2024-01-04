/*
Bind C++ objects to python using pybind11.
https://pybind11.readthedocs.io/en/stable/index.html
*/

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
#include <pybind11/stl.h>   // https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html

#include <string>
#include <vector>
#include <map>

#include "forces.hpp"
#include "mesh.hpp"
#include "system.hpp"
#include "tools.hpp"
// #include "pickle.hpp"
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
//         // pickle
//         .def(pybind11::pickle(
//             &pybind11_getstate_MultiIntKeyDict<T>,
//             &pybind11_setstate_MultiIntKeyDict<T>))
        ;
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
        .def_property_readonly("type",
            &Vertex::getType)
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
//         // pickle
//         .def(pybind11::pickle(
//             &pybind11_getstate<Vertex>,
//             &pybind11_setstate<Vertex>))
        ;

    pybind11::class_<HalfEdge>(m, "HalfEdge")
        // attributes
        .def_property_readonly("index",
            &HalfEdge::getIndex)
        .def_property_readonly("type",
            &HalfEdge::getType)
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
//         // pickle
//         .def(pybind11::pickle(
//             &pybind11_getstate<HalfEdge>,
//             &pybind11_setstate<HalfEdge>))
        ;

    pybind11::class_<Mesh>(m, "Mesh")
//         // pickle
//         .def(pybind11::pickle(
//             &pybind11_getstate<Mesh>,
//             &pybind11_setstate<Mesh>))
        ;

    /*
     *  [base_forces.hpp]
     *
     */

    pybind11::class_<BaseForce<ForcesType>,
        std::shared_ptr<BaseForce<ForcesType>>>(
        m, "BaseForce",
        "Python wrapper around C++ force computation object.")
        // attributes
        .def_property_readonly("type",
            &BaseForce<ForcesType>::getType)
        .def_property_readonly("parameters",
            &BaseForce<ForcesType>::getParameters);

    pybind11::class_<HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<HalfEdgeForce<ForcesType>>>(
        m, "HalfEdgeForce",
        "Python wrapper around C++ vertex-based force computation object.");

    pybind11::class_<VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<VertexForce<ForcesType>>>(
        m, "VertexForce",
        "Python wrapper around C++ half-edge-based force computation object.");

    /*
     *  [system.hpp]
     *
     *  Also defines some wrappers for functions inherited from Mesh in
     *  mesh.hpp.
     *
     */

    pybind11::class_<VertexModel, Mesh>(m, "VertexModel",
        "Python wrapper around C++ vertex model simulation object.")
        // constructor
        .def(pybind11::init
            <long int const&>(),
            "Parameters\n"
            "----------\n"
            "seed : int\n"
            "    Random number generator seed. (default: 0)",
            pybind11::arg("seed")=0)
        // attributes [Mesh (mesh.hpp)]
        .def_property_readonly("vertices",
            &VertexModel::getVertices)
        .def_property_readonly("halfEdges",
            &VertexModel::getHalfEdges)
        .def_property_readonly("systemSize",
            [](VertexModel const& self) {
                std::vector<double> const systemSize = self.getSystemSize();
                return pybind11::array_t<double>({2}, &(systemSize[0]));
            })
        // attributes [VertexModel (system.hpp)]
        .def_property_readonly("halfEdgeForces",
            [](VertexModel const& self)
                { return (self.getHalfEdgeForces()).map(); })
        .def_property_readonly("vertexForces",
            [](VertexModel const& self)
                { return (self.getVertexForces()).map(); })
        .def_property_readonly("seed",
            &VertexModel::getSeed)
        .def_property_readonly("time",
            &VertexModel::getTime)
        .def_property_readonly("nT1",
            &VertexModel::getnT1)
        // methods [Mesh (mesh.hpp)]
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
        .def("clear",
            &VertexModel::clear,
            "Clear all data.")
        .def("checkMesh",
            &VertexModel::checkMesh,
            "Check that the vertices and half-edges define a planar mesh,\n"
            "with anticlockwise triangles.")
        // methods [VertexModel (system.hpp)]
        .def("nintegrate",
            [](VertexModel& self,
                long int const& niter, double const& dt,
                double const& delta, double const& epsilon) {
                for (long int i=0; i < niter; i++) {
                    self.integrate(dt, delta, epsilon);
                }
            },
            "Compute forces and integrate dynamics. T1s are performed at\n"
            "each time step.\n"
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
            [](VertexModel& self)
                { self.computeForces(); return self.getForces(); },
            "Compute forces on vertices."
            "\n"
            "Returns\n"
            "-------\n"
            "forces : {int: list} dict\n"
            "    Dictionary which associates vertex indices to force applied\n"
            "    on the vertex.")
        .def("removeHalfEdgeForce",
            &VertexModel::removeHalfEdgeForce,
            "Remove half-edge force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Name of the force to remove.",
            pybind11::arg("name"))
        .def("removeVertexForce",
            &VertexModel::removeVertexForce,
            "Remove vertex force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Name of the force to remove.",
            pybind11::arg("name"))
        // add forces methods [HalfEdgeForce, VertexForce (base_forces.hpp, forces.hpp, forces.cpp)]
        .def("addPerimeterForce",
            &VertexModel::addVertexForce<PerimeterForce,
                double const&, double const&>,
            "Add perimeter force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "kP : float\n"
            "    Perimeter elasticity.\n"
            "P0 : float\n"
            "    Target perimeter.",
            pybind11::arg("name"),
            pybind11::arg("kP"),
            pybind11::arg("P0"))
        .def("addAreaForce",
            &VertexModel::addVertexForce<AreaForce,
                double const&, double const&>,
            "Add area force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "kA : float\n"
            "    Area elasticity.\n"
            "A0 : float\n"
            "    Target area.",
            pybind11::arg("name"),
            pybind11::arg("kA"),
            pybind11::arg("A0"))
        .def("addActiveBrownianForce",
            &VertexModel::addVertexForce<ActiveBrownianForce,
                double const&, double const&>,
            "Add perimeter force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "v0 : float\n"
            "    Self-propulsion velocity.\n"
            "taup : float\n"
            "    Persistence time.",
            pybind11::arg("name"),
            pybind11::arg("v0"),
            pybind11::arg("taup"))
        // initialisation methods [VertexModel (initialisation.cpp)]
        .def("initRegularTriangularLattice",
            &VertexModel::initRegularTriangularLattice,
            "Initialise a regular triangular lattice.\n"
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
        .def("initOpenRegularTriangularLattice",
            &VertexModel::initOpenRegularTriangularLattice,
            "Initialises a regular triangular lattice with a cell replaced\n"
            "with a hole. (TESTING)\n"
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
        .def("initOpenRegularHexagonalLattice",
            &VertexModel::initOpenRegularHexagonalLattice,
            "Initialise a regular square lattice with open outer bondary.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "nCells : int\n"
            "    Number of cells. (default: 1)\n"
            "junctionLength : float\n"
            "    Length of nearest neighbour junctions. (default: 1)",
            pybind11::arg("nCells")=1,
            pybind11::arg("junctionLength")=1)
//         // pickle
//         .def(pybind11::pickle(
//             &pybind11_getstate<VertexModel>,
//             &pybind11_setstate<VertexModel>))
            ;
}

