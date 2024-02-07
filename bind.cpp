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
#include "forces_pickle.hpp"
#include "mesh.hpp"
#include "plot.hpp"
#include "system_pickle.hpp"
#include "system.hpp"
#include "tools.hpp"

PYBIND11_MODULE(bind, m) {
    m.doc() =
        "Module bind wraps C++ objects and functions to integrate and plot\n"
        "vertex model.\n";

    /*
     *  [tools.hpp]
     *
     */

    m.def("angle2",
        [](double const& x, double const& y) { return angle2(x, y); },
        "Angle of 2D vector (x, y) with respect to horizontax axis.",
        pybind11::arg("x"),
        pybind11::arg("y"));

    /*
     *  [plot.hpp]
     *
     */

    m.def("getLinesHalfEdge", &getLinesHalfEdge);
    m.def("getLinesJunction", &getLinesJunction);
    m.def("getPolygonsCell", &getPolygonsCell);

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
                return pybind11::array_t<double>(2, &(position[0]));
            })
        .def_property_readonly("uposition",
            [](Vertex const& self) {
                std::vector<double> const uposition = self.getUPosition();
                return pybind11::array_t<double>(2, &(uposition[0]));
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
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<HalfEdge>,
            &pybind11_setstate<HalfEdge>))
        ;

    pybind11::class_<Mesh>(m, "Mesh")
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Mesh>,
            &pybind11_setstate<Mesh>))
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
            &BaseForce<ForcesType>::getParameters)
        // methods
        .def("integrate",
            &BaseForce<ForcesType>::integrate,
            "Integrate internal degrees of freedom.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "dt : float\n"
            "    Integration time step.",
            pybind11::arg("dt"));

    pybind11::class_<HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<HalfEdgeForce<ForcesType>>>(
        m, "HalfEdgeForce",
        "Python wrapper around C++ vertex-based force computation object.");

    pybind11::class_<VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<VertexForce<ForcesType>>>(
        m, "VertexForce",
        "Python wrapper around C++ half-edge-based force computation object.");

    /*
     *  [forces.hpp]
     *
     */

    pybind11::class_<ActiveBrownianForce,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<ActiveBrownianForce>>(
        m, "ActiveBrownianForce",
        "Python wrapper around C++ active Brownian force computation object.")
        .def_property("theta",
            &ActiveBrownianForce::getTheta,
            &ActiveBrownianForce::setTheta);

    pybind11::class_<OrnsteinUhlenbeckTension,
        HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<OrnsteinUhlenbeckTension>>(
        m, "OrnsteinUhlenbeckTension",
        "Python wrapper around C++ Ornstein-Uhlenbeck tension computation\n"
        "object.")
        .def_property("tension",
            &OrnsteinUhlenbeckTension::getTension,
            &OrnsteinUhlenbeckTension::setTension);

    /*
     *  MODELS 0-4
     *
     */

    pybind11::class_<Model0,
        HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<Model0>>(
        m, "Model0",
        "Python wrapper around C++ model 0 computation object.")
        .def_property("perimeter",
            &Model0::getPerimeter,
            &Model0::setPerimeter)
        .def_property("noise",
            &Model0::getNoise,
            &Model0::setNoise)
        .def_property("tension",
            &Model0::getTension,
            &Model0::setTension);

    pybind11::class_<Model1,
        HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<Model1>>(
        m, "Model1",
        "Python wrapper around C++ model 1 computation object.")
        .def_property("perimeter",
            &Model1::getPerimeter,
            &Model1::setPerimeter)
        .def_property("tension",
            &Model1::getTension,
            &Model1::setTension);

    pybind11::class_<Model2,
        HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<Model2>>(
        m, "Model2",
        "Python wrapper around C++ model 2 computation object.")
        .def_property("length",
            &Model2::getLength,
            &Model2::setLength)
        .def_property("restLength",
            &Model2::getRestLength,
            &Model2::setRestLength)
        .def_property("noise",
            &Model2::getNoise,
            &Model2::setNoise)
        .def_property("tension",
            &Model2::getTension,
            &Model2::setTension);

    pybind11::class_<Model3,
        HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<Model3>>(
        m, "Model3",
        "Python wrapper around C++ model 3 computation object.")
        .def_property("length",
            &Model3::getLength,
            &Model3::setLength)
        .def_property("restLength",
            &Model3::getRestLength,
            &Model3::setRestLength)
        .def_property("tension",
            &Model3::getTension,
            &Model3::setTension);

    pybind11::class_<Model4,
        HalfEdgeForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<Model4>>(
        m, "Model4",
        "Python wrapper around C++ model 4 computation object.")
        .def_property("length",
            &Model4::getLength,
            &Model4::setLength)
        .def_property("restLength",
            &Model4::getRestLength,
            &Model4::setRestLength)
        .def_property("tension",
            &Model4::getTension,
            &Model4::setTension);

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
        .def(pybind11::init
            <VertexModel const&>(),
            "Copy constructor for VertexModel.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vm : VertexModel\n"
            "    Instance of VertexModel to copy.",
            pybind11::arg("vm"))
        // attributes [Mesh (mesh.hpp)]
        .def_property_readonly("vertices",
            &VertexModel::getVertices)
        .def_property_readonly("halfEdges",
            &VertexModel::getHalfEdges)
        .def_property_readonly("systemSize",
            [](VertexModel const& self) {
                std::vector<double> const systemSize = self.getSystemSize();
                return pybind11::array_t<double>(2, &(systemSize[0]));
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
                return pybind11::array_t<double>(2, &(wpos[0]));
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
                return pybind11::array_t<double>(2, &(disp[0]));
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
        .def("getNeighbourVertices",
            [](VertexModel const& self, long int const& vertexIndex) {
                std::vector<std::vector<long int>> neighbours =
                    self.getNeighbourVertices(vertexIndex);
                long int const nNeigh = neighbours[0].size();
                return pybind11::make_tuple(
                    pybind11::array_t<long int>(nNeigh, &(neighbours[0][0])),
                    pybind11::array_t<long int>(nNeigh, &(neighbours[1][0])));
            },
            "Indices of neighbouring vertices and indices of half-edges\n"
            "towards them.\n"
            "\n"
            "NOTE: If the half-edge construction is correct, these should be\n"
            "      in anti-clockwise order.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex : int\n"
            "    Index of the vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "neighbourVerticesIndices : (*,) int numpy array\n"
            "    Indices of vertices neighbouring this vertex.\n"
            "halfEdgesToNeighboursIndices : (*,) int numpy array\n"
            "    Indices of half-edges from this vertex towards neighbour\n"
            "    vertices.",
            pybind11::arg("vertexIndex"))
        .def("getVertexToNeighboursArea",
            &VertexModel::getVertexToNeighboursArea,
            "Area encapsulated by the neighbours of a vertex (shoelace\n"
            "formula).\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex : int\n"
            "    Index of vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "area : double\n"
            "    Area encapsulated by neighbours.",
            pybind11::arg("vertexIndex"))
        .def("getVertexToNeighboursPerimeter",
            &VertexModel::getVertexToNeighboursPerimeter,
            "Perimeter encapsulated by the neighbours of a vertex.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex : int\n"
            "    Index of vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "perimeter : double\n"
            "    Perimeter encapsulated by neighbours.",
            pybind11::arg("vertexIndex"))
        .def("clear",
            &VertexModel::clear,
            "Clear all data.")
        .def("checkMesh",
            &VertexModel::checkMesh,
            "Check that the vertices and half-edges define a planar mesh,\n"
            "with anticlockwise triangles.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeTypes : list of str\n"
            "    These types should not be defined identically for both\n"
            "    half-edges of a single pair.",
            pybind11::arg("helfEdgeTypes")=std::vector<std::string>())
        .def("getVertexIndicesByType",
            [](VertexModel& self, std::string const type) {
                if (type == "") { throw std::runtime_error("No type."); }
                std::map<long int, Vertex> const vertices =
                    self.getVertices();
                std::vector<long int> vertexIndices;
                for (auto it=vertices.begin(); it != vertices.end(); ++it) {
                    if ((it->second).getType() == type) {
                        vertexIndices.push_back(it->first);
                    }
                }
                return vertexIndices;
            },
            "Return vertex indices corresponding to type.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "type : str\n"
            "    Type of vertices.\n"
            "    NOTE: This cannot be empty.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "vertexIndices : list of int\n"
            "    Vertex indices corresponding to type.",
            pybind11::arg("type"))
        .def("getHalfEdgeIndicesByType",
            [](VertexModel& self, std::string const type) {
                if (type == "") { throw std::runtime_error("No type."); }
                std::map<long int, HalfEdge> const halfEdges =
                    self.getHalfEdges();
                std::vector<long int> halfEdgeIndices;
                for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
                    if ((it->second).getType() == type) {
                        halfEdgeIndices.push_back(it->first);
                    }
                }
                return halfEdgeIndices;
            },
            "Return half-edge indices corresponding to type.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "type : str\n"
            "    Type of half-edges.\n"
            "    NOTE: This cannot be empty.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "halfEdgeIndices : list of int\n"
            "    Half-edge indices corresponding to type.",
            pybind11::arg("type"))
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
        .def("getPositions",
            [](VertexModel& self, bool const wrapped=true) {
                std::map<long int, std::vector<double>> positions;
                std::map<long int, Vertex> const vertices = self.getVertices();
                for (auto it=vertices.begin(); it != vertices.end(); ++it) {
                    if (!(it->second).getBoundary()) {
                        if (wrapped) {
                            positions.emplace(
                                it->first, (it->second).getPosition());
                        }
                        else {
                            positions.emplace(
                                it->first, (it->second).getUPosition());
                        }
                    }
                }
                return positions;
            },
            "Return positions of all vertices excluding boundary vertices.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "wrapped : bool\n"
            "    Wrap the position around periodic boundary conditions.\n"
            "    (default: True)\n"
            "\n"
            "Returns\n"
            "-------\n"
            "positions : {int: list} dict\n"
            "    Dictionary which associates vertex indices to position of\n"
            "    vertex.",
            pybind11::arg("wrapped")=true)
        .def("getForces",
            [](VertexModel& self)
                { self.computeForces(); return self.getForces(); },
            "Compute forces on vertices.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "forces : {int: list} dict\n"
            "    Dictionary which associates vertex indices to force applied\n"
            "    on the vertex.")
        .def("getForcesCentre",
            [](VertexModel& self) {
                std::map<long int, std::vector<double>> const forces =
                    self.getForces();
                std::map<long int, std::vector<double>> forcesCentre;
                std::map<long int, Vertex> const vertices = self.getVertices();
                for (auto it=vertices.begin(); it != vertices.end(); ++it) {
                    if ((it->second).getType() == "centre") {       // loop over cell centre
                        forcesCentre[it->first] = {0, 0};
                        std::vector<long int> const neighbours =
                            self.getNeighbourVertices(it->first)[0];
                        long int const nNeighbours = neighbours.size();
                        for (long int i : neighbours) {             // loop over neighbours
                            for (int dim=0; dim < 2; dim++) {
                                if (inMap(forces, i)) {
                                    forcesCentre[it->first][dim] += // average force on neighbours
                                        (forces.at(i)).at(dim)/nNeighbours;
                                }
                            }
                        }
                    }
                }
                return forcesCentre;
            },
            "Compute forces on cell centres.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "forces : {int: list} dict\n"
            "    Dictionary which associates cell centre vertex indices to\n"
            "    force applied on the vertex.")
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
            "Add cell perimeter restoring force.\n"
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
            "Add cell area restoring force.\n"
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
            "Add active Brownian force on vertices.\n"
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
        .def("addOrnsteinUhlenbeckTension",
            &VertexModel::addHalfEdgeForce<OrnsteinUhlenbeckTension,
                double const&, double const&, double const&>,
            "Add Ornstein-Uhlenbeck tension on junctions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "t0 : float\n"
            "    Standard deviation of tension.\n"
            "taup : float\n"
            "    Persistence time.",
            pybind11::arg("name"),
            pybind11::arg("t0"),
            pybind11::arg("st0"),
            pybind11::arg("taup"))
        /*
         *  MODELS 0-4
         *
         */
        .def("addModel0",
            &VertexModel::addHalfEdgeForce<Model0,
                double const&, double const&, double const&, double const&>,
            "Add model 0 for active junctions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "Gamma : float\n"
            "    Perimeter elasticity constant.\n"
            "P0 : float\n"
            "    Target parameter.\n"
            "sigma : float\n"
            "    Tension noise strength.\n"
            "taup : float\n"
            "    Noise persistence time.",
            pybind11::arg("name"),
            pybind11::arg("Gamma"),
            pybind11::arg("P0"),
            pybind11::arg("sigma"),
            pybind11::arg("taup"))
        .def("addModel1",
            &VertexModel::addHalfEdgeForce<Model1,
                double const&, double const&, double const&, double const&>,
            "Add model 1 for active junctions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "Gamma : float\n"
            "    Perimeter elasticity constant.\n"
            "P0 : float\n"
            "    Target parameter.\n"
            "sigma : float\n"
            "    Tension noise strength.\n"
            "taup : float\n"
            "    Noise persistence time.",
            pybind11::arg("name"),
            pybind11::arg("Gamma"),
            pybind11::arg("P0"),
            pybind11::arg("sigma"),
            pybind11::arg("taup"))
        .def("addModel2",
            &VertexModel::addHalfEdgeForce<Model2,
                double const&, double const&, double const&, double const&>,
            "Add model 2 for active junctions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "Gamma : float\n"
            "    Junction elasticity constant.\n"
            "taur : float\n"
            "    Rest length relaxation time.\n"
            "sigma : float\n"
            "    Tension noise strength.\n"
            "taup : float\n"
            "    Noise persistence time.",
            pybind11::arg("name"),
            pybind11::arg("Gamma"),
            pybind11::arg("taur"),
            pybind11::arg("sigma"),
            pybind11::arg("taup"))
        .def("addModel3",
            &VertexModel::addHalfEdgeForce<Model3,
                double const&, double const&, double const&>,
            "Add model 3 for active junctions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "Gamma : float\n"
            "    Junction elasticity constant.\n"
            "sigma : float\n"
            "    Rest length noise strength.\n"
            "taup : float\n"
            "    Noise persistence time.",
            pybind11::arg("name"),
            pybind11::arg("Gamma"),
            pybind11::arg("sigma"),
            pybind11::arg("taup"))
        .def("addModel4",
            &VertexModel::addHalfEdgeForce<Model4,
                double const&, double const&, double const&, double const&>,
            "Add model 2 for active junctions.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "Gamma : float\n"
            "    Junction elasticity constant.\n"
            "taur : float\n"
            "    Rest length relaxation time.\n"
            "sigma : float\n"
            "    Tension noise strength.\n"
            "taup : float\n"
            "    Noise persistence time.",
            pybind11::arg("name"),
            pybind11::arg("Gamma"),
            pybind11::arg("taur"),
            pybind11::arg("sigma"),
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
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<VertexModel>,
            &pybind11_setstate<std::unique_ptr<VertexModel>>));

}

