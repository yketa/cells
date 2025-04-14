/*
Bind C++ objects to python using pybind11.
https://pybind11.readthedocs.io/en/stable/index.html
*/

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
#include <pybind11/stl.h>   // https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "analysis.hpp"
#include "forces.hpp"
#include "forces_pickle.hpp"
#include "integrators.hpp"
#include "integrators_pickle.hpp"
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
        pybind11::vectorize([](double const& x, double const& y) {
            return pmod(angle2(x, y)
                + std::numbers::pi, 2*std::numbers::pi)
                - std::numbers::pi;
        }),
        "Angle of vector (x, y) with respect to horizontal axis in [-pi, pi].",
        pybind11::arg("x"),
        pybind11::arg("y"));

    m.def("hexagonEdgeLength",
        pybind11::vectorize(hexagonEdgeLength),
        "Edge length of a regular hexagon defined by its area.",
        pybind11::arg("hexagonArea"));

    /*
     *  [plot.hpp]
     *
     */

    m.def("getLinesHalfEdge", &getLinesHalfEdge);
    m.def("getLinesJunction", &getLinesJunction);
    m.def("getPolygonsCell", &getPolygonsCell);
    m.def("getMaximumFeretAxesCell", &getMaximumFeretAxesCell);

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

    pybind11::class_<VolumeForce,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<VolumeForce>>(
        m, "VolumeForce",
        "Python wrapper around C++ volume restoring force computation object.")
        .def_property("height",
            &VolumeForce::getHeight,
            &VolumeForce::setHeight)
        .def_property_readonly("heightVelocity",
            &VolumeForce::getHeightVelocity);

    pybind11::class_<LinearVolumeForce,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<LinearVolumeForce>>(
        m, "LinearVolumeForce",
        "Python wrapper around C++ volume restoring force which is linear in\n"
        "edge length computation object.")
        .def_property("height",
            &LinearVolumeForce::getHeight,
            &LinearVolumeForce::setHeight)
        .def_property_readonly("heightVelocity",
            &LinearVolumeForce::getHeightVelocity);

    pybind11::class_<GrowingAreaPerimeterForce,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<GrowingAreaPerimeterForce>>(
        m, "GrowingAreaPerimeterForce",
        "Python wrapper around C++ area and perimeter restoring force with\n"
        "growing target area computation object.\n")
        .def_property("targetArea",
            &GrowingAreaPerimeterForce::getTargetArea,
            &GrowingAreaPerimeterForce::setTargetArea);

    pybind11::class_<ActiveBrownianForce,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<ActiveBrownianForce>>(
        m, "ActiveBrownianForce",
        "Python wrapper around C++ active Brownian force on vertices\n"
        "computation object.")
        .def_property("theta",
            &ActiveBrownianForce::getTheta,
            &ActiveBrownianForce::setTheta);

    pybind11::class_<ActiveBrownianCellForce,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<ActiveBrownianCellForce>>(
        m, "ActiveBrownianCellForce",
        "Python wrapper around C++ active Brownian force on cell centres\n"
        "computation object.")
        .def_property("theta",
            &ActiveBrownianCellForce::getTheta,
            &ActiveBrownianCellForce::setTheta);

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
     *  [forces.hpp]
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
     *  [forces.hpp]
     *  KERATIN
     *
     */

    pybind11::class_<KeratinModel,
        VertexForce<ForcesType>, BaseForce<ForcesType>,
        std::shared_ptr<KeratinModel>>(
        m, "KeratinModel",
        "Python wrapper around C++ keratin model computation object.")
        .def_property("keratin",
            &KeratinModel::getKeratin,
            &KeratinModel::setKeratin)
        .def_property("targetArea",
            &KeratinModel::getTargetArea,
            &KeratinModel::setTargetArea)
        .def_property_readonly("pressure",
            &KeratinModel::getPressure)
        .def_property_readonly("area",
            &KeratinModel::getArea)
        .def_property_readonly("tension",
            &KeratinModel::getTension)
        .def_property_readonly("tension_junction",
            &KeratinModel::getTensionJunction);

    /*
     *  [base_integrators.hpp]
     *
     */

    pybind11::class_<BaseIntegrator<ForcesType, VelocitiesType>,
        std::shared_ptr<BaseIntegrator<ForcesType, VelocitiesType>>>(
        m, "BaseIntegrator",
        "Python wrapper around C++ integrator object.")
        // attributes
        .def_property_readonly("parameters",
            &BaseIntegrator<ForcesType, VelocitiesType>::getParameters)
        // methods
        .def("integrate",
            &BaseIntegrator<ForcesType, VelocitiesType>::integrate,
            "Integrate velocities.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "dt : float\n"
            "    Integration time step.",
            pybind11::arg("dt"));

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
        .def_property_readonly("centre",
            [](VertexModel const& self) {
                std::vector<double> const centre = self.getCentre();
                return pybind11::array_t<double>(2, &(centre[0]));
            })
        // attributes [VertexModel (system.hpp)]
        .def_property_readonly("forces",
            &VertexModel::getForces)
        .def_property_readonly("halfEdgeForces",
            [](VertexModel const& self)
                { return (self.getHalfEdgeForces()).map(); })
        .def_property_readonly("vertexForces",
            [](VertexModel const& self)
                { return (self.getVertexForces()).map(); })
        .def_property_readonly("velocities",
            &VertexModel::getVelocities)
        .def_property_readonly("integrator",
            &VertexModel::getIntegrator)
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
                pybind11::array_t<double> const& toPos,
                bool const& unit=false) {
                double const* fromPosPTR =
                    (double const*) fromPos.request().ptr;
                double const* toPosPTR =
                    (double const*) toPos.request().ptr;
                std::vector<double> const disp =
                    self.wrapDiff(
                        std::vector<double>(fromPosPTR, fromPosPTR + 2),
                        std::vector<double>(toPosPTR, toPosPTR + 2),
                        unit);
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
            "unit : bool\n"
            "    Return unitary vector. (default: False)\n"
            "\n"
            "Returns\n"
            "-------\n"
            "disp : (2,) float numpy array\n"
            "    Difference vector.",
            pybind11::arg("fromPos"),
            pybind11::arg("toPos"),
            pybind11::arg("unit")=false)
        .def("getHalfEdgeIndex",
            &VertexModel::getHalfEdgeIndex,
            "Index of the half-edge between two vertices.\n"
            "Throws error if this half-edge does not exist.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "fromVertexIndex : int\n"
            "    Index of origin vertex.\n"
            "toVertexIndex : int\n"
            "    Index of destination vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "halfEdgeIndex : int\n"
            "    Index of half-edge.",
            pybind11::arg("fromVertexIndex"),
            pybind11::arg("toVertexIndex"))
        .def("getHalfEdgeBetweenIndex",
            &VertexModel::getHalfEdgeBetweenIndex,
            "Index of a half-edge within a triangle to which vertex\n"
            "`vertexIndex0' belongs and which is paired to a half-edge in a\n"
            "triangle to which vertex `vertexIndex1' belongs.\n"
            "Throws error if this half-edge does not exist.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex0 : int\n"
            "    Index of first vertex.\n"
            "vertexIndex1 : int\n"
            "    Index of second vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "halfEdgeIndex : int\n"
            "    Index of half-edge.",
            pybind11::arg("vertexIndex0"),
            pybind11::arg("vertexIndex1"))
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
        .def("moveToNeighboursBarycentre",
            &VertexModel::moveToNeighboursBarycentre,
            "Move vertex to centre of mass (= barycentre) of its neighbours.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex : int\n"
            "    Index of vertex.",
            pybind11::arg("vertexIndex"))
        .def("moveToNeighboursCentroid",
            &VertexModel::moveToNeighboursCentroid,
            "Move vertex to centroid of (= arithmetic mean of all the points\n"
            "on) the surface enclosed by its neighbours.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex : int\n"
            "    Index of vertex.",
            pybind11::arg("vertexIndex"))
        .def("deleteEdge",
            &VertexModel::deleteEdge,
            "Delete edge by merging two vertices.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIntex : int\n"
            "    Index of half-edge linking two vertices to be merged.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "deletedVertexIndex : int\n"
            "    Index of deleted vertex.\n"
            "deletedHalfEdgeIndices : tuple of int\n"
            "    Indices of deleted half-edges.",
            pybind11::arg("halfEdgeIntex"))
        .def("createEdge",
            &VertexModel::createEdge,
            "Create edge by splitting one vertex.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIndex0 : int\n"
            "    Index of first half-edge going out of a vertex, and from\n"
            "    whose pair half-edge it will be separated after the\n"
            "    introduction of a new junction.\n"
            "halfEdgeIndex1 : int\n"
            "    Index of second half-edge going out of the same vertex, and\n"
            "    from whose pair half-edge it will be separated after the\n"
            "    introduction of a new junction.\n"
            "angle : float\n"
            "    Angle of the new junction with respect to the horizontal\n"
            "    axis.\n"
            "length : float\n"
            "    Length to set for the new junction. (default: 1)\n"
            "type0 : str\n"
            "    Name of the type of half-edge for the half-edge going from\n"
            "    the splitted vertex to the created vertex. (default: \"\")\n"
            "type1 : str\n"
            "    Name of the type of half-edge for the half-edge going from\n"
            "    the created vertex to the splitted vertex. (default: \"\")\n"
            "NOTE: Other created half-edges will inherit the type of their\n"
            "      pair half-edge, and the created vertex inheretit the type\n"
            "      of the splitted vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "createdVertexIndex : int\n"
            "    Index of created vertex.\n"
            "createdHalfEdgeIndices : tuple of int\n"
            "    Indices of created half-edges.",
            pybind11::arg("halfEdgeIndex0"),
            pybind11::arg("halfEdgeIndex1"),
            pybind11::arg("angle"),
            pybind11::arg("length")=1,
            pybind11::arg("type0")="",
            pybind11::arg("type1")="")
        .def("mergeVertices",
            &VertexModel::mergeVertices,
            "Merge two vertices sharing an edge.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIndex : int\n"
            "    Index of half-edge from the first to second vertex where\n"
            "    the first will be merged into the second.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "deletedVertexIndex : int\n"
            "    Index of vertex which was merged into the other.\n"
            "    NOTE: This is fromVertexIndex.\n"
            "deletedHalfEdgeIndices : tuple of int\n"
            "    Indices of half-edges deleted by this operation.",
            pybind11::arg("halfEdgeIndex"))
        .def("splitVertices",
            &VertexModel::splitVertices,
            "Split a vertex into two vertices along an edge.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIntex : int\n"
            "    Index of half-edge in the middle of which a vertex will be\n"
            "    created.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "createdVertexIndex : int\n"
            "    Index of created vertex.\n"
            "createdHalfEdgeIndices : tuple of int\n"
            "    Indices of created half-edges.",
            pybind11::arg("halfEdgeIntex"))
        .def("swapEdge",
            &VertexModel::swapEdge,
            "Remove an edge and replace it with an other edge between the\n"
            "two triangle corners which the original edge did not link.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIndex : int\n"
            "    Index of half-edge belonging to the edge which will be\n"
            "    removed.\n"
            "type0 : str\n"
            "    Type to give to the first of the created half-edges.\n"
            "    (default: \"\")\n"
            "type1 : str\n"
            "    Type to give to the second of the created half-edges.\n"
            "    (default: \"\")\n"
            "\n"
            "Returns\n"
            "-------\n"
            "newHalfEdgeIndex0 : int\n"
            "    Index of the first of the created half-edges.\n"
            "newHalfEdgeIndex1 : int\n"
            "    Index of the second of the created half-edges.",
            pybind11::arg("halfEdgeIndex"),
            pybind11::arg("type0")="",
            pybind11::arg("type1")="")
        .def("changeToBoundary",
            &VertexModel::changeToBoundary,
            "Change vertex to a boundary vertex with identical attributes\n"
            "(including index) and empty type.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex :\n"
            "    Vertex to change to a boundary vertex.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "deletedVertexIndex :\n"
            "    Index of vertex which was changed to boundary vertex.\n"
            "    NOTE: This is `vertexIndex'.",
            pybind11::arg("vertexIndex"))
        .def("scale",
            &VertexModel::scale,
            "Rescale all lengths.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "scalingFactor : float\n"
            "    Factor by which to multiply all lengths.",
            pybind11::arg("scalingFactor"))
        .def("setSystemSize",
            [](VertexModel& self, pybind11::array_t<double> const& L) {
                auto l = L.unchecked<>();
                self.setSystemSize({l[0], l[1]});
            },
            "Set system size and move vertices according to displacement of\n"
            "system centre.\n"
            "\n"
            "WARNING: This should only be performed on systems which are NOT\n"
            "         periodic.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "systemSize : (2,) float array-like\n"
            "    New system size.",
            pybind11::arg("L"))
        .def("getVertexIndicesByType",
            &VertexModel::getVertexIndicesByType,
            "Return vertex indices corresponding to type.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "type : str\n"
            "    Type of vertices.\n"
            "    NOTE: if type == \"\" then return all vertices.\n"
            "exclude_boundary : bool\n"
            "    Exclude boundary vertices. (default: True)\n"
            "\n"
            "Returns\n"
            "-------\n"
            "vertexIndices : list of int\n"
            "    Vertex indices corresponding to type.",
            pybind11::arg("type"),
            pybind11::arg("exclude_boundary")=true)
        .def("getHalfEdgeIndicesByType",
            &VertexModel::getHalfEdgeIndicesByType,
            "Return half-edge indices corresponding to type.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "type : str\n"
            "    Type of half-edges.\n"
            "    NOTE: if type == \"\" then return all half-edges.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "halfEdgeIndices : list of int\n"
            "    Half-edge indices corresponding to type.",
            pybind11::arg("type"))
        .def("getCentreHalfEdges",
            [](VertexModel const& self) {
                std::vector<long int> centreHalfEdgeIndices(0);
                std::map<long int, HalfEdge> const& halfEdges =
                    self.getHalfEdges();
                std::vector<long int> const vertexIndices =
                    self.getVertexIndicesByType("centre");
                for (long int vertexIndex : vertexIndices) {
                    std::vector<long int> halfEdgeToNeighbourIndices =
                        self.getNeighbourVertices(vertexIndex)[1];
                    for (long int halfEdgeIndex : halfEdgeToNeighbourIndices) {
                        centreHalfEdgeIndices.push_back(
                            halfEdgeIndex);
                        centreHalfEdgeIndices.push_back(
                            (halfEdges.at(halfEdgeIndex)).getPairIndex());
                    }
                }
                return centreHalfEdgeIndices;
            },
            "Return all half-edges connected to a cell centre.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "halfEdgeIndices : list of int\n"
            "    Half-edge indices.\n")
        // methods [VertexModel (system.hpp)]
        .def("setSeed",
            &VertexModel::setSeed,
            "Change seed and re-initialise random number generator.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "seed : int\n"
            "    Random number generator seed.",
            pybind11::arg("seed"))
        .def("getNeighbouringCellIndices",
            &VertexModel::getNeighbouringCellIndices,
            "Compute neighbouring cell indices.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "vertexIndex : int\n"
            "    Index of vertex from which to compute neighbouring cell\n"
            "    indices.\n"
            "    NOTE: This vertex should be surrounded by cells from which\n"
            "          it is separated by single edges.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "neighbourCellIndices : list of int\n"
            "    Indices of neighbouring cells.",
            pybind11::arg("vertexIndex"))
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
        .def("splitCell",
            &VertexModel::splitCell,
            "Split one cell into two by adding vertices in the middle of\n"
            "edges with half-edges with index `halfEdgeIndex0' and\n"
            "`halfEdgeIndex1' and adding a junction between these vertices.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIndex0 : int\n"
            "    First half-edge on the middle of which to add a vertex.\n"
            "halfEdgeIndex1 : int\n"
            "    Second half-edge on the middle of which to add a vertex.\n"
            "NOTE: `halfEdgeIndex0' and `halfEdgeIndex1' must belong to the\n"
            "      boundary of the same cell.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "newCellVertexIndex : int\n"
            "    Index of newly created cell centre.",
            pybind11::arg("halfEdgeIndex0"),
            pybind11::arg("halfEdgeIndex1"))
        .def("splitCellAtMax",
            &VertexModel::splitCellAtMax,
            "Split cell across the two most separated boundary half-edge\n"
            "centres.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "cellVertexIndex : int\n"
            "    Vertex index of centre of cell to split.\n"
            "avoidThreeEdgeCells : bool\n"
            "    Ignore divisions which would result in cells with three\n"
            "    edges. (default: False)\n"
            "\n"
            "Returns\n"
            "-------\n"
            "newCellVertexIndex : int\n"
            "    Index of newly created cell centre.",
            pybind11::arg("cellVertexIndex"),
            pybind11::arg("avoidThreeEdgeCells")=false)
        .def("mergeCell",
            &VertexModel::mergeCell,
            "Merge two cells whose centres are separated by the half-edge\n"
            "with index `halfEdgeIndex'.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeIndex : int\n"
            "    Half-edge separating the two cell centres to merge.\n"
            "    NOTE: Two cell centres should not be connected by any\n"
            "          half-edge. This half-edge should belong to the\n"
            "          junction between the cells.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "toCellVertexIndex : int\n"
            "    Index of cell centre to which the origin cell centre was\n"
            "    merged."
            "neighbouringCellIndices : list of int\n"
            "    Indices of neighbouring cell centres whose edges have been\n"
            "     removed.",
            pybind11::arg("halfEdgeIndex"))
        .def("mergeCellAtMin",
            &VertexModel::mergeCellAtMin,
            "Merge cell with its neighbouring cell of minimum area.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "cellVertexIndex : int\n"
            "    Index of centre of cell to merge with its neighbour.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "toCellVertexIndex : int\n"
            "    Index of cell centre to which the origin cell centre was\n"
            "    merged."
            "neighbouringCellIndices : list of int\n"
            "    Indices of neighbouring cell centres whose edges have been\n"
            "     removed.",
            pybind11::arg("mergeCellAtMin"))
        .def("checkMesh",
            &VertexModel::checkMesh,
            "Check that the vertices and half-edges define a planar mesh.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "halfEdgeTypes : list of str\n"
            "    These types should not be defined identically for both\n"
            "    half-edges of a single pair (default: []).\n"
            "checkOrientations : bool\n"
            "    Check that triangles have anticlockwise orientation.\n"
            "    (default: True)",
            pybind11::arg("helfEdgeTypes")=std::vector<std::string>(),
            pybind11::arg("checkOrientations")=true)
        .def("getPositions",
            [](VertexModel const& self, bool const wrapped=true) {
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
        // remove forces methods [VertexModel (system.hpp)]
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
        .def("addSurfaceForce",
            &VertexModel::addVertexForce<SurfaceForce,
                double const&, double const&>,
            "Add cell surface restoring force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "Lambda : float\n"
            "    Surface tension.\n"
            "V0 : float\n"
            "    Cell volume.",
            pybind11::arg("name"),
            pybind11::arg("Lambda"),
            pybind11::arg("V0"))
        .def("addVolumeForce",
            &VertexModel::addVertexForce<VolumeForce,
                double const&, double const&, double const&>,
            "Add cell volume restoring force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "kV : float\n"
            "    Volume elasticity.\n"
            "H0 : float\n"
            "    Target height.\n"
            "A0 : float\n"
            "    Target area.",
            pybind11::arg("name"),
            pybind11::arg("kV"),
            pybind11::arg("H0"),
            pybind11::arg("A0"))
        .def("addLinearVolumeForce",
            &VertexModel::addVertexForce<LinearVolumeForce,
                double const&, double const&,
                double const&, double const&,
                double const&, double const&, double const&>,
            "Add cell volume restoring force which is linear in edge length.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "kA : float\n"
            "    Area elasticity.\n"
            "A0 : float\n"
            "    Target area."
            "kP : float\n"
            "    Perimeter elasticity.\n"
            "P0 : float\n"
            "    Target perimeter."
            "taur : float\n"
            "    Volume relaxation time.\n"
            "H0 : float\n"
            "    Target area-to-volume ratio."
            "taua : float\n"
            "    Height differences relaxation time.\n",
            pybind11::arg("name"),
            pybind11::arg("kA"),
            pybind11::arg("A0"),
            pybind11::arg("kP"),
            pybind11::arg("P0"),
            pybind11::arg("taur"),
            pybind11::arg("H0"),
            pybind11::arg("taua"))
        .def("addPressureForce",
            &VertexModel::addVertexForce<PressureForce,
                double const&, bool const&>,
            "Add pressure force.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "F : float\n"
            "    Force scale.\n"
            "fixedForce : bool\n"
            "    Force norm is constant and equal to F for all vertices\n"
            "    rather than derives from a pressure term. (default: True)\n"
            "    NOTE: if not(fixedForce) then F sets the product of the\n"
            "          pressure, the number of boundary vertices and the\n"
            "          boundary perimeter.",
            pybind11::arg("name"),
            pybind11::arg("F"),
            pybind11::arg("fixedForce")=true)
        .def("addGrowingAreaPerimeterForce",
            &VertexModel::addVertexForce<GrowingAreaPerimeterForce,
                double const&, double const&,
                double const&, double const&>,
            "Add area and perimeter force with growing target area.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "kA : float\n"
            "    Area elasticity.\n"
            "s0 : float\n"
            "    Shape inde.\n"
            "A0 : float\n"
            "    Reference area.\n"
            "tauA : float\n"
            "    Time scale of increase of target area.\n"
            "    NOTE: `A0'/`tauA' is the target area increase rate.",
            pybind11::arg("name"),
            pybind11::arg("kA"),
            pybind11::arg("s0"),
            pybind11::arg("A0"),
            pybind11::arg("tauA"))
        .def("addBoundaryTension",
            &VertexModel::addVertexForce<BoundaryTension,
                double const&>,
            "Add line tension on open boundaries.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "gamma : float\n"
            "    Line tension.",
            pybind11::arg("name"),
            pybind11::arg("gamma"))
        .def("addActiveBrownianForce",
            &VertexModel::addVertexForce<ActiveBrownianCellForce,
                double const&, double const&>,
            "Add active Brownian force on cell centres.\n"
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
        // Model0-4
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
        // KeratinModel
        .def("addKeratinModel",
            &VertexModel::addVertexForce<KeratinModel,
                double const&, double const&, double const&,
                double const&, double const&,
                double const&, double const&, double const&,
                double const&, double const&, double const&>,
            "Add keratin model.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "name : str\n"
            "    Unique name for the force.\n"
            "K : float\n"
            "    Area elasticity.\n"
            "A0 : float\n"
            "    Minimum target area.\n"
            "taur : float\n"
            "    Target area relaxation time scale.\n"
            "Gamma : float\n"
            "    Perimeter elasticity.\n"
            "p0 : float\n"
            "    Target shape index.\n"
            "alpha : float\n"
            "    Keratin pressure scale.\n"
            "beta : float\n"
            "    Keratin to area elasticity and relaxation parameter.\n"
            "kth : float\n"
            "    Keratin concentration threshold.\n"
            "tau : float\n"
            "    Keratin concentration evolution time scale.\n"
            "sigma : float\n"
            "    Keratin concentration noise standard deviation.\n"
            "ron : float\n"
            "    Keratin concentration on-rate evolution time rate\n"
            "    (=1/tauon).\n",
            pybind11::arg("name"),
            pybind11::arg("K"),
            pybind11::arg("A0"),
            pybind11::arg("taur"),
            pybind11::arg("Gamma"),
            pybind11::arg("p0"),
            pybind11::arg("alpha"),
            pybind11::arg("beta"),
            pybind11::arg("kth"),
            pybind11::arg("tau"),
            pybind11::arg("sigma"),
            pybind11::arg("ron"))
        // set integrator [BaseIntegrator (integrators.hpp, integrators.hpp, integrators.cpp)]
        .def("setUnitOverdampedIntegrator",
            &VertexModel::setIntegrator<UnitOverdamped>,
            "Unit vertex-substrate drag coefficient overdamped integrator.")
        .def("setOverdampedIntegrator",
            &VertexModel::setIntegrator<Overdamped, double const&>,
            "Custom vertex-substrate drag coefficient overdamped integrator.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "xi : float\n"
            "    Vertex drag coefficient.",
            pybind11::arg("xi"))
        .def("setPairFrictionIntegrator",
            &VertexModel::setIntegrator<PairFriction, double const&>,
            "Cell corner pair friction and unit vertex-substrate drag\n"
            "coefficient overdamped integrator.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "eta : float\n"
            "    Ratio of vertex-vertex over vertex-substrate friction\n"
            "    friction coefficient.",
            pybind11::arg("eta"))
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
            "hexagonArea : float\n"
            "    Area of regular hexagons. (default: 1)",
            pybind11::arg("size")=6,
            pybind11::arg("hexagonArea")=1)
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
            "hexagonArea : float\n"
            "    Area of regular hexagons. (default: 1)",
            pybind11::arg("size")=6,
            pybind11::arg("hexagonArea")=1)
        .def("initOpenRegularHexagonalLattice",
            &VertexModel::initOpenRegularHexagonalLattice,
            "Initialise a regular square lattice with open outer boundary.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "nCells : int\n"
            "    Number of cells. (default: 1)\n"
            "    NOTE: This must be a squared integer.\n"
            "hexagonArea : float\n"
            "    Area of regular hexagons. (default: 1)\n"
            "boxLength : float\n"
            "    Length of the box in units of sqrt(`nCells')*junctionLength\n"
            "    where junctionLength is the length of a side of a regular\n"
            "    hexagon. (default: 3)",
            pybind11::arg("nCells")=1,
            pybind11::arg("hexagonArea")=1,
            pybind11::arg("boxLength")=3)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<VertexModel>,
            &pybind11_setstate<std::unique_ptr<VertexModel>>))
        // analysis
        .def("getCentreVelocities", &getCentreVelocities,
            "Compute velocities of the centroid of each cell.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "centreVelocities : {int: list} dict\n"
            "    Dictionary which associates cell centre vertex indices to\n"
            "    the velocity of the centroid of the cell.")
        .def("getVectorsToNeighbouringCells", &getVectorsToNeighbouringCells,
            "Compute vectors to neighbouring cell centres.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "vectorsToNeighbours : {int: list} dict\n"
            "    Dictionary which associates cell centre vertex indices to\n"
            "    list of vectors to neighbouring cells.")
        .def("getPAticOrderParameters", &getPAticOrderParameters,
            "Compute p-atic order parameters of cell centres.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "p : int\n"
            "    Symmetry parameter. (default: 6)\n"
            "\n"
            "Returns\n"
            "-------\n"
            "psip : {int: complex} dict\n"
            "    Dictionary which associates cell centre vertex indices to\n"
            "    p-atic order parameter.",
            pybind11::arg("p")=6);

    /*
     *  Miscellaneous
     *
     */

    m.def("getMaxLengthCells", &getMaxLengthCells,
        "Return maximum length between two cell corners in each cell,\n"
        "considering periodic boundary conditions.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "vm : cells.bind.VertexModel\n"
        "    Vertex model object.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "maxLength : {int: float}\n"
        "    Dictionary which associates cell centre vertex indices to the\n"
        "    maximum length between two cell corners in this cell.",
        pybind11::arg("vm"));

    m.def("getMaxLengthBoundaries", &getMaxLengthBoundaries,
        "Return maximum length between two cell corners around each boundary\n"
        "vertex, ignoring periodic boundary conditions.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "vm : cells.bind.VertexModel\n"
        "    Vertex model object.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "maxLength : {int: float}\n"
        "    Dictionary which associates boundary vertex indices to the\n"
        "    maximum length between two cell corners around this boundary\n"
        "    vertex.",
        pybind11::arg("vm"));

    m.def("getPercentageKeptNeighbours", &getPercentageKeptNeighbours,
        "Return percentage of kept cell neighbours between two states of a\n"
        "vertex model.\n"
        "\n"
        "NOTE: It is assumed that cell centre indices match.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "vm0 : cells.bind.VertexModel\n"
        "    Initial vertex model object.\n"
        "vm1 : cells.bind.VertexModel\n"
        "    Final vertex model object.\n"
        "a : float\n"
        "    Maximum distance between cell centres in unwrapped coordinates\n"
        "    for them to be considered neighbours, in addition to sharing a\n"
        "    junction. (default: inf)\n"
        "\n"
        "Returns\n"
        "-------\n"
        "pct : {int: float} dict\n"
        "    Dictionary which associates cell centre vertex indices to the\n"
        "    percentage of kept cell neighbours.",
        pybind11::arg("vm0"),
        pybind11::arg("vm1"),
        pybind11::arg("a")=std::numeric_limits<double>::infinity());

    m.def("getMaximumFeretAnglesCell", &getMaximumFeretAnglesCell,
        "Return maximum Feret axes and angles of Feret axes with respect to\n"
        "the centre of the boundary.\n"
        "\n"
        "NOTE: There should be exactly 1 boundary vertex.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "vm : VertexModel\n"
        "    Vertex model object from which to compute Feret axes and\n"
        "    angles.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "feretAxes : (*, 2) float list\n"
        "    Vector of unit vectors corresponding to maximum Feret diameter\n"
        "    of each cell.\n"
        "feretAngles : (*, 2) float list\n"
        "    Vector of tuple of norm of the radius between the cell centre\n"
        "    and the centre of boundary vertices, and angle between the\n"
        "    maximum Feret diameter and this radius.\n",
        pybind11::arg("vm"));

    m.def("getAllWaveVectors2D", &getAllWaveVectors2D,
        "Return wave vectors associated to rectangular box.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "L : float or (1,)- or (2,) float array-like\n"
        "    Size of the box.\n"
        "qmin : float\n"
        "    Minimum wave vector norm.\n"
        "qmax : float\n"
        "    Maximum wave vector norm.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "wv : (*, 2) float numpy array\n"
        "    Array of (2\\pi/L nx, 2\\pi/L ny) wave vectors corresponding to\n"
        "    to the target interval [`qmin', `qmax'].\n"
        "    NOTE: Only a single vector of each pair of opposite wave\n"
        "          vectors is returned. Here it is chosen such that ny >= 0.",
        pybind11::arg("L"),
        pybind11::arg("qmin"),
        pybind11::arg("qmax"));

    m.def("getAllFT2D", &getAllFT2D,
        "Return 2D Fourier transform of delta-peaked values.\n"
        "\n"
        ".. math::"
        "V(k_l) = \\sum_i \\exp(-1i k_l \\cdot r_i) v_i\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "positions : (*, 2) float array-like\n"
        "    Positions r_i of delta-peaked values.\n"
        "L : float or (1,)- or (2,) float array-like\n"
        "    Size of the box.\n"
        "values : (*,) complex array-like\n"
        "    Delta-peaked values v_i.\n"
        "qmin : float\n"
        "    Minimum wave vector norm.\n"
        "qmax : float\n"
        "    Maximum wave vector norm.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "ft : (**,) complex numpy array\n"
        "    Fourier transform of `values' for each wave vector in the\n"
        "    target  norm interval [`qmin', `qmax'].\n"
        "    NOTE: These are given by getWaveVectors2D.\n",
        pybind11::arg("positions"),
        pybind11::arg("L"),
        pybind11::arg("values"),
        pybind11::arg("qmin"),
        pybind11::arg("qmax"));

}

