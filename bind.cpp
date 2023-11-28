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
#include "plot.hpp"

/*
Functions to save to and load from pickle.
*/

template<class T> pybind11::tuple pybind11_getstate(T const& obj);  // python __getstate__
template<class T> T pybind11_setstate(pybind11::tuple const& t);    // python __setstate__

void checkSize(pybind11::tuple const& t, pybind11::size_t const& correctSize) {
/*
Check that pybind11::tuple has correct size to generate data.
*/
    if (t.size() != correctSize) {
        throw std::runtime_error("Invalid state (size does not match).");
    }
}

/*
 *  std::map<long int, TT>
 *
 */

template<class TT> std::map<long int, pybind11::tuple>
    pybind11_getstate_map(std::map<long int, TT> const& map) {
    // extract data from TT
    std::map<long int, pybind11::tuple> data;
    for (auto it=map.begin(); it != map.end(); ++it) {
        data.emplace(it->first, pybind11_getstate<TT>(it->second));
    }
    // return data
    return data;
}
template<class TT> std::map<long int, TT>
    pybind11_setstate_map(std::map<long int, pybind11::tuple> const& data) {
    // convert data from TT
    std::map<long int, TT> map;
    for (auto it=data.begin(); it != data.end(); ++it) {
        map.emplace(it->first, pybind11_setstate<TT>(it->second));
    }
    // return map
    return map;
}

/*
 *  [tools.hpp]
 *
 */

// MultiIntKeyDict<TT>
template<class TT> pybind11::tuple pybind11_getstate_MultiIntKeyDict(
    MultiIntKeyDict<TT> const& mikd) {
    // copy data from MultiIntKeyDict
    std::map<long int, long int> const mikd_keys = mikd.getKeys();
    std::map<long int, TT> const mikd_data = mikd.getData();
    long int const mikd_maxIndex = mikd.getMaxIndex();
    // extract data from TT
    std::map<long int, pybind11::tuple> data =
        pybind11_getstate_map<TT>(mikd_data);
    // return ensemble of data
    return pybind11::make_tuple(
        mikd_keys, data, mikd_maxIndex);
}
template<class TT> MultiIntKeyDict<TT> pybind11_setstate_MultiIntKeyDict(
    pybind11::tuple const& t) {
    // check
    checkSize(t, 3);
    // get data
    std::map<long int, long int> const mikd_keys =
        t[0].cast<std::map<long int, long int>>();
    std::map<long int, pybind11::tuple> const data =
        t[1].cast<std::map<long int, pybind11::tuple>>();
    std::map<long int, TT> const mikd_data =
        pybind11_setstate_map<TT>(data);
    long int const mikd_maxIndex =
        t[2].cast<long int>();
    // create MultiIntKeyDict with data
    MultiIntKeyDict<TT> mikd(mikd_keys, mikd_data, mikd_maxIndex);
    return mikd;
}

/*
 *  [mesh.hpp]
 *
 */

// Vertex
template<>
pybind11::tuple pybind11_getstate<Vertex>(Vertex const& v) {
    return pybind11::make_tuple(
        v.getIndex(),
        v.getPosition(),
        v.getHalfEdgeIndex());
}
template<>
Vertex pybind11_setstate<Vertex>(pybind11::tuple const& t) {
    checkSize(t, 3);
    long int const index =
        t[0].cast<long int>();
    std::vector<double>position =
        t[1].cast<std::vector<double>>();
    long int const halfEdgeIndex =
        t[2].cast<long int>();
    return Vertex(index, position, halfEdgeIndex);
}

// HalfEdge
template<>
pybind11::tuple pybind11_getstate<HalfEdge>(HalfEdge const& he) {
    return pybind11::make_tuple(
        he.getIndex(),
        he.getFromIndex(),
        he.getToIndex(),
        he.getPreviousIndex(),
        he.getNextIndex(),
        he.getPairIndex());
}
template<>
HalfEdge pybind11_setstate<HalfEdge>(pybind11::tuple const& t) {
    checkSize(t, 6);
    long int const index = t[0].cast<long int>();
    long int const fromIndex = t[1].cast<long int>();
    long int const toIndex = t[2].cast<long int>();
    long int const previousIndex = t[3].cast<long int>();
    long int const nextIndex = t[4].cast<long int>();
    long int const pairIndex = t[5].cast<long int>();
    return HalfEdge(
        index, fromIndex, toIndex, previousIndex, nextIndex, pairIndex);
}

// Mesh
template<>
pybind11::tuple pybind11_getstate<Mesh>(Mesh const& m) {
    return pybind11::make_tuple(
        pybind11_getstate_map<Vertex>(m.getVertices()),
        pybind11_getstate_map<HalfEdge>(m.getHalfEdges()),
        m.getSystemSize());
}
template<>
Mesh pybind11_setstate<Mesh>(pybind11::tuple const& t) {
    checkSize(t, 3);
    std::map<long int, Vertex> const vertices =
        pybind11_setstate_map<Vertex>
            (t[0].cast<std::map<long int, pybind11::tuple>>());
    std::map<long int, HalfEdge> const halfEdges =
        pybind11_setstate_map<HalfEdge>
            (t[1].cast<std::map<long int, pybind11::tuple>>());
    std::vector<double> const systemSize =
        t[2].cast<std::vector<double>>();
    return Mesh(vertices, halfEdges, systemSize);
}

/*
 *  [system.hpp]
 *
 */

// SPVertex
template<>
pybind11::tuple pybind11_getstate<SPVertex>(SPVertex const& spv) {
    return pybind11::make_tuple(
        spv.getVertexIndex(),
        spv.gettheta(),
        spv.getv0(),
        spv.getDr());
}
template<>
SPVertex pybind11_setstate<SPVertex>(pybind11::tuple const& t) {
    checkSize(t, 4);
    long int const vertexIndex = t[0].cast<long int>();
    double const theta = t[1].cast<double>();
    double const v0 = t[2].cast<double>();
    double const Dr = t[3].cast<double>();
    return SPVertex(vertexIndex, theta, v0, Dr);
}

// Cell
template<>
pybind11::tuple pybind11_getstate<Cell>(Cell const& cell) {
    return pybind11::make_tuple(
        cell.getVertexIndex(),
        cell.getArea(),
        cell.getkA(),
        cell.getA0(),
        cell.getPerimeter(),
        cell.getkP(),
        cell.getp0());
}
template<>
Cell pybind11_setstate<Cell>(pybind11::tuple const& t) {
    checkSize(t, 7);
    long int const vertexIndex = t[0].cast<long int>();
    double const area = t[1].cast<double>();
    double const kA = t[2].cast<double>();
    double const A0 = t[3].cast<double>();
    double const perimeter = t[4].cast<double>();
    double const kP = t[5].cast<double>();
    double const p0 = t[6].cast<double>();
    Cell cell(vertexIndex, A0, p0, kA, kP);
    cell.setArea(area);
    cell.setPerimeter(perimeter);
    return cell;
}

// Face
template<>
pybind11::tuple pybind11_getstate<Face>(Face const& face) {
    return pybind11::make_tuple(
        face.getHalfEdgeIndex());
}
template<>
Face pybind11_setstate<Face>(pybind11::tuple const& t) {
    checkSize(t, 1);
    long int const halfEdgeIndex = t[0].cast<long int>();
    return Face(halfEdgeIndex);
}

// Junction
template<>
pybind11::tuple pybind11_getstate<Junction>(Junction const& j) {
    return pybind11::make_tuple(
        j.getHalfEdgeIndex());
}
template<>
Junction pybind11_setstate<Junction>(pybind11::tuple const& t) {
    checkSize(t, 1);
    long int const halfEdgeIndex = t[0].cast<long int>();
    return Junction(halfEdgeIndex);
}

// VertexModel
template<>
pybind11::tuple pybind11_getstate<VertexModel>(VertexModel const& vm) {
    return pybind11::make_tuple(
        pybind11_getstate<Mesh>(vm),
        pybind11_getstate_MultiIntKeyDict<SPVertex>(vm.getSPVertices()),
        pybind11_getstate_MultiIntKeyDict<Cell>(vm.getCells()),
        pybind11_getstate_MultiIntKeyDict<Face>(vm.getFaces()),
        pybind11_getstate_MultiIntKeyDict<Junction>(vm.getJunctions()),
        vm.getv0(),
        vm.getDr(),
        vm.getp0(),
        vm.getSeed(),
        vm.getTime(),
        vm.getnT1());
}
template<>
VertexModel pybind11_setstate<VertexModel>(pybind11::tuple const& t) {
    checkSize(t, 11);
    Mesh const mesh =
        pybind11_setstate<Mesh>
            (t[0].cast<pybind11::tuple>());
    MultiIntKeyDict<SPVertex> const sPVertices =
        pybind11_setstate_MultiIntKeyDict<SPVertex>
            (t[1].cast<pybind11::tuple>());
    MultiIntKeyDict<Cell> const cells =
        pybind11_setstate_MultiIntKeyDict<Cell>
            (t[2].cast<pybind11::tuple>());
    MultiIntKeyDict<Face> const faces =
        pybind11_setstate_MultiIntKeyDict<Face>
            (t[3].cast<pybind11::tuple>());
    MultiIntKeyDict<Junction> const junctions =
        pybind11_setstate_MultiIntKeyDict<Junction>
            (t[4].cast<pybind11::tuple>());
    double const v0 = t[5].cast<double>();
    double const Dr = t[6].cast<double>();
    double const p0 = t[7].cast<double>();
    long int const seed = t[8].cast<long int>();
    double const time = t[9].cast<double>();
    long int const nT1 = t[10].cast<long int>();
    return VertexModel(mesh, sPVertices, cells, faces, junctions,
        v0, Dr, p0, seed, time, nT1);
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
            &MultiIntKeyDict<T>::size)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate_MultiIntKeyDict<T>,
            &pybind11_setstate_MultiIntKeyDict<T>));
}

PYBIND11_MODULE(bind, m) {

    m.def("getLinesHalfEdge", &getLinesHalfEdge);
    m.def("getLinesJunction", &getLinesJunction);
    m.def("getPolygonsCell", &getPolygonsCell);

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
            [](Vertex const& self) {
                std::vector<double> const position = self.getPosition();
                return pybind11::array_t<double>({2}, &(position[0]));
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

    pybind11::class_<SPVertex>(m, "SPVertex")
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

    pybind11::class_<Cell>(m, "Cell")
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

    pybind11::class_<Face>(m, "Face")
        // attributes
        .def_property_readonly("halfEdgeIndex",
            &Face::getHalfEdgeIndex)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Face>,
            &pybind11_setstate<Face>));

    pybind11::class_<Junction>(m, "Junction")
        // attributes
        .def_property_readonly("halfEdgeIndex",
            &Junction::getHalfEdgeIndex)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<Junction>,
            &pybind11_setstate<Junction>));

    pybind11::class_<VertexModel, Mesh>(m, "VertexModel")
        // constructor
        .def(pybind11::init
            <long int const&, double const&, double const&, double const&>(),
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
            })
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
            &VertexModel::checkMesh)
        // pickle
        .def(pybind11::pickle(
            &pybind11_getstate<VertexModel>,
            &pybind11_setstate<VertexModel>));
}

