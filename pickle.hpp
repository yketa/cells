/*
Provide pickling ability for C++ objects via pybind11.
*/

#ifndef PICKLE_HPP
#define PICKLE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>   // https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html

#include <vector>
#include <map>

#include "mesh.hpp"
#include "system.hpp"
#include "tools.hpp"

/*
Functions to save to and load from pickle.
*/

template<class T> pybind11::tuple pybind11_getstate(T const& obj);  // python __getstate__ (converts data from object to pickleable tuple)
template<class T> T pybind11_setstate(pybind11::tuple const& t);    // python __setstate__ (converts tuple back to object)

void checkSize(pybind11::tuple const& t, pybind11::size_t const& correctSize) {
/*
Check that pybind11::tuple has correct size to generate data.
*/
    if (t.size() != correctSize) {
        throw std::runtime_error("Invalid state (size does not match).");
    }
}

/*
 *  long int
 *
 */

template<>
pybind11::tuple pybind11_getstate<long int>(long int const& i)
    { return pybind11::make_tuple(i); }
template<>
long int pybind11_setstate<long int>(pybind11::tuple const& t)
    { checkSize(t, 1); return t[0].cast<long int>(); }

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
        v.getBoundary(),
        v.getPosition(),
        v.getUPosition(),
        v.getHalfEdgeIndex());
}
template<>
Vertex pybind11_setstate<Vertex>(pybind11::tuple const& t) {
    checkSize(t, 5);
    long int const index =
        t[0].cast<long int>();
    bool const boundary =
        t[1].cast<bool>();
    std::vector<double> position =
        t[2].cast<std::vector<double>>();
    std::vector<double> uposition =
        t[3].cast<std::vector<double>>();
    long int const halfEdgeIndex =
        t[4].cast<long int>();
    return Vertex(index, boundary, position, uposition, halfEdgeIndex);
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

#endif

