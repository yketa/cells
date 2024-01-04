/*
Provide pickling ability for C++ VertexModel and associated objects via
pybind11.
*/

#ifndef SYSTEM_PICKLE_HPP
#define SYSTEM_PICKLE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>
#include <map>

#include "pickle.hpp"
#include "mesh.hpp"
#include "system.hpp"
#include "tools.hpp"

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
        v.getPosition(),
        v.getUPosition(),
        v.getHalfEdgeIndex(),
        v.getBoundary(),
        v.getType());
}

template<>
Vertex pybind11_setstate<Vertex>(pybind11::tuple const& t) {
    checkSize(t, 6);
    long int const index =
        t[0].cast<long int>();
    std::vector<double> position =
        t[1].cast<std::vector<double>>();
    std::vector<double> uposition =
        t[2].cast<std::vector<double>>();
    long int const halfEdgeIndex =
        t[3].cast<long int>();
    bool const boundary =
        t[4].cast<bool>();
    std::string const type =
        t[5].cast<std::string>();
    return Vertex(index, position, uposition, halfEdgeIndex, boundary, type);
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
        he.getPairIndex(),
        he.getType());
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
    std::string const type = t[6].cast<std::string>();
    return HalfEdge(
        index, fromIndex, toIndex, previousIndex, nextIndex, pairIndex, type);
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

