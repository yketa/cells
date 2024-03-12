/*
Provide pickling ability for C++ VertexModel and associated objects via
pybind11.
*/

#ifndef SYSTEM_PICKLE_HPP
#define SYSTEM_PICKLE_HPP

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <sstream>
#include <vector>

#include "base_pickle.hpp"
#include "mesh.hpp"
#include "random.hpp"
#include "system.hpp"
#include "tools.hpp"

/*
 *  [random.hpp]
 *
 */

// Random

template<>
pybind11::tuple pybind11_getstate<Random>(Random const& rnd) {
    // extract to string stream
    std::stringstream ss;
    ss << rnd;
    // save as bytes
    // https://github.com/htm-community/htm.core/issues/160
    return pybind11::make_tuple(pybind11::bytes(ss.str()));
}

template<>
Random pybind11_setstate<Random>(pybind11::tuple const& t) {
    // check
    checkSize(t, 1);
    // convert data to string stream
    pybind11::bytes const b = t[0].cast<pybind11::bytes>();
    std::istringstream ss(b);
    // create and set random generator from string stream
    Random random;
    ss >> random;
    // return object
    return random;
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
    checkSize(t, 7);
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
        pybind11_getstate_map<long int, Vertex>(m.getVertices()),
        pybind11_getstate_map<long int, HalfEdge>(m.getHalfEdges()),
        m.getSystemSize());
}

template<>
Mesh pybind11_setstate<Mesh>(pybind11::tuple const& t) {
    checkSize(t, 3);
    std::map<long int, Vertex> const vertices =
        pybind11_setstate_map<long int, Vertex>(
            t[0].cast<pybind11::tuple>());
    std::map<long int, HalfEdge> const halfEdges =
        pybind11_setstate_map<long int, HalfEdge>(
            t[1].cast<pybind11::tuple>());
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
        vm.getSeed(),
        pybind11_getstate<Random>(vm.getRandom()),
        vm.getTime(),
        vm.getnT1(),
        pybind11_getstate_force_class_factory<VertexForce<ForcesType>>(
            vm.getVertexForces()),
        pybind11_getstate_force_class_factory<HalfEdgeForce<ForcesType>>(
            vm.getHalfEdgeForces()));
}

template<>
std::unique_ptr<VertexModel> pybind11_setstate<std::unique_ptr<VertexModel>>(   // using unique_ptr to conserve pointers in force computation objects
    pybind11::tuple const& t) {
    // check
    checkSize(t, 7);
    // extract data
    Mesh const mesh =
        pybind11_setstate<Mesh>(t[0].cast<pybind11::tuple>());
    long int const seed =
        t[1].cast<long int>();
    Random const random =
        pybind11_setstate<Random>(t[2].cast<pybind11::tuple>());
    double const time =
        t[3].cast<double>();
    long int const nT1 =
        t[4].cast<long int>();
    // initialise simulation object
    std::unique_ptr<VertexModel> vm =
        std::make_unique<VertexModel>(mesh, seed, time, nT1);
    // set forces
    pybind11_setstate_force_class_factory<VertexForce<ForcesType>>(
        *vm, t[5].cast<std::map<std::string, pybind11::tuple>>());
    pybind11_setstate_force_class_factory<HalfEdgeForce<ForcesType>>(
        *vm, t[6].cast<std::map<std::string, pybind11::tuple>>());
    // copy random generator
    // (to be done after forces initialisation which may call the generator)
    vm->setRandom(random);
    return vm;
}

#endif

