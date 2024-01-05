/*
Provide templates for pickling operations on C++ objects via pybind11.
*/

#ifndef PICKLE_HPP
#define PICKLE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <string>

#include "base_forces.hpp"
#include "class_factory.hpp"
#include "system.hpp"

/*
General functions to save to and load from pickle.
*/

void checkSize(pybind11::tuple const& t, pybind11::size_t const& correctSize) {
/*
Check that pybind11::tuple has correct size to generate data.
*/
    if (t.size() != correctSize) {
        throw std::runtime_error("Invalid state (size does not match).");
    }
}

template<class T> pybind11::tuple pybind11_getstate(T const& obj)
/*
python __getstate__ (converts data from object to pickleable tuple)
*/
    { return pybind11::make_tuple(obj); }

template<class T> T pybind11_setstate(pybind11::tuple const& t)
/*
python __setstate__ (converts tuple back to object)
*/
    { checkSize(t, 1); return t[0].cast<T>(); }

/*
Specific functions to save to and load force computation objects from pickle.
*/

// VertexForce<ForcesType>

template<class DerivedType>
pybind11::tuple pybind11_getstate_vertex_forces(                // python __getstate__
    std::shared_ptr<VertexForce<ForcesType>> const& obj) {      // default behaviour
    return pybind11::make_tuple(obj->getType(), obj->getParameters());
}

template<class DerivedType>                                        
void pybind11_setstate_vertex_forces(                           // python __setstate__
    pybind11::tuple const& t, VertexModel& vm) {                // default behaviour
    // check
    checkSize(t, 2);
    // get data
    std::string const type = t[0].cast<std::string>();
    std::map<std::string, double> const parameters = t[1].cast<std::string>();
    // add force
    vm.addVertexForce<DerivedType, ParametersType const&>(
        type, parameters);
}

// HalfEdgeForce<ForcesType>

template<class DerivedType>
pybind11::tuple pybind11_getstate_half_edge_forces(             // python __getstate__
    std::shared_ptr<HalfEdgeForce<ForcesType>> const& obj) {    // default behaviour
    return pybind11::make_tuple(obj->getType(), obj->getParameters());
}

template<class DerivedType>                                        
void pybind11_setstate_half_edge_forces(                        // python __setstate__
    pybind11::tuple const& t, VertexModel& vm) {                // default behaviour
    // check
    checkSize(t, 2);
    // get data
    std::string const type = t[0].cast<std::string>();
    std::map<std::string, double> const parameters = t[1].cast<std::string>();
    // add force
    vm.addHalfEdgeForce<DerivedType, ParametersType const&>(
        type, parameters);
}

#endif

