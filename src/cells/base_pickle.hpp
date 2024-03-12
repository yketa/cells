/*
Provide templates for pickling operations on C++ objects via pybind11.
*/

#ifndef PICKLE_HPP
#define PICKLE_HPP

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <string>

#include "base_forces.hpp"
#include "class_factory.hpp"
#include "system.hpp"

/*
 *  General functions to save to and load from pickle.
 *
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
Python __getstate__ (converts data from object to pickleable tuple).
*/
    { return pybind11::make_tuple(obj); }

template<class T> T pybind11_setstate(pybind11::tuple const& t)
/*
Python __setstate__ (converts tuple back to object).
*/
    { checkSize(t, 1); return t[0].cast<T>(); }

/*
 *  Specific functions to save and load maps from pickle.
 *
 */

template<class T, class TT>
pybind11::tuple pybind11_getstate_map(std::map<T, TT> const& obj) {
    // convert data
    std::map<pybind11::tuple, pybind11::tuple> data;
    for (auto it=obj.begin(); it != obj.end(); ++it) {
        data.emplace(
            pybind11_getstate<T>(it->first),
            pybind11_getstate<TT>(it->second));
    }
    // return data
    return pybind11::make_tuple(data);
}

template<class T, class TT>
std::map<T, TT> pybind11_setstate_map(pybind11::tuple const& t) {
    // check
    checkSize(t, 1);
    // convert data
    std::map<pybind11::tuple, pybind11::tuple> const data =
        t[0].cast<std::map<pybind11::tuple, pybind11::tuple>>();
    std::map<T, TT> obj;
    for (auto it=data.begin(); it != data.end(); ++it) {
        obj.emplace(
            pybind11_setstate<T>(it->first),
            pybind11_setstate<TT>(it->second));
    }
    // return object
    return obj;
}

/*
 *  Specific functions to save and load force computation object from pickle.
 *
 */

template<class Force>
std::map<std::string, pybind11::tuple> pybind11_getstate_force_class_factory(
    ClassFactory<Force> const& classFactory) {
    std::map<std::string, pybind11::tuple> stateMap;
    for (auto it=classFactory.cbegin(); it != classFactory.cend(); ++it)
        { stateMap.emplace(it->first, (it->second)->pybind11_getstate()); }
    return stateMap;
}

// template specialisation
template std::map<std::string, pybind11::tuple>
pybind11_getstate_force_class_factory<HalfEdgeForce<ForcesType>>
    (ClassFactory<HalfEdgeForce<ForcesType>> const&);
template std::map<std::string, pybind11::tuple>
pybind11_getstate_force_class_factory<VertexForce<ForcesType>>
    (ClassFactory<VertexForce<ForcesType>> const&);

template<class Force>
void pybind11_setstate_force_class_factory(
    VertexModel& vm, std::map<std::string, pybind11::tuple> const& stateMap);

#endif

