/*
Provide pickling ability for C++ force computation objects via pybind11.
*/

#ifndef FORCES_PICKLE_HPP
#define FORCES_PICKLE_HPP

#include <pybind11/pybind11.h>                                                  
#include <pybind11/stl.h>

#include "forces.hpp"
#include "pickle.hpp"
#include "system.hpp"

/*
 *  PerimeterForce
 *
 */

template<> pybind11::tuple pybind11_getstate_forces<
    VertexForce<ForcesType>, PerimeterForce>(
        std::shared_ptr<VertexForce<ForcesType>> const& obj);               // use default behaviour

template<> void pybind11_setstate_forces<
    VertexForce<ForcesType>, PerimeterForce>(
        pybind11::tuple const& t, VertexModel& vm);                         // use default behaviour

template<>
void VertexModel::addVertexForce<PerimeterForce, ParametersType const&>(    // forces binding to VertexModel from parameters
    std::string const& name, ParametersType const& parameters) {
    addVertexForce<PerimeterForce, double const&, double const&>(
        name, parameters.at("kP"), parameters.at("P0"));
}

/*
 *  AreaForce
 *
 */

template<> pybind11::tuple pybind11_getstate_forces<
    VertexForce<ForcesType>, AreaForce>(
        std::shared_ptr<VertexForce<ForcesType>> const& obj);               // use default behaviour

template<> void pybind11_setstate_forces<
    VertexForce<ForcesType>, AreaForce>(
        pybind11::tuple const& t, VertexModel& vm);                         // use default behaviour

template<>
void VertexModel::addVertexForce<AreaForce, ParametersType const&>(         // forces binding to VertexModel from parameters
    std::string const& name, ParametersType const& parameters) {
    addVertexForce<AreaForce, double const&, double const&>(
        name, parameters.at("kA"), parameters.at("A0"));
}

/*
 *  ActiveBrownianForce
 *
 */

template<> pybind11::tuple pybind11_getstate_forces<
    VertexForce<ForcesType>, ActiveBrownianForce>(
        std::shared_ptr<VertexForce<ForcesType>> const& obj) {
    return pybind11::make_tuple(
        obj->getType(), obj->getParameters(), obj->getTheta());
}

template<> void pybind11_setstate_forces<
    VertexForce<ForcesType>, ActiveBrownianForce>(
        pybind11::tuple const& t, VertexModel& vm) {
    // check
    checkSize(t, 3);
    // get data
    std::string const type =
        t[0].cast<std::string>();
    std::map<std::string, double> const parameters =
        t[1].cast<std::string>();
    std::map<std::string, double> const theta =
        t[2].cast<std::string>();
    // add force
    vm.addVertexForce<DerivedType, ParametersType const&,
        std::map<std::string, double> const&>(
        type, parameters, theta);
}

template<>
void VertexModel::addVertexForce<ActiveBrownianForce, ParametersType const&,    // forces binding to VertexModel from parameters
    std::map<std::string, double> const&>(
    std::string const& name, ParametersType const& parameters,
        std::map<std::string, double> const& theta) {
    addVertexForce<AreaForce, double const&, double const&>(
        name, parameters.at("v0"), parameters.at("taup"));
    halfEdgeForces.get(name).setTheta(theta);
    halfEdgeForces.get(name).integrate(0);
}

#endif

