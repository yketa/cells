/*
Provide pickling ability for C++ integrator objects via pybind11.
*/

#ifndef INTEGRATORS_PICKLE_HPP
#define INTEGRATORS_PICKLE_HPP

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base_pickle.hpp"
#include "integrators.hpp"
#include "system.hpp"

/*
 *  UnitOverdamped
 *
 */

// save state
pybind11::tuple UnitOverdamped::pybind11_getstate() const {
    return pybind11::make_tuple(
        // unique identifying string for the integrator object
        "UnitOverdamped"
        // state
        );
}

// load state
template<> void
VertexModel::setIntegrator<UnitOverdamped, pybind11::tuple const&>(
    pybind11::tuple const& t) {
    // check
    checkSize(t, 1);
    assert(t[0].cast<std::string>() == "UnitOverdamped");
    // initialise force
    setIntegrator<UnitOverdamped>();
}

/*
 *  pickle -> VertexModel::setIntegrator
 *
 */

void pybind11_setstate_integrator(VertexModel& vm, pybind11::tuple const& t) {
    std::string const integratorName = t[0].cast<std::string>();
    // templated lambda function to set integrator in VertexModel
    auto setIntegrator = [&vm, &t]<class Integrator>()
        { vm.setIntegrator<Integrator, pybind11::tuple const&>(t); };
    // --- MANUAL ASSOCIATION TO INTEGRATORS ---
    if (integratorName == "UnitOverdamped") {
        setIntegrator.template operator()<UnitOverdamped>();
    }
    // throw error if integrator not recognised
    else {
            throw std::runtime_error(
                "Integrator object '" + integratorName + "' does not exist.");
    }
}

#endif

