/*
Base types for integrators.
*/

#ifndef BASE_INTEGRATORS_HPP
#define BASE_INTEGRATORS_HPP

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base.hpp"

// GENERAL BASE

template<class T, class TT>
class BaseIntegrator {
/*
All integrators.
*/

    protected:

        ParametersType const parameters;    // force parameters
        T* const forces;                    // pointer to forces container
        TT* velocities;                     // pointer to velocities container

    public:

        // GETTERS

        ParametersType const getParameters() const
            { return parameters; }
        T* const& getForces() const
            { return forces; }
        TT* getVelocities()
            { return velocities; }

        // CONSTRUCTORS AND DESTRUCTORS

        BaseIntegrator(
            ParametersType const& parameters_,
            T* const& forces_, TT* velocities_) :
            parameters(parameters_),
            forces(forces_), velocities(velocities_) {}

        BaseIntegrator(BaseIntegrator<T, TT>& bI) : // copy constructor
            BaseIntegrator(
                bI.getParameters(), bI.getForces(), bI.getVelocities()) {}

        virtual ~BaseIntegrator() {}

        // METHODS

        virtual void integrate(double const& dt) {}
        /*
        Integrate velocities.
        */

        virtual pybind11::tuple pybind11_getstate() const
            { return pybind11::make_tuple(); }
        /*
        Tuple with data to regenerate force computation object.
        First element MUST be an unique identifying string for the force object
        type.
        */

};

#endif

