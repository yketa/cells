/*
Integrators definitions.
*/

#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include "base_integrators.hpp"

class UnitOverdamped : public BaseIntegrator<ForcesType, VelocitiesType> {
/*
Unit vertex-substrate drag coefficient overdamped integrator.
*/

    public:

        UnitOverdamped(
            ForcesType* const forces_, VelocitiesType* velocities_) :
            BaseIntegrator<ForcesType, VelocitiesType>(
                {}, forces_, velocities_) {}

        void integrate(double const& dt) override { *velocities = *forces; }

        pybind11::tuple pybind11_getstate() const override;

};

#endif

