/*
Integrators binding to VertexModel (system.hpp).
*/

#include "integrators.hpp"
#include "system.hpp"

template<> void VertexModel::setIntegrator<
    // derived integrator class
    UnitOverdamped
    // argument types
    >(
    // user-defined arguments
    ) {
    // set integrator
    integrator = std::make_shared<UnitOverdamped>(
                                // user-defined parameters
        &forces, &velocities);  // VertexModel attributes
}

template<> void VertexModel::setIntegrator<
    // derived integrator class
    PairFriction,
    // argument types
    double const&>(
    // user-defined arguments
    double const& eta) {
    // set integrator
    integrator = std::make_shared<PairFriction>(
        eta,                            // user-defined parameters
        this, &forces, &velocities);    // VertexModel attributes
}

