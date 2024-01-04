/*
Forces binding to VertexModel (system.hpp).
*/

#include "forces.hpp"
#include "system.hpp"

template<> void VertexModel::addVertexForce<
    // derived force class
    PerimeterForce,
    // argument types
    double const&, double const&>(
    // user-defined arguments
    std::string const& name, double const& kP, double const& P0) {
    // set force
    vertexForces.add<PerimeterForce>(
        name,                       // (unique) user-defined name for forcess
        kP, P0,                     // user-defined parameters
        this, &forces, &vertices);  // VertexModel attributes
}

template<> void VertexModel::addVertexForce<
    // derived force class
    AreaForce,
    // argument types
    double const&, double const&>(
    // user-defined arguments
    std::string const& name, double const& kA, double const& A0) {
    // set force
    vertexForces.add<AreaForce>(
        name,                       // (unique) user-defined name for forcess
        kA, A0,                     // user-defined parameters
        this, &forces, &vertices);  // VertexModel attributes
}

template<> void VertexModel::addVertexForce<
    // derived force class
    ActiveBrownianForce,
    // argument types
    double const&, double const&>(
    // user-defined arguments
    std::string const& name, double const& v0, double const& taup) {
    // set force
    vertexForces.add<ActiveBrownianForce>(
        name,                           // (unique) user-defined name for forcess
        v0, taup,                       // user-defined parameters
        &random, &forces, &vertices);   // VertexModel attributes
}

