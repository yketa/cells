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

template<> void VertexModel::addHalfEdgeForce<
    // derived force class
    OrnsteinUhlenbeckTension,
    // argument types
    double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& t0, double const& st0, double const& taup) {
    // set force
    halfEdgeForces.add<OrnsteinUhlenbeckTension>(
        name,                                   // (unique) user-defined name for forcess
        t0, st0, taup,                          // user-defined parameters
        this, &random, &forces, &halfEdges);    // VertexModel attributes
}

/*
 *  MODELS 0-4
 *
 */

template<> void VertexModel::addHalfEdgeForce<
    // derived force class
    Model0,
    // argument types
    double const&, double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& Gamma, double const& P0,
    double const& sigma, double const& taup) {
    // set force
    halfEdgeForces.add<Model0>(
        name,                                           // (unique) user-defined name for forces
        Gamma, P0, sigma, taup,                         // user-defined parameters
        this, &vertices, &random, &forces, &halfEdges); // VertexModel attributes
}

template<> void VertexModel::addHalfEdgeForce<
    // derived force class
    Model1,
    // argument types
    double const&, double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& Gamma, double const& P0,
    double const& sigma, double const& taup) {
    // set force
    halfEdgeForces.add<Model1>(
        name,                                           // (unique) user-defined name for forces
        Gamma, P0, sigma, taup,                         // user-defined parameters
        this, &vertices, &random, &forces, &halfEdges); // VertexModel attributes
}

template<> void VertexModel::addHalfEdgeForce<
    // derived force class
    Model2,
    // argument types
    double const&, double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& Gamma, double const& taur,
    double const& sigma, double const& taup) {
    // set force
    halfEdgeForces.add<Model2>(
        name,                                           // (unique) user-defined name for forces
        Gamma, taur, sigma, taup,                       // user-defined parameters
        this, &random, &forces, &halfEdges);            // VertexModel attributes
}

