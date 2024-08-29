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
    VolumeForce,
    // argument types
    double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& kV, double const& H0, double const& A0) {
    // set force
    vertexForces.add<VolumeForce>(
        name,                       // (unique) user-defined name for forcess
        kV, H0, A0,                 // user-defined parameters
        this, &forces, &vertices);  // VertexModel attributes
}

template<> void VertexModel::addVertexForce<
    // derived force class
    LinearVolumeForce,
    // argument types
    double const&, double const&,
    double const&, double const&,
    double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& kA, double const& A0,
    double const& kP, double const& P0,
    double const& taur, double const& H0, double const& taua) {
    // set force
    vertexForces.add<LinearVolumeForce>(
        name,                           // (unique) user-defined name for forcess
        kA, A0, kP, P0, taur, H0, taua, // user-defined parameters
        this, &forces, &vertices);      // VertexModel attributes
}

template<> void VertexModel::addVertexForce<
    // derived force class
    PressureForce,
    // argument types
    double const&, bool const&>(
    // user-defined arguments
    std::string const& name,
    double const& F, bool const& fixedPerimeterForce) {
    // set force
    vertexForces.add<PressureForce>(
        name,                       // (unique) user-defined name for forcess
        F, fixedPerimeterForce,     // user-defined parameters
        this, &forces, &vertices);  // VertexModel attributes
}

template<> void VertexModel::addVertexForce<
    // derived force class
    BoundaryTension,
    // argument types
    double const&>(
    // user-defined arguments
    std::string const& name, double const& gamma) {
    // set force
    vertexForces.add<BoundaryTension>(
        name,                       // (unique) user-defined name for forcess
        gamma,                      // user-defined parameters
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

template<> void VertexModel::addHalfEdgeForce<
    // derived force class
    Model3,
    // argument types
    double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& Gamma,
    double const& sigma, double const& taup) {
    // set force
    halfEdgeForces.add<Model3>(
        name,                                           // (unique) user-defined name for forces
        Gamma, sigma, taup,                             // user-defined parameters
        this, &random, &forces, &halfEdges);            // VertexModel attributes
}

template<> void VertexModel::addHalfEdgeForce<
    // derived force class
    Model4,
    // argument types
    double const&, double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& Gamma, double const& taur,
    double const& sigma, double const& taup) {
    // set force
    halfEdgeForces.add<Model4>(
        name,                                           // (unique) user-defined name for forces
        Gamma, taur, sigma, taup,                       // user-defined parameters
        this, &random, &forces, &halfEdges);            // VertexModel attributes
}

/*
 *  KERATIN
 *
 */

template<> void VertexModel::addVertexForce<    // initial time as argument (for copies)
    // derived force class
    KeratinModel,
    // argument types
    double const&, double const&, double const&,
    double const&, double const&, double const&,
    double const&, double const&, double const&,
    double const&, double const&, double const&>(
    // user-defined arguments
    std::string const& name,
    double const& K, double const& A0, double const& taur,
    double const& Gamma, double const& p0, double const& T,
    double const& alpha, double const& beta, double const& kth,
    double const& tau, double const& sigma, double const& ron) {
    // set force
    vertexForces.add<KeratinModel>(
        name,                                                           // (unique) user-defined name for forces
        K, A0, taur, Gamma, p0, T, alpha, beta, kth, tau, sigma, ron,   // user-defined parameters
        this, &random, &forces, &vertices);                             // VertexModel attributes and initial time
}

