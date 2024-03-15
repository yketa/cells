/*
Integrators definitions.
*/

#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <assert.h>
#include <map>
#include <vector>

#include <Eigen/Sparse>

#include "base_integrators.hpp"
#include "mesh.hpp"

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

class PairFriction : public BaseIntegrator<ForcesType, VelocitiesType> {
/*
Cell corner pair friction and unit vertex-substrate drag coefficient overdamped
integrator.

WARNING: This explicitly ignores centre vertices.
*/

    protected:

        Mesh* const mesh;   // mesh object

        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> sol;      // solver for inversion

    public:

        PairFriction(
            double const& eta_,
            Mesh* const& mesh_,
            ForcesType* const forces_, VelocitiesType* velocities_) :
            BaseIntegrator<ForcesType, VelocitiesType>(
                {{"eta", eta_}},
                forces_, velocities_),
            mesh(mesh_) {}

        void integrate(double const& dt) override {

            /*
             *  INTERPARTICLE FRICTION
             *  ----------------------
             *  \dot{X} = velocity, F = force, \underline{C} = contact matrix
             *  C_{iajb} = -\delta_{ab},                          if i and j are in contact
             *             \delta_{ab} {number of contacts of i}, if i == j
             *             0,                                     otherwise
             *  \dot{X} = F - \xi \underline{C} \dot{X}
             *      <=> \dot{X} = \underline{M}^{-1} F
             *  \underline{M} = modified contact matrix
             *  \underline{M} = \mathbbm{1} + \xi \underline{C}
             */

            // IDENTIFY VERTICES

            long int count = 0;
            std::map<long int, long int> vertexIndices; // associate contiguous indices to vertices
            std::vector<long int> invVertexIndices;     // inverse association of vertex indices
            std::vector<long int> sumContacts;          // number of contact per vertex
            VerticesType const& vertices = mesh->getVertices();
            for (auto it=vertices.begin(); it != vertices.end(); ++it) {
                if (
                    (it->second).getType() != "centre"  // not center vertices
                    && !(it->second).getBoundary()) {   // not boundary vertices
                    vertexIndices.emplace(it->first, count);    // forward map
                    invVertexIndices.push_back(it->first);      // backward map
                    sumContacts.push_back(0);
                    count++;
                }
            }
            long int const N = count;                   // number of vertices

            // COMPUTE CONTACT MATRIX

            std::vector<Eigen::Triplet<double>> modContactTrip(0);              // modified contact matrix elements
            Eigen::SparseMatrix<double, Eigen::RowMajor> modContact(2*N, 2*N);  // modified contact matrix
            Eigen::VectorXd generalVelocities(2*N);                             // general velocity vector
            Eigen::VectorXd generalForces(2*N);                                 // general force vector

            HalfEdgesType const& halfEdges = mesh->getHalfEdges();
            for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {      // loop over all half-edges therefore each pair of connected vertices will appear twice in forward and backward direction
                Vertex const v0 = vertices.at((it->second).getFromIndex());     // origin vertex
                Vertex const v1 = vertices.at((it->second).getToIndex());       // destination vertex
                if (
                    v0.getType() != "centre" && v1.getType() != "centre"        // not centre vertices
                    && !v0.getBoundary() && !v1.getBoundary()) {                // not boundary vertices
                    long int iv0 = vertexIndices[v0.getIndex()];
                    long int iv1 = vertexIndices[v1.getIndex()];

                    for (int dim=0; dim < 2; dim++) {   // off-diagonal terms
                        modContactTrip.push_back(
                            Eigen::Triplet<double>(
                                2*iv0 + dim, 2*iv1 + dim,
                                -parameters.at("eta")));
                    }

                    sumContacts[iv0]++;                 // compute number of contacts
                }
            }

            for (long int i=0; i < N; i++) {            // loop over all vertices
                for (int dim=0; dim < 2; dim++) {       // diagonal terms
                    modContactTrip.push_back(
                        Eigen::Triplet<double>(
                            2*i + dim, 2*i + dim,
                            1 + parameters.at("eta")*sumContacts[i]));
                }
            }

            // SOLVE LINEAR SYSTEM FOR VELOCITIES

            modContact.setZero();                               // reset modified contact matrix
            modContact.setFromTriplets(                         // set modified contact matrix
                modContactTrip.begin(), modContactTrip.end());

            for (long int i=0; i < N; i++) {
                long int const vertexIndex = invVertexIndices[i];
                for (int dim=0; dim < 2; dim++) {
                    generalForces[2*i + dim] =          // set (general) forces vector
                        (forces->at(vertexIndex)).at(dim);
                }
            }

            sol.compute(modContact);                            // compute decomposition
            assert(sol.info() == Eigen::Success);               // check decomposition worked
            generalVelocities = sol.solve(generalForces);       // invert matrix to compute velocities
            assert(sol.info() == Eigen::Success);               // check inversion worked

            velocities->clear();
            for (long int i=0; i < N; i++) {
                long int const vertexIndex = invVertexIndices[i];
                velocities->emplace(vertexIndex, std::vector<double>({0, 0}));
                for (int dim=0; dim < 2; dim++) {
                    (*velocities)[vertexIndex][dim] =   // set velocities vector
                        generalVelocities[2*i + dim];
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

#endif

