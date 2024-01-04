/*
Forces definitions.
*/

#ifndef FORCES_HPP
#define FORCES_HPP

#include <math.h>
#include <numbers>
#include <iostream>

#include "base_forces.hpp"
#include "mesh.hpp"
#include "random.hpp"
#include "tools.hpp"

typedef std::map<long int, std::vector<double>> ForcesType;

class PerimeterForce : public VertexForce<ForcesType> {
/*
Cell perimeter restoring force.
*/

    protected:

        Mesh* const mesh;

    public:

        PerimeterForce(
            double const& kP_, double const& P0_,
            Mesh* mesh_, ForcesType* forces_, VerticesType* vertices_) :
            VertexForce<ForcesType>("centre",
                {{"kP", kP_}, {"P0", P0_}},
                forces_, vertices_),
            mesh(mesh_) {}

        void addForce(Vertex const& vertex) override {
//             std::cerr <<
//                 "Perimeter force on vertex " << vertex.getIndex() << "."
//                 << std::endl;
            double const perimeter = mesh->getVertexToNeighboursPerimeter(
                vertex.getIndex());
            std::vector<long int> neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            int numberNeighbours = neighbourVerticesIndices.size();
            std::vector<double> toPreviousNeighbour, toNextNeighbour;
            for (int i=0; i < numberNeighbours; i++) {
                assert(vertices->at(neighbourVerticesIndices[i]).getType()
                    != "centre");
                toPreviousNeighbour = mesh->wrapTo(
                    neighbourVerticesIndices[i],
                    neighbourVerticesIndices[pmod(i - 1, numberNeighbours)],
                    true);
                toNextNeighbour = mesh->wrapTo(
                    neighbourVerticesIndices[i],
                    neighbourVerticesIndices[pmod(i + 1, numberNeighbours)],
                    true);
                for (int dim=0; dim < 2; dim++) {
                    (*forces)[neighbourVerticesIndices[i]][dim] +=
                        parameters.at("kP")*(
                            perimeter - parameters.at("P0"))*(
                            toPreviousNeighbour[dim] + toNextNeighbour[dim]);
                }
            }
        }

};

class AreaForce : public VertexForce<ForcesType> {
/*
Cell area restoring force.
*/

    protected:

        Mesh* mesh;

    public:

        AreaForce(
            double const& kA_, double const& A0_,
            Mesh* mesh_, ForcesType* forces_, VerticesType* vertices_) :
            VertexForce<ForcesType>("centre",
                {{"kA", kA_}, {"A0", A0_}},
                forces_, vertices_),
            mesh(mesh_) {}

        void addForce(Vertex const& vertex) override {
//             std::cerr <<
//                 "Area force on vertex " << vertex.getIndex() << "."
//                 << std::endl;
            double const area = mesh->getVertexToNeighboursArea(
                vertex.getIndex());
            std::vector<long int> neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            int numberNeighbours = neighbourVerticesIndices.size();
            std::vector<double> crossToPreviousNeighbour, crossToNextNeighbour;
            for (int i=0; i < numberNeighbours; i++) {
                assert(vertices->at(neighbourVerticesIndices[i]).getType()
                    != "centre");
                crossToPreviousNeighbour = cross2z(mesh->wrapTo(
                    neighbourVerticesIndices[i],
                    neighbourVerticesIndices[pmod(i - 1, numberNeighbours)],
                    false));
                crossToNextNeighbour = cross2z(mesh->wrapTo(
                    neighbourVerticesIndices[i],
                    neighbourVerticesIndices[pmod(i + 1, numberNeighbours)],
                    false));
                for (int dim=0; dim < 2; dim++) {
                    (*forces)[neighbourVerticesIndices[i]][dim] +=
                        (parameters.at("kA")/2.)*(
                            area - parameters.at("A0"))*(
                            crossToPreviousNeighbour[dim]
                                - crossToNextNeighbour[dim]);
                }
            }
        }

};

class ActiveBrownianForce : public VertexForce<ForcesType> {
/*
Active Brownian self-propulsion force acting on vertices.
*/

    protected:

        Random* random;                                     // random number generator
        std::map<long int, double> theta;                   // orientation of self-propulsion force
        std::map<long int, std::vector<double>> abForces;   // actual forces
        bool const woCOM = true;                            // substract centre of mass force

    public:

        ActiveBrownianForce(
            double const& v0_, double const& taup_,
            Random* random_, ForcesType* forces_, VerticesType* vertices_) :
            VertexForce<ForcesType>("vertex",
                {{"v0", v0_}, {"taup", taup_}},
                forces_, vertices_),
            random(random_)
            { integrate(0); }

        std::map<long int, double> const& getTheta() const
            { return theta; }
        std::map<long int, std::vector<double>> const& getAbForces() const
            { return abForces; }

        void setTheta(std::map<long int, double> const& theta_)
            { theta = theta_; }

        void addForce(Vertex const& vertex) override {
//             std::cerr <<
//                 "Active Brownian force on vertex " << vertex.getIndex() << "."
//                 << std::endl;
            for (int dim=0; dim < 2; dim++) {
                (*forces)[vertex.getIndex()][dim] +=
                    abForces[vertex.getIndex()][dim];
            }
        }

        void integrate(double const& dt) override {

            // index correspondence between vertices, orientations, and forces
            for (auto it=theta.begin(); it != theta.end(); ++it) {
                if (!inMap(*vertices, it->first)) { // vertex index not present anymore
                    theta.erase(it->first);
                    abForces.erase(it->first);
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(theta, it->first)        // new vertex index
                    && (it->second).getType() == type) {
                    theta[it->first] = 2.*std::numbers::pi*random->random01();
                    abForces[it->first] = {0, 0};
                }
            }

            // force computation
            double avForce[2] = {0, 0};
            for (auto it=theta.begin(); it != theta.end(); ++it) {
                theta[it->first] +=
                    sqrt(2.*dt/parameters.at("taup"))*random->gauss();
                for (int dim=0; dim < 2; dim++) {
                    abForces[it->first][dim] =  // active Brownian force
                        parameters.at("v0")*cos(theta[it->first]
                            - dim*std::numbers::pi/2.);
                    if (woCOM) {                // centre of mass force
                        avForce[dim] += abForces[it->first][dim];
                    }
                }
            }
            if (woCOM) {                        // remove centre of mass force
                for (int dim=0; dim < 2; dim++) {
                    avForce[dim] /= theta.size();
                    for (auto it=theta.begin(); it != theta.end(); ++it) {
                        abForces[it->first][dim] -= avForce[dim];
                    }
                }
            }
        }

};

#endif

