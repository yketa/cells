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

        Random* random;                     // random number generator
        std::map<long int, double> theta;   // orientation of self-propulsion force

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
        void setTheta(std::map<long int, double> const& theta_)
            { theta = theta_; }

        void addForce(Vertex const& vertex) override {
//             std::cerr <<
//                 "Active Brownian force on vertex " << vertex.getIndex() << "."
//                 << std::endl;
            for (int dim=0; dim < 2; dim++) {
                (*forces)[vertex.getIndex()][dim] +=
                    parameters.at("v0")*cos(theta[vertex.getIndex()]
                        - dim*std::numbers::pi/2.);
            }
        }

        void integrate(double const& dt) override {
            // index correspondence between vertices and orientations
            for (auto it=theta.begin(); it != theta.end(); ++it) {
                if (!inMap(*vertices, it->first)) { // vertex index not present anymore
                    theta.erase(it->first);
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(theta, it->first)        // new vertex index
                    && (it->second).getType() == type) {
                    theta[it->first] = 2.*std::numbers::pi*random->random01();
                }
            }
            // integration
            double const amp = sqrt(2.*dt/parameters.at("taup"));
            for (auto it=theta.begin(); it != theta.end(); ++it) {
                theta[it->first] += amp*random->gauss();
            }
        }

};

class OrnsteinUhlenbeckTension : public HalfEdgeForce<ForcesType> {
/*
Ornstein-Uhlenbeck tension acting on junctions.
http://arxiv.org/abs/2309.04818
*/

    protected:

        Mesh* mesh;                         // mesh object
        Random* random;                     // random number generator
        std::map<long int, double> tension; // orientation of self-propulsion force

    public:

        OrnsteinUhlenbeckTension(
            double const& t0_, double const& st0_, double const& taup_,
            Mesh* mesh_, Random* random_,
            ForcesType* forces_, HalfEdgesType* halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"t0", t0_}, {"st0", st0_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getTension() const
            { return tension; }
        void setTension(std::map<long int, double> const& tension_)
            { tension = tension_; }

        void addForce(HalfEdge const& halfEdge) override {
            long int const fromIndex = halfEdge.getFromIndex();     // origin vertex
            long int const toIndex = halfEdge.getToIndex();         // destination vertex
            std::vector<double> const fromTo =
                mesh->getHalfEdgeVector(halfEdge.getIndex(), true); // unit half-edge vector
            double force;
            for (int dim=0; dim < 2; dim++) {
                force = tension[halfEdge.getIndex()]*fromTo[dim];
                (*forces)[fromIndex][dim] += force;
                (*forces)[toIndex][dim] -= force;
            }
        }

        void integrate(double const& dt) override {
            // index correspondence between half-edges and tensions
            for (auto it=tension.begin(); it != tension.end(); ++it) {
                if (!inMap(*halfEdges, it->first)) {    // half-edge index not present anymore
                    tension.erase(it->first);
                }
            }
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if (!inMap(tension, it->first)          // new half-edge index
                    && (it->second).getType() == type) {
                    tension[it->first] = parameters.at("t0")
                        + parameters.at("st0")*random->gauss();
                }
            }
            // integration
            double const dt_ = dt/parameters.at("taup");
            double const amp = parameters.at("st0")*sqrt(2.*dt_);
            for (auto it=tension.begin(); it != tension.end(); ++it) {
                tension[it->first] =
                    (1 - dt_)*tension[it->first]
                        + dt_*parameters.at("t0")
                        + amp*random.gauss();
            }
        }

};

#endif

