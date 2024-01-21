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
            Mesh* const mesh_, ForcesType* forces_, VerticesType* vertices_) :
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

        Mesh* const mesh;

    public:

        AreaForce(
            double const& kA_, double const& A0_,
            Mesh* const mesh_, ForcesType* forces_, VerticesType* vertices_) :
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
            for (auto it=theta.begin(); it != theta.end();) {
                if (!inMap(*vertices, it->first)) {                     // vertex index not present any more
                    it = theta.erase(it);
                }
                else if ((vertices->at(it->first)).getType() != type) { // not a vertex any more
                    it = theta.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(theta, it->first)                            // new vertex index
                    && (it->second).getType() == type) {
                    theta[it->first] = 2.*std::numbers::pi*random->random01();
                }
            }
            // integration
            double const amp = sqrt(2.*dt/parameters.at("taup"));
            for (auto it=theta.begin(); it != theta.end(); ++it) {
                it->second += amp*random->gauss();
            }
        }

};

class OrnsteinUhlenbeckTension : public HalfEdgeForce<ForcesType> {
/*
Ornstein-Uhlenbeck tension acting on junctions.
http://arxiv.org/abs/2309.04818
*/

    protected:

        Mesh* const mesh;                   // mesh object
        Random* random;                     // random number generator
        std::map<long int, double> tension; // tension in each junction

    public:

        OrnsteinUhlenbeckTension(
            double const& t0_, double const& st0_, double const& taup_,
            Mesh* const mesh_, Random* random_,
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
            for (auto it=tension.begin(); it != tension.end();) {
                if (!inMap(*halfEdges, it->first)) {                        // half-edge index not present any more
                    it = tension.erase(it);
                }
                else if ((halfEdges->at(it->first)).getType() != type) {    // half-edge is not a junction any more
                    it = tension.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if (!inMap(tension, it->first)                              // new half-edge index
                    && (it->second).getType() == type) {
                    tension[it->first] = parameters.at("t0")
                        + parameters.at("st0")*random->gauss();
                }
            }
            // integration
            double const dt_ = dt/parameters.at("taup");
            double const amp = parameters.at("st0")*sqrt(2.*dt_);
            for (auto it=tension.begin(); it != tension.end(); ++it) {
                it->second =
                    (1 - dt_)*tension[it->first]
                        + dt_*parameters.at("t0")
                        + amp*random->gauss();
            }
        }

};

/*
 *  MODELS 0-4
 *
 */

class Model0 : public HalfEdgeForce<ForcesType> {

    protected:

        Mesh* const mesh;                       // mesh object
        VerticesType* const vertices;           // vertices
        Random* random;                         // random number generator
        std::map<long int, double> perimeter;   // cell perimeters
        std::map<long int, double> noise;       // noise in tension
        std::map<long int, double> tension;     // junction tension

    public:

        Model0(
            double const& Gamma_, double const& P0_,
            double const& sigma_, double const& taup_,
            Mesh* const mesh_, VerticesType* const vertices_, Random* random_,
            ForcesType* forces_, HalfEdgesType* halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_}, {"P0", P0_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), vertices(vertices_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getPerimeter() const
            { return perimeter; }

        std::map<long int, double> const& getNoise() const
            { return noise; }
        void setNoise(std::map<long int, double> const& noise_)
            { noise = noise_; }

        std::map<long int, double> const& getTension() const
            { return tension; }

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
            // index correspondence between cell centre vertices and perimeters
            perimeter.clear();
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == "centre") {                   // compute perimeter around cell centres
                    perimeter[it->first] =
                        mesh->getVertexToNeighboursPerimeter(it->first);
                }
            }
            // index correspondence between half-edges and noises
            for (auto it=noise.begin(); it != noise.end();) {
                if (!inMap(*halfEdges, it->first)) {                        // half-edge index not present any more
                    it = noise.erase(it);
                }
                else if ((halfEdges->at(it->first)).getType() != type) {    // half-edge is not a junction any more
                    it = noise.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if (!inMap(noise, it->first)                                // new half-edge index
                    && (it->second).getType() == type) {
                    noise[it->first] = parameters.at("sigma")*random->gauss();
                }
            }
            // clear tension
            tension.clear();
            // noise integration and tension setting
            double const dt_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dt_);
            for (auto it=noise.begin(); it != noise.end(); ++it) {
                // noise
                noise[it->first] =
                    (1 - dt_)*noise[it->first] + amp*random->gauss();
                // tension
                long int const cellIndices[2] =                             // indices of cell of the junction
                    {halfEdges->at(
                        halfEdges->at(it->first)
                            .getNextIndex()).getToIndex(),
                    halfEdges->at(
                        halfEdges->at(halfEdges->at(it->first).getPairIndex())
                            .getNextIndex()).getToIndex()};
                tension[it->first] = noise[it->first];                      // stochastic part
                for (long int cellIndex : cellIndices) {
                    if (inMap(perimeter, cellIndex)) {
                        tension[it->first] += parameters.at("Gamma")*(      // deterministic part
                            perimeter[cellIndex] - parameters.at("P0"));
                    }
                }
            }
        }

};

class Model1 : public HalfEdgeForce<ForcesType> {

    protected:

        Mesh* const mesh;                       // mesh object
        VerticesType* const vertices;           // vertices
        Random* random;                         // random number generator
        std::map<long int, double> perimeter;   // cell perimeters
        std::map<long int, double> tension;     // junction tension

    public:

        Model1(
            double const& Gamma_, double const& P0_,
            double const& sigma_, double const& taup_,
            Mesh* const mesh_, VerticesType* const vertices_, Random* random_,
            ForcesType* forces_, HalfEdgesType* halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_}, {"P0", P0_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), vertices(vertices_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getPerimeter() const
            { return perimeter; }

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
            // index correspondence between cell centre vertices and perimeters
            perimeter.clear();
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == "centre") {                   // compute perimeter around cell centres
                    perimeter[it->first] =
                        mesh->getVertexToNeighboursPerimeter(it->first);
                }
            }
            // index correspondence between half-edges and tensions
            for (auto it=tension.begin(); it != tension.end();) {
                if (!inMap(*halfEdges, it->first)) {                        // half-edge index not present any more
                    it = tension.erase(it);
                }
                else if ((halfEdges->at(it->first)).getType() != type) {    // half-edge is not a junction any more
                    it = tension.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if (!inMap(tension, it->first)                              // new half-edge index
                    && (it->second).getType() == type) {
                    tension[it->first] =
                        parameters.at("sigma")*random->gauss();             // noise
                    long int const cellIndices[2] =                         // indices of cell of the junction
                        {halfEdges->at(
                            halfEdges->at(it->first)
                                .getNextIndex()).getToIndex(),
                        halfEdges->at(
                            halfEdges->at(halfEdges->at(it->first).getPairIndex())
                                .getNextIndex()).getToIndex()};
                    for (long int cellIndex : cellIndices) {
                        if (inMap(perimeter, cellIndex)) {
                            tension[it->first] += parameters.at("Gamma")*(  // target tension
                                perimeter[cellIndex] - parameters.at("P0"));
                        }
                    }
                }
            }
            // tension integration
            double const dt_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dt_);
            for (auto it=tension.begin(); it != tension.end(); ++it) {
                // tension
                long int const cellIndices[2] =                             // indices of cell of the junction
                    {halfEdges->at(
                        halfEdges->at(it->first)
                            .getNextIndex()).getToIndex(),
                    halfEdges->at(
                        halfEdges->at(halfEdges->at(it->first).getPairIndex())
                            .getNextIndex()).getToIndex()};
                tension[it->first] = (1 - dt_)*tension[it->first]           // zero-mean deterministic part
                    + amp*random->gauss();                                  // stochastic part
                for (long int cellIndex : cellIndices) {
                    if (inMap(perimeter, cellIndex)) {
                        tension[it->first] += dt_*parameters.at("Gamma")*(  // target deterministic part
                            perimeter[cellIndex] - parameters.at("P0"));
                    }
                }
            }
        }

};

#endif

