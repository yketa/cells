/*
Forces definitions.
*/

#ifndef FORCES_HPP
#define FORCES_HPP

#include <cmath>
#include <numbers>
#include <iostream>

#include "base_forces.hpp"
#include "mesh.hpp"
#include "random.hpp"
#include "tools.hpp"

class PerimeterForce : public VertexForce<ForcesType> {
/*
Cell perimeter restoring force.
*/

    protected:

        Mesh* const mesh;

    public:

        PerimeterForce(
            double const& kP_, double const& P0_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
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
            long int const numberNeighbours = neighbourVerticesIndices.size();
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

        pybind11::tuple pybind11_getstate() const override;

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
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
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
            long int const numberNeighbours = neighbourVerticesIndices.size();
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

        pybind11::tuple pybind11_getstate() const override;

};

class VolumeForce : public VertexForce<ForcesType> {
/*
Cell volume restoring force.
*/

    protected:

        Mesh* const mesh;                           // mesh object
        std::map<long int, double> height;          // height of each cell

        std::map<long int, double> heightVelocity;  // height velocity of each cell

    public:

        VolumeForce(
            double const& kV_, double const& h0_, double const& A0_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"kV", kV_}, {"h0", h0_}, {"A0", A0_}},
                forces_, vertices_),
            mesh(mesh_)
            { integrate(0); }

        std::map<long int, double> const& getHeight() const
            { return height; }
        void setHeight(std::map<long int, double> const& height_)
            { height = height_; }

        std::map<long int, double> const& getHeightVelocity() const
            { return heightVelocity; }

        void addForce(Vertex const& vertex) override {
            double const volume = height[vertex.getIndex()]*
                mesh->getVertexToNeighboursArea(vertex.getIndex());
            std::vector<long int> neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            long int const numberNeighbours = neighbourVerticesIndices.size();
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
                        (parameters.at("kV")/2.)*height[vertex.getIndex()]*(
                            volume - parameters.at("h0")*parameters.at("A0"))*(
                            crossToPreviousNeighbour[dim]
                                - crossToNextNeighbour[dim]);
                }
            }
            heightVelocity[vertex.getIndex()] += -parameters.at("kV")*
                (volume - parameters.at("h0")*parameters.at("A0"))*volume;
        }

        void integrate(double const& dt) override {
            // index correspondence between vertices, heights, and height velocities
            for (auto it=height.begin(); it != height.end();) {
                if (!inMap(*vertices, it->first)) {                     // vertex index not present any more
                    it = height.erase(it);
                }
                else if ((vertices->at(it->first)).getType() != type) { // not a vertex any more
                    it = height.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=heightVelocity.begin();
                it != heightVelocity.end();) {
                if (!inMap(*vertices, it->first)) {                     // vertex index not present any more
                    it = heightVelocity.erase(it);
                }
                else if ((vertices->at(it->first)).getType() != type) { // not a vertex any more
                    it = heightVelocity.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(heightVelocity, it->first)                 // new vertex index
                    && (it->second).getType() == type) {
                    height[it->first] = parameters.at("h0");
                    heightVelocity[it->first] = 0;
                }
            }
            // integrate heights
            for (auto it=height.begin(); it != height.end(); ++it) {
                it->second += dt*heightVelocity[it->first];
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

class EdgePullForce : public VertexForce<ForcesType> {
/*
Pulling on Vertices at the ounter border with a constant radial tension.
*/
    protected:

        Mesh* const mesh;

    public:

        EdgePullForce(
            double const& Fpull_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("", // type is "" therefore addForce is called for every vertex
                                        // (meaning that VertexForce::addAllForces filters nothing)
                                        // we select the special boundary vertices in addForce
                {{"Fpull", Fpull_}},
                forces_, vertices_),
            mesh(mesh_) {}

        void addForce(Vertex const& vertex) override {

            if (vertex.getBoundary()) { // force is computed for neighbours of boundary vertices
//                 std::cout << "Found boundary vertex " << vertex.getIndex() << std::endl;

                // compute centre of mass of vertices
                std::vector<double> posCM = {0, 0};             // position of centre of mass
                long int N = 0;                                 // number of vertices
                for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                    if ((it->second).getType() == "vertex") {   // we consider only "vertex" vertices
                        std::vector<double> const pos =
                            (it->second).getPosition();
                        for (int dim=0; dim < 2; dim++) {
                            posCM[dim] += pos.at(dim);
                        }
                        N++;
                    }
                }
                for (int dim=0; dim < 2; dim++) { posCM[dim] /= N; }

                // compute force
                std::vector <long int> bindices =
                    mesh->getNeighbourVertices(vertex.getIndex())[0];
                for (long int idx: bindices) {
//                     std::cout << "Adding boundary force to " << idx << std::endl;

                    std::vector<double> const dir = mesh->wrapDiff( // vector going...
                        posCM,                                      // ... from centre of mass...
                        (vertices->at(idx)).getPosition(),          // ... to the edge vertex...
                        true);                                      // ... and normalised
//                     std::cout << dir[0] << " " << dir[1] << std::endl;
                    for (int dim=0; dim < 2; dim++) {
                        (*forces)[idx][dim] += parameters.at("Fpull")*dir[dim];
                    }
                }

            }
        }

        pybind11::tuple pybind11_getstate() const override;
};

class BoundaryTension : public VertexForce<ForcesType> {
/*
Line tension on open boundaries.
*/

    protected:

        Mesh* const mesh;

    public:

        BoundaryTension(
            double const& gamma_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("",
                {{"gamma", gamma_}},
                forces_, vertices_),
            mesh(mesh_) {}

        void addForce(Vertex const& vertex) override {

            if (vertex.getBoundary()) {
                std::vector<long int> neighbours =                      // indices of vertices on the boundary
                    mesh->getNeighbourVertices(vertex.getIndex())[0];
                long int const numberNeighbours = neighbours.size();    // number of vertices on the boundary
                std::vector<double> toPreviousNeighbour, toNextNeighbour;
                for (int i=0; i < numberNeighbours; i++) {
                    assert(vertices->at(neighbours[i]).getType() != "centre");
                    toPreviousNeighbour = mesh->wrapTo(
                        neighbours[i],
                        neighbours[pmod(i - 1, numberNeighbours)],
                        true);
                    toNextNeighbour = mesh->wrapTo(
                        neighbours[i],
                        neighbours[pmod(i + 1, numberNeighbours)],
                        true);
                    for (int dim=0; dim < 2; dim++) {
                        (*forces)[neighbours[i]][dim] +=
                            parameters.at("gamma")*(
                                toPreviousNeighbour[dim]
                                + toNextNeighbour[dim]);
                    }
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;
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
            Random* random_,
            ForcesType* forces_, VerticesType* const vertices_) :
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

        pybind11::tuple pybind11_getstate() const override;

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
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
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

        void deleteEdge(TopoChangeEdgeInfoType const& del) override {
            for (long int halfEdgeIndex : std::get<1>(del)) {
                if (inMap(tension, halfEdgeIndex)) {
                    tension.erase(halfEdgeIndex);
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

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
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_}, {"P0", P0_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), vertices(vertices_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getPerimeter() const
            { return perimeter; }
        void setPerimeter(std::map<long int, double> const& perimeter_)
            { perimeter = perimeter_; }

        std::map<long int, double> const& getNoise() const
            { return noise; }
        void setNoise(std::map<long int, double> const& noise_)
            { noise = noise_; }

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
            // noise integration
            double const dt_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dt_);
            for (auto it=noise.begin(); it != noise.end(); ++it) {
                // noise
                noise[it->first] =
                    (1 - dt_)*noise[it->first] + amp*random->gauss();
            }
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
            // tension setting
            for (auto it=noise.begin(); it != noise.end(); ++it) {
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

        void deleteEdge(TopoChangeEdgeInfoType const& del) override {
            for (long int halfEdgeIndex : std::get<1>(del)) {
                if (inMap(noise, halfEdgeIndex)) {
                    noise.erase(halfEdgeIndex);
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

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
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_}, {"P0", P0_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), vertices(vertices_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getPerimeter() const
            { return perimeter; }
        void setPerimeter(std::map<long int, double> const& perimeter_)
            { perimeter = perimeter_; }

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
            // tension integration
            double const dt_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dt_);
            for (auto it=tension.begin(); it != tension.end();) {
                // tension
                long int const cellIndices[2] =                                 // indices of cell of the junction
                    {halfEdges->at(
                        halfEdges->at(it->first)
                            .getNextIndex()).getToIndex(),
                    halfEdges->at(
                        halfEdges->at(halfEdges->at(it->first).getPairIndex())
                            .getNextIndex()).getToIndex()};
                if (!inMap(perimeter, cellIndices[0])
                    && !inMap(perimeter, cellIndices[1])) {                     // this junction does not belong to any cell
                    it = tension.erase(it);
                }
                else {
                    tension[it->first] = (1 - dt_)*tension[it->first]           // zero-mean deterministic part
                        + amp*random->gauss();                                  // stochastic part
                    for (long int cellIndex : cellIndices) {
                        if (inMap(perimeter, cellIndex)) {
                            tension[it->first] += dt_*parameters.at("Gamma")*(  // target deterministic part
                                perimeter[cellIndex] - parameters.at("P0"));
                        }
                    }
                    ++it;
                }
            }
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
        }

        void deleteEdge(TopoChangeEdgeInfoType const& del) override {
            for (long int halfEdgeIndex : std::get<1>(del)) {
                if (inMap(tension, halfEdgeIndex)) {
                    tension.erase(halfEdgeIndex);
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

class Model2 : public HalfEdgeForce<ForcesType> {

    protected:

        Mesh* const mesh;                       // mesh object
        Random* random;                         // random number generator
        std::map<long int, double> length;      // junction lengths
        std::map<long int, double> restLength;  // junction rest lengths
        std::map<long int, double> noise;       // noise in tension
        std::map<long int, double> tension;     // junction tension

    public:

        Model2(
            double const& Gamma_, double const& taur_,
            double const& sigma_, double const& taup_,
            Mesh* const mesh_, Random* random_,
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_}, {"taur", taur_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getLength() const
            { return length; }
        void setLength(std::map<long int, double> const& length_)
            { length = length_; }

        std::map<long int, double> const& getRestLength() const
            { return restLength; }
        void setRestLength(std::map<long int, double> const& restLength_)
            { restLength = restLength_; }

        std::map<long int, double> const& getNoise() const
            { return noise; }
        void setNoise(std::map<long int, double> const& noise_)
            { noise = noise_; }

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
            // rest length and noise integration
            double const dtr_ = dt/parameters.at("taur");
            double const dtp_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dtp_);
            for (auto it=length.begin(); it != length.end(); ++it) {
                // rest length
                restLength[it->first] =
                    (1 - dtr_)*restLength[it->first] + dtr_*length[it->first];
                // noise
                noise[it->first] =
                    (1 - dtp_)*noise[it->first] + amp*random->gauss();
            }
            // index correspondence between half-edges and junction lengths
            length.clear();
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if ((it->second).getType() == type) {                       // compute junction length
                    length[it->first] = mesh->getEdgeLength(it->first);
                }
            }
            // index correspondence between half-edges, rest lengths, and noises
            for (auto it=restLength.begin(); it != restLength.end();) {
                if (!inMap(length, it->first)) {                            // half-edge index not present any more or is not a junction
                    it = restLength.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=noise.begin(); it != noise.end();) {
                if (!inMap(length, it->first)) {                            // half-edge index not present any more or is not a junction
                    it = noise.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=length.begin(); it != length.end(); ++it) {        // new half-edge index
                if (!inMap(restLength, it->first)) {
                    restLength[it->first] = length[it->first];
                }
                if (!inMap(restLength, it->first)) {
                    noise[it->first] = parameters.at("sigma")*random->gauss();
                }
            }
            // clear tension
            tension.clear();
            // tension setting
            for (auto it=length.begin(); it != length.end(); ++it) {
                tension[it->first] =
                    parameters.at("Gamma")*(                                // deterministic part
                        length[it->first] - restLength[it->first])
                    + noise[it->first];                                     // noise part
            }
        }

        void deleteEdge(TopoChangeEdgeInfoType const& del) override {
            for (long int halfEdgeIndex : std::get<1>(del)) {
                if (inMap(restLength, halfEdgeIndex)) {
                    restLength.erase(halfEdgeIndex);
                }
                if (inMap(noise, halfEdgeIndex)) {
                    noise.erase(halfEdgeIndex);
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

class Model3 : public HalfEdgeForce<ForcesType> {

    protected:

        Mesh* const mesh;                       // mesh object
        Random* random;                         // random number generator
        std::map<long int, double> length;      // junction lengths
        std::map<long int, double> restLength;  // junction rest lengths
        std::map<long int, double> tension;     // junction tension

    public:

        Model3(
            double const& Gamma_,
            double const& sigma_, double const& taup_,
            Mesh* const mesh_, Random* random_,
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getLength() const
            { return length; }
        void setLength(std::map<long int, double> const& length_)
            { length = length_; }

        std::map<long int, double> const& getRestLength() const
            { return restLength; }
        void setRestLength(std::map<long int, double> const& restLength_)
            { restLength = restLength_; }

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
            // rest length integration
            double const dtp_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dtp_);
            for (auto it=restLength.begin(); it != restLength.end(); ++it) {
                restLength[it->first] =
                    (1 - dtp_)*restLength[it->first] + dtp_*length[it->first]   // deterministic part
                    + amp*parameters.at("sigma")*random->gauss();               // noise part
            }
            // index correspondence between half-edges and junction lengths
            length.clear();
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if ((it->second).getType() == type) {                           // compute junction length
                    length[it->first] = mesh->getEdgeLength(it->first);
                }
            }
            // index correspondence between half-edges and rest lengths
            for (auto it=restLength.begin(); it != restLength.end();) {
                if (!inMap(length, it->first)) {                                // half-edge index not present any more or is not a junction
                    it = restLength.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=length.begin(); it != length.end(); ++it) {            // new half-edge index
                if (!inMap(restLength, it->first)) {
                    restLength[it->first] = length[it->first]
                        + parameters.at("sigma")*random->gauss();
                }
            }
            // clear tension
            tension.clear();
            // tension setting
            for (auto it=length.begin(); it != length.end(); ++it) {
                tension[it->first] =
                    parameters.at("Gamma")*(
                        length[it->first] - restLength[it->first]);
            }
        }

        void deleteEdge(TopoChangeEdgeInfoType const& del) override {
            for (long int halfEdgeIndex : std::get<1>(del)) {
                if (inMap(restLength, halfEdgeIndex)) {
                    restLength.erase(halfEdgeIndex);
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

class Model4 : public HalfEdgeForce<ForcesType> {

    protected:

        Mesh* const mesh;                       // mesh object
        Random* random;                         // random number generator
        std::map<long int, double> length;      // junction lengths
        std::map<long int, double> restLength;  // junction rest lengths
        std::map<long int, double> tension;     // junction tension

    public:

        Model4(
            double const& Gamma_, double const& taur_,
            double const& sigma_, double const& taup_,
            Mesh* const mesh_, Random* random_,
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Gamma", Gamma_}, {"taur", taur_},
                    {"sigma", sigma_}, {"taup", taup_}},
                forces_, halfEdges_),
            mesh(mesh_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getLength() const
            { return length; }
        void setLength(std::map<long int, double> const& length_)
            { length = length_; }

        std::map<long int, double> const& getRestLength() const
            { return restLength; }
        void setRestLength(std::map<long int, double> const& restLength_)
            { restLength = restLength_; }

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
            // tension and rest length integration
            double const dtp_ = dt/parameters.at("taup");
            double const amp = parameters.at("sigma")*sqrt(2.*dtp_);
            double const dtr_ = dt/parameters.at("taur");
            for (auto it=tension.begin(); it != tension.end(); ++it) {
                // tension
                tension[it->first] =
                    (1 - dtp_)*tension[it->first]                               // zero-mean deterministic part
                    + dtp_*parameters.at("Gamma")*(                             // target deterministic part
                        length[it->first] - restLength[it->first])
                    + amp*parameters.at("sigma")*random->gauss();               // noise part
                // rest length
                restLength[it->first] =
                    (1 - dtr_)*restLength[it->first] + dtr_*length[it->first];
            }
            // index correspondence between half-edges and junction lengths
            length.clear();
            for (auto it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if ((it->second).getType() == type) {                           // compute junction length
                    length[it->first] = mesh->getEdgeLength(it->first);
                }
            }
            // index correspondence between half-edges and rest lengths
            for (auto it=restLength.begin(); it != restLength.end();) {
                if (!inMap(length, it->first)) {                                // half-edge index not present any more or is not a junction
                    it = restLength.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=length.begin(); it != length.end(); ++it) {            // new half-edge index
                if (!inMap(restLength, it->first)) {
                    restLength[it->first] = length[it->first]
                        + parameters.at("sigma")*random->gauss();
                }
            }
            // clear tension
            tension.clear();
            // tension setting
            for (auto it=length.begin(); it != length.end(); ++it) {
                tension[it->first] =
                    parameters.at("Gamma")*(
                        length[it->first] - restLength[it->first]);
            }
        }

        void deleteEdge(TopoChangeEdgeInfoType const& del) override {
            for (long int halfEdgeIndex : std::get<1>(del)) {
                if (inMap(restLength, halfEdgeIndex)) {
                    restLength.erase(halfEdgeIndex);
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

/*
 *  KERATIN
 *
 */

class KeratinModel : public VertexForce<ForcesType> {
/*
Keratin model.
*/

    protected:

        Mesh* const mesh;                       // mesh object
        Random* random;                         // random number generator

        double* const time;                     // integrated time
        double const time0;                     // time at which force was initialised

        std::map<long int, double> keratin;     // cell keratin concentration
        std::map<long int, double> targetArea;  // target area

        // degrees of freedom computed in KeratinModel::addAllForces
        std::map<long int, double> pressure;    // cell pressure
        std::map<long int, double> area;        // cell area
        std::map<long int, double> tension;     // bonds keratin tension

    public:

        KeratinModel(
            double const& K_, double const& taur_,
            double const& Gamma_, double const& p0_,
            double const& l0_, double const& alpha_, double const& kth_,
            double const& tau_, double const& sigma_,
            double const& ron_, double const& k0_, double const& pr0_,
            Mesh* const mesh_, Random* random_,
            double* const time_, double const& time0_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"K", K_}, {"taur", taur_},
                {"Gamma", Gamma_}, {"p0", p0_},
                {"l0", l0_}, {"alpha", alpha_}, {"kth", kth_},
                {"tau", tau_}, {"sigma", sigma_},
                {"ron", ron_}, {"k0", k0_}, {"pr0", pr0_}},
                forces_, vertices_),
            mesh(mesh_), random(random_), time(time_), time0(time0_) {

            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == "centre") {   // loop over all cell centres
                    keratin.emplace(it->first,              // zero keratin initial value
                        0);
                    area.emplace(it->first,                 // initial area
                        mesh->getVertexToNeighboursArea(it->first));
                    targetArea.emplace(it->first,           // initial area as initial target area
                        area[it->first]);
                }
            }
        }

        double const& getTime0() const
            { return time0; }

        std::map<long int, double> const& getKeratin() const
            { return keratin; }
        void setKeratin(std::map<long int, double> const& keratin_)
            { keratin = keratin_; }

        std::map<long int, double> const& getTargetArea() const
            { return targetArea; }
        void setTargetArea(std::map<long int, double> const& targetArea_)
            { targetArea = targetArea_; }

        std::map<long int, double> const& getPressure() const
            { return pressure; }

        std::map<long int, double> const& getArea() const
            { return area; }

        std::map<long int, double> const& getTension() const
            { return tension; }

        void addForce(Vertex const& vertex) override {

            // reset cell pressure
            pressure[vertex.getIndex()] = 0;

            // cell area and perimeter
            area.emplace(vertex.getIndex(),
                mesh->getVertexToNeighboursArea(vertex.getIndex()));
            double const A0 = targetArea[vertex.getIndex()];
            double const perimeter =
                mesh->getVertexToNeighboursPerimeter(vertex.getIndex());
            double const P0 = parameters.at("p0")*sqrt(A0);

            std::vector<long int> neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            long int const numberNeighbours = neighbourVerticesIndices.size();
            std::vector<double> radiusToCorner;
            double distRadiusToCorner;
            std::vector<double> toPreviousNeighbour, toNextNeighbour;
            std::vector<double> crossToPreviousNeighbour, crossToNextNeighbour;
            double distToPreviousNeighbour, distToNextNeighbour;
            long int halfEdgeIndex;
            double force;
            for (int i=0; i < numberNeighbours; i++) {  // loop over neighbours
                assert(vertices->at(neighbourVerticesIndices[i]).getType()
                    != "centre");
                // from cell centre to cell corner
                radiusToCorner = mesh->wrapTo(
                    vertex.getIndex(),
                    neighbourVerticesIndices[i]);
                distRadiusToCorner = sqrt(
                    radiusToCorner[0]*radiusToCorner[0]
                    + radiusToCorner[1]*radiusToCorner[1]);
                // towards previous neighbour
                toPreviousNeighbour = mesh->wrapTo(
                    neighbourVerticesIndices[i],
                    neighbourVerticesIndices[pmod(i - 1, numberNeighbours)],
                    false);
                crossToPreviousNeighbour = cross2z(toPreviousNeighbour);
                distToPreviousNeighbour = sqrt(
                    toPreviousNeighbour[0]*toPreviousNeighbour[0]
                    + toPreviousNeighbour[1]*toPreviousNeighbour[1]);
                // towards next neighbour
                toNextNeighbour = mesh->wrapTo(
                    neighbourVerticesIndices[i],
                    neighbourVerticesIndices[pmod(i + 1, numberNeighbours)],
                    false);
                crossToNextNeighbour = cross2z(toNextNeighbour);
                distToNextNeighbour = sqrt(
                    toNextNeighbour[0]*toNextNeighbour[0]
                    + toNextNeighbour[1]*toNextNeighbour[1]);
                for (int dim=0; dim < 2; dim++) {
                    // area force
                    force = (parameters.at("K")/2.)*(
                        area[vertex.getIndex()] - A0)*( // cell area
                            crossToPreviousNeighbour[dim]
                                - crossToNextNeighbour[dim]);
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                    pressure[vertex.getIndex()] -=      // pressure = -Tr(stress)
                        radiusToCorner[dim]             // radius towards vertex
                            *(-force);                  // force applied on cell centre
                    // perimeter force
                    force = parameters.at("Gamma")*(
                        perimeter - P0)*(               // cell perimeter
                            toPreviousNeighbour[dim]        // force from previous neighbour
                                /distToPreviousNeighbour);
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                    force = parameters.at("Gamma")*(
                        perimeter - P0)*(               // cell perimeter
                            toNextNeighbour[dim]            // force from next neighbour
                                /distToNextNeighbour);
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                    pressure[vertex.getIndex()] -=      // pressure = -Tr(stress)
                        toNextNeighbour[dim]            // radius from cell corner to next neighbour
                            *(force);                   // force applied on cell corner
                    // keratin force between cell centre and cell corner
                    halfEdgeIndex = mesh->getHalfEdgeIndex(
                        vertex.getIndex(),
                        neighbourVerticesIndices[i]);
                    if (dim == 0) { // tension is defined once per half-edge
                        tension.emplace(halfEdgeIndex,
                            2*parameters.at("alpha")*std::max(0.,
                                keratin[vertex.getIndex()]  // cell keratin concentration
                                    - parameters.at("kth"))*(
                                        distRadiusToCorner  // length between vertices
                                            - parameters.at("l0")));
                    }
                    force = -tension[halfEdgeIndex]                 // tension
                        *radiusToCorner[dim]/distRadiusToCorner;    // normalised vector;
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                    (*forces)[vertex.getIndex()][dim] -= force;             // force on cell centre
                    pressure[vertex.getIndex()] -=      // pressure = -Tr(stress)
                        radiusToCorner[dim]             // radius from cell centre to cell corner
                            *(-force);                  // force applied on cell centre
                    // keratin force between cell corner and next neighbour
                    halfEdgeIndex = mesh->getHalfEdgeIndex(
                        neighbourVerticesIndices[i],
                        neighbourVerticesIndices[
                            pmod(i + 1, numberNeighbours)]);
                    if (dim == 0) { // tension is defined once per half-edge
                        tension.emplace(halfEdgeIndex,
                            parameters.at("alpha")*std::max(0.,
                                keratin[vertex.getIndex()]  // cell keratin concentration
                                    - parameters.at("kth"))*(
                                        distToNextNeighbour // length between vertices
                                            - parameters.at("l0")));
                    }
                    force = tension[halfEdgeIndex]                  // tension
                        *toNextNeighbour[dim]/distToNextNeighbour;  // normalised vector
                    (*forces)[neighbourVerticesIndices[                     // force on vertex i
                        i]][dim] += force;
                    (*forces)[neighbourVerticesIndices[                     // force on next neighbour
                        pmod(i + 1, numberNeighbours)]][dim] -= force;
                    pressure[vertex.getIndex()] -=      // pressure = -Tr(stress)
                        toNextNeighbour[dim]            // radius from cell corner to next neighbour
                            *(force);                   // force applied on cell centre
                }
            }

            pressure[vertex.getIndex()] /= area[vertex.getIndex()];
        }

        void addAllForces() override {
            area.clear(); tension.clear();
            VertexForce<ForcesType>::addAllForces();
        }

        void integrate(double const& dt) override {

            // integrate keratin concentration
            double const dt_ = dt/parameters.at("tau");
            double const amp = parameters.at("sigma")*sqrt(2.*dt_);
            for (auto it=keratin.begin(); it != keratin.end(); ++it) {
                double const kon = 1 + (*time - time0)*parameters.at("ron");    // time-dependent on rate
                double const koff = 1 + exp(parameters.at("k0")*(               // keratin-dependent off rate
                    pressure[it->first] - parameters.at("pr0")));               // addAllForces sets pressure
                keratin[it->first] +=
                    dt_*(kon - koff*keratin[it->first]) // deterministic part
                    + amp*random->gauss();              // stochastic part
            }
            // integrate target area
            double const dtr_ = dt/parameters.at("taur");
            for (auto it=targetArea.begin(); it != targetArea.end(); ++it) {
                targetArea[it->first] +=
                    -dtr_*(targetArea[it->first] - area[it->first]);
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

#endif

