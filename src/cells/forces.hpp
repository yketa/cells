/*
Forces definitions.
*/

#ifndef FORCES_HPP
#define FORCES_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>

#include "assert.hpp"
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

class SurfaceForce : public VertexForce<ForcesType> {
/*
Cell surface minimising force at constant volume.
*/

    protected:

        Mesh* const mesh;

        std::map<long int, double> volume;

    public:

        SurfaceForce(
            double const& Lambda_, double const& V0_, double const& tauV_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"Lambda", Lambda_}, {"V0", V0_}, {"tauV", tauV_}},
                forces_, vertices_),
            mesh(mesh_) {
            // initialise volumes
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == type) {
                    volume.emplace(it->first, parameters.at("V0"));
                }
            }
        }

        std::map<long int, double> const& getVolume() const
            { return volume; }
        void setVolume(std::map<long int, double> const& volume_)
            { volume = volume_; }

        std::map<long int, double> getHeight() const {
            std::map<long int, double> height;
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == type) {
                    double const area =
                        mesh->getVertexToNeighboursArea(it->first);
                    height.emplace(it->first, volume.at(it->first)/area);
                }
            }
            return height;
        }

        void addForce(Vertex const& vertex) override {
            double const area = mesh->getVertexToNeighboursArea(
                vertex.getIndex());
            double const perimeter = mesh->getVertexToNeighboursPerimeter(
                vertex.getIndex());
            double const V0 = volume.at(vertex.getIndex());
            std::vector<long int> neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            long int const numberNeighbours = neighbourVerticesIndices.size();
            std::vector<double> crossToPreviousNeighbour, crossToNextNeighbour;
            std::vector<double> toPreviousNeighbour, toNextNeighbour;
            for (int i=0; i < numberNeighbours; i++) {
                assert(vertices->at(neighbourVerticesIndices[i]).getType()
                    != "centre");
                // force analogous to area force
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
                        (parameters.at("Lambda")/2)*(
                            2 - V0*perimeter/(area*area))*(
                                crossToPreviousNeighbour[dim]
                                    - crossToNextNeighbour[dim]);
                }
                // force analogous to perimeter force
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
                        parameters.at("Lambda")*(V0/area)*(
                            toPreviousNeighbour[dim] + toNextNeighbour[dim]);
                }
            }
        }

        void integrate(double const& dt) override {

            // index correspondence between vertices and volumes
            for (auto it=volume.begin(); it != volume.end();) {
                if (!inMap(*vertices, it->first)) {                     // vertex index not present any more
                    it = volume.erase(it);
                }
                else if ((vertices->at(it->first)).getType() != type) { // not a vertex any more
                    it = volume.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(volume, it->first)
                    && (it->second).getType() == type) {                // new vertex index
                    volume[it->first] = parameters.at("V0");
                }
            }

            if (parameters.at("tauV") == 0) return;

            // integrate volume
            for (auto it=volume.begin(); it != volume.end(); ++it) {
//                 // exponential relaxation
//                 it->second += (parameters.at("V0") - it->second)
//                     *dt/parameters.at("tauV");
                // linear variation
                it->second += (parameters.at("V0")/parameters.at("tauV"))*dt;
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
            double const& kV_, double const& H0_, double const& A0_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"kV", kV_}, {"H0", H0_}, {"A0", A0_}},
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
                            volume - parameters.at("H0")*parameters.at("A0"))*(
                            crossToPreviousNeighbour[dim]
                                - crossToNextNeighbour[dim]);
                }
            }
            heightVelocity[vertex.getIndex()] += -parameters.at("kV")*
                (volume - parameters.at("H0")*parameters.at("A0"))*volume;
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
                if (!inMap(heightVelocity, it->first)                   // new vertex index
                    && (it->second).getType() == type) {
                    height[it->first] = parameters.at("H0");
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

class LinearVolumeForce : public VertexForce<ForcesType> {
/*
Cell volume restoring force which is linear in edge length.
*/

    protected:

        Mesh* const mesh;                           // mesh object
        std::map<long int, double> height;          // height of each cell

        std::map<long int, double> heightVelocity;  // height velocity of each cell

    public:

        LinearVolumeForce(
            double const& kA_, double const& A0_,
            double const& kP_, double const& P0_,
            double const& taur_, double const H0_, double const& taua_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"kA", kA_}, {"A0", A0_},
                {"kP", kP_}, {"P0", P0_},
                {"taur", taur_}, {"H0", H0_}, {"taua", taua_}},
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
            // cell dimensions
            double const perimeter =
                mesh->getVertexToNeighboursPerimeter(vertex.getIndex());
            double const area =
                mesh->getVertexToNeighboursArea(vertex.getIndex());
            double const volume =
                height[vertex.getIndex()]*area;
            // neighbours
            HalfEdgesType const& halfEdges = mesh->getHalfEdges();
            std::vector<std::vector<long int>> const neighbours =
                mesh->getNeighbourVertices(vertex.getIndex());
            int const numberNeighbours = neighbours.at(0).size();
            // forces on vertices
            for (int i=0; i < numberNeighbours; i++) {
                assert(vertices->at(neighbours.at(0).at(i)).getType()
                    != "centre");
                std::vector<double> const toPreviousVertex = mesh->wrapTo(
                    neighbours.at(0).at(i),
                    neighbours.at(0).at(pmod(i - 1, numberNeighbours)),
                    true);
                std::vector<double> const toNextVertex = mesh->wrapTo(
                    neighbours.at(0).at(i),
                    neighbours.at(0).at(pmod(i + 1, numberNeighbours)),
                    true);
                std::vector<double> const toInterior =
                    cross2z(mesh->wrapTo(
                        neighbours.at(0).at(pmod(i + 1, numberNeighbours)),
                        neighbours.at(0).at(pmod(i - 1, numberNeighbours)),
                        true));
                for (int dim=0; dim < 2; dim++) {
                    // perimeter force
                    (*forces)[neighbours.at(0).at(i)][dim] +=   // perimeter conservation
                        parameters.at("kP")
                            *(perimeter - parameters.at("P0"))
                            *(toPreviousVertex[dim] + toNextVertex[dim]);
                    // volume force
                    (*forces)[neighbours.at(0).at(i)][dim] +=   // area conservation
                        parameters.at("kA")*parameters.at("A0")
                            *(sqrt(area) - sqrt(parameters.at("A0")))
                            *toInterior[dim];
                }
            }
            // forces on heights
            heightVelocity[vertex.getIndex()] +=                // volume conservation
                -(1./parameters.at("taur"))*pow(parameters.at("A0"), -1./4.)
                    *(sqrt(volume)
                        - sqrt(parameters.at("H0")*parameters.at("A0")));
            for (long int halfEdgeIndex : neighbours.at(1)) {
                // neighbour cell index
                long int const neighbourVertexIndex =
                    halfEdges.at(
                        halfEdges.at(
                            halfEdges.at(
                                halfEdges.at(halfEdgeIndex).getNextIndex()
                            ).getPairIndex()
                        ).getNextIndex()
                    ).getToIndex();
                if (vertices->at(neighbourVertexIndex).getBoundary())   // ignore boundary vertices
                    { continue; }
                assert(                                                 // the vertex should be a cell centre
                    vertices->at(neighbourVertexIndex).getType()
                        == "centre");
                // force
                heightVelocity[vertex.getIndex()] +=            // penalise height differences
                    (1./parameters.at("taua"))
                        *(height[neighbourVertexIndex]
                            - height[vertex.getIndex()]);
            }
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
                if (!inMap(heightVelocity, it->first)                   // new vertex index
                    && (it->second).getType() == type) {
                    height[it->first] =
                        parameters.at("H0")*parameters.at("A0")
                            /mesh->getVertexToNeighboursArea(it->first);
                    heightVelocity[it->first] = 0;
                }
            }
            // integrate heights
            for (auto it=height.begin(); it != height.end(); ++it) {
                it->second += dt*heightVelocity[it->first];
                it->second = std::max(it->second, 0.);
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

class PressureForce : public VertexForce<ForcesType> {
/*
Force acting on boundary vertices outward from their barycentre. Forces on each
boundary vertex are either in the direction of the line going from the
barycentre to the vertex with a norm equal to the force scale (fixedForce=true)
or are derived from a constant inner pressure equal to the force scale times
the number of boundary vertices divided by the perimeter of the boundary
(fixedForce=false).
*/
    protected:

        Mesh* const mesh;

    public:

        PressureForce(
            double const& F_, bool const& fixedForce_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("", // type is "" therefore addForce is called for every vertex
                                        // (meaning that VertexForce::addAllForces filters nothing)
                                        // we select the special boundary vertices in addForce
                {{"F", F_}, {"fixedForce", fixedForce_}},
                forces_, vertices_),
            mesh(mesh_) {}

        void addForce(Vertex const& vertex) override {

            if (vertex.getBoundary()) { // force is computed for neighbours of boundary vertices

                // compute centre of mass of vertices
                mesh->moveToNeighboursBarycentre(vertex.getIndex());    // move to barycentre to compute radial vectors
                std::vector<double> const posCM =                       // position of centre of mass (= barycentre)
                    (vertices->at(vertex.getIndex())).getPosition();

                // compute neighbours
                std::vector<long int> const neighbours =    // boundary vertices are neighbours of the boundary vertex
                    mesh->getNeighbourVertices(vertex.getIndex())[0];
                int const numberNeighbours = neighbours.size();
                std::map<long int, std::vector<double>> crossFromPrevious;
                std::map<long int, std::vector<double>> crossToNext;
                if (!parameters.at("fixedForce")) {
                    for (int i=0; i < numberNeighbours; i++) {
                        assert(vertices->at(neighbours[i]).getType()
                            != "centre");
                        assert(
                            cross2(
                                // CAREFUL: neighbours of boundary are in reverse order
                                mesh->wrapTo(
                                    vertex.getIndex(),
                                    neighbours[pmod(i + 1, numberNeighbours)]),
                                mesh->wrapTo(
                                    vertex.getIndex(),
                                    neighbours[i]))
                            > 0);
                        crossFromPrevious.emplace(neighbours[i],
                            cross2z(
                                // CAREFUL: neighbours of boundary are in reverse order
                                mesh->wrapTo(
                                    neighbours[pmod(i + 1, numberNeighbours)],
                                    neighbours[i],
                                    false)));
                        crossToNext.emplace(neighbours[i],
                            cross2z(
                                // CAREFUL: neighbours of boundary are in reverse order
                                mesh->wrapTo(
                                    neighbours[i],
                                    neighbours[pmod(i - 1, numberNeighbours)],
                                    false)));
                    }
                }

                // compute force scale
                double const scale =                                        // scale used for force
                    [this, &crossFromPrevious, &crossToNext,
                        &numberNeighbours]() {
                        if (parameters.at("fixedForce")) {
                            return (this->parameters).at("F");              // dimension of a force
                        }
                        else {
                            double sumLength = 0;
                            for (auto it=crossFromPrevious.begin();
                                it != crossFromPrevious.end(); ++it) {
                                sumLength += sqrt(
                                    pow(
                                        crossFromPrevious[it->first][0]
                                            + crossToNext[it->first][0],
                                        2)
                                    + pow(
                                        crossFromPrevious[it->first][1]
                                            + crossToNext[it->first][1],
                                        2));
                            }
                            return numberNeighbours                         // ensures average norm if F
                                *(this->parameters).at("F")/sumLength;      // dimension of a force per length
                        }
                    }();

                // set force
                for (long int neighbourIndex : neighbours) {
                    if (parameters.at("fixedForce")) {
                        std::vector<double> const dir = mesh->wrapDiff(     // vector going...
                            posCM,                                          // ... from centre of mass...
                            (vertices->at(neighbourIndex)).getPosition(),    // ... to the edge vertex...
                            true);                                          // ... and normalised
                        for (int dim=0; dim < 2; dim++) {
                            (*forces)[neighbourIndex][dim] += scale
                                *dir[dim];
                        }
                    }
                    else {
                        for (int dim=0; dim < 2; dim++) {
                            (*forces)[neighbourIndex][dim] += scale
                                *(crossFromPrevious[neighbourIndex][dim]
                                    + crossToNext[neighbourIndex][dim]);
                        }
                    }
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;
};

class TensionForce : public HalfEdgeForce<ForcesType> {
/*
Constant tension on all junctions.
*/

    protected:

        Mesh* const mesh;

    public:

        TensionForce(
            double const& Lambda_,
            Mesh* const mesh_,
            ForcesType* forces_, HalfEdgesType* const halfEdges_) :
            HalfEdgeForce<ForcesType>("junction",   // acts on half-edges of type "junction"
                {{"Lambda", Lambda_}},
                forces_, halfEdges_),
            mesh(mesh_) {}

        void addForce(HalfEdge const& halfEdge) override {

            std::vector<double> const junctionDir = // unit vector in the direction of the junction
                mesh->getHalfEdgeVector(halfEdge.getIndex(), true);

            for (int dim=0; dim < 2; dim++) {
                (*forces)[halfEdge.getFromIndex()][dim] +=
                    parameters.at("Lambda")*junctionDir[dim];
                (*forces)[halfEdge.getToIndex()][dim] +=
                    -parameters.at("Lambda")*junctionDir[dim];
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

class GrowingAreaPerimeterForce : public VertexForce<ForcesType> {
/*
Cell area and perimeter restoring force with growing target area.
*/

    protected:

        Mesh* const mesh;

        std::map<long int, double> targetArea;  // target area

    public:

        GrowingAreaPerimeterForce(
            double const& kA_, double const& s0_,
                double const& A0_, double const& tauA_,
            Mesh* const mesh_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"kA", kA_}, {"s0", s0_},
                    {"A0", A0_}, {"tauA", tauA_}},
                forces_, vertices_),
            mesh(mesh_) {

            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == "centre") {
                    targetArea.emplace(it->first,   // initial area as initial target area
                        mesh->getVertexToNeighboursArea(it->first));
                }
            }
        }

        std::map<long int, double> const& getTargetArea() const
            { return targetArea; }
        void setTargetArea(std::map<long int, double> const& targetArea_)
            { targetArea = targetArea_; }

        void addForce(Vertex const& vertex) override {

            double const area = mesh->getVertexToNeighboursArea(
                vertex.getIndex());
            double const perimeter = mesh->getVertexToNeighboursPerimeter(
                vertex.getIndex());

            double const A0 = targetArea.at(vertex.getIndex());
            double const P0 = parameters.at("s0")*sqrt(A0);

            std::vector<long int> neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            long int const nNeighbours = neighbourVerticesIndices.size();

            for (int i=0; i < nNeighbours; i++) {
                assert(vertices->at(neighbourVerticesIndices[i]).getType()
                    != "centre");
                // area force
                std::vector<double> const crossToPreviousNeighbour =
                    cross2z(mesh->wrapTo(
                        neighbourVerticesIndices[i],
                        neighbourVerticesIndices[pmod(i - 1, nNeighbours)],
                        false));
                std::vector<double> const crossToNextNeighbour =
                    cross2z(mesh->wrapTo(
                        neighbourVerticesIndices[i],
                        neighbourVerticesIndices[pmod(i + 1, nNeighbours)],
                        false));
                for (int dim=0; dim < 2; dim++) {
                    (*forces)[neighbourVerticesIndices[i]][dim] +=
                        (parameters.at("kA")/2.)*(
                            area - A0)*(
                            crossToPreviousNeighbour[dim]
                                - crossToNextNeighbour[dim]);
                }
                // perimeter force
                std::vector<double> const toPreviousNeighbour =
                    mesh->wrapTo(
                        neighbourVerticesIndices[i],
                        neighbourVerticesIndices[pmod(i - 1, nNeighbours)],
                        true);
                std::vector<double> const toNextNeighbour =
                    mesh->wrapTo(
                        neighbourVerticesIndices[i],
                        neighbourVerticesIndices[pmod(i + 1, nNeighbours)],
                        true);
                for (int dim=0; dim < 2; dim++) {
                    (*forces)[neighbourVerticesIndices[i]][dim] +=
                        parameters.at("kA")*A0*(
                            perimeter - P0)*(
                            toPreviousNeighbour[dim] + toNextNeighbour[dim]);
                }
            }
        }

        void integrate(double const& dt) override {

            // index correspondence between vertices and target areas
            for (auto it=targetArea.begin(); it != targetArea.end();) {
                if (!inMap(*vertices, it->first)) {                     // vertex index not present any more
                    it = targetArea.erase(it);
                }
                else if ((vertices->at(it->first)).getType() != type) { // not a vertex any more
                    it = targetArea.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(targetArea, it->first)
                    && (it->second).getType() == type) {                // new vertex index
                    targetArea[it->first] = mesh->getVertexToNeighboursArea(
                        it->first);
                }
            }

            // integrate target areas
            for (auto it=targetArea.begin(); it != targetArea.end(); ++it) {
                it->second += dt*parameters.at("A0")/parameters.at("tauA");
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
            if (parameters.at("taup") != 0) {
                double const amp = sqrt(2.*dt/parameters.at("taup"));
                for (auto it=theta.begin(); it != theta.end(); ++it) {
                    it->second += amp*random->gauss();
                }
            }
            else {
                for (auto it=theta.begin(); it != theta.end(); ++it) {
                    it->second += 2.*std::numbers::pi*random->random01();
                }
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

class ActiveBrownianCellForce : public VertexForce<ForcesType> {
/*
Active Brownian self-propulsion force acting on cell centres.
*/

    protected:

        Mesh* const mesh;                   // mesh object
        Random* random;                     // random number generator
        std::map<long int, double> theta;   // orientation of self-propulsion force

    public:

        ActiveBrownianCellForce(
            double const& v0_, double const& taup_,
            Mesh* const mesh_, Random* random_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("vertex",
                {{"v0", v0_}, {"taup", taup_}},
                forces_, vertices_),
            mesh(mesh_), random(random_)
            { integrate(0); }

        std::map<long int, double> const& getTheta() const
            { return theta; }
        void setTheta(std::map<long int, double> const& theta_)
            { theta = theta_; }

        void addForce(Vertex const& vertex) override {
            std::vector<long int> const neighbours =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            std::vector<double> neighbourThetas;
            for (long int neighbour : neighbours) {
                if (inMap(theta, neighbour)) {
                    neighbourThetas.push_back(theta.at(neighbour));
                }
            }
            long int const nNeighbours = neighbourThetas.size();
            if (nNeighbours > 0) {
                for (double neighbourTheta : neighbourThetas) {
                    for (int dim=0; dim < 2; dim++) {
                        (*forces)[vertex.getIndex()][dim] +=
                            parameters.at("v0")*cos(neighbourTheta
                                - dim*std::numbers::pi/2.)
                                /nNeighbours;
                    }
                }
            }
        }

        void integrate(double const& dt) override {
            // index correspondence between vertices and orientations
            for (auto it=theta.begin(); it != theta.end();) {
                if (!inMap(*vertices, it->first)) {                         // vertex index not present any more
                    it = theta.erase(it);
                }
                else if ((vertices->at(it->first)).getType() != "centre") { // not a cell centre any more
                    it = theta.erase(it);
                }
                else {
                    ++it;
                }
            }
            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if (!inMap(theta, it->first)                                // new vertex index
                    && (it->second).getType() == "centre") {
                    theta[it->first] = 2.*std::numbers::pi*random->random01();
                }
            }
            // integration
            if (parameters.at("taup") != 0) {
                double const amp = sqrt(2.*dt/parameters.at("taup"));
                for (auto it=theta.begin(); it != theta.end(); ++it) {
                    it->second += amp*random->gauss();
                }
            }
            else {
                for (auto it=theta.begin(); it != theta.end(); ++it) {
                    it->second += 2.*std::numbers::pi*random->random01();
                }
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
        HalfEdgesType const& halfEdges;         // half edges

        std::map<long int, double> keratin;     // cell keratin concentration
        std::map<long int, double> targetArea;  // target area

        // degrees of freedom computed in KeratinModel::addAllForces
        std::map<long int, double> area;        // cell area
        std::map<long int, double> pressure;    // cell pressure
        std::map<long int, double> tension;     // cell tension

        double keffect(long int const& index) const {
            // keratin effect on area elastic constant and target area relaxation time
            return 1 +
                parameters.at("beta")*std::max(0.,
                    keratin.at(index) - parameters.at("kth"));
        }

    public:

        KeratinModel(
            double const& K_, double const& A0, double const& taur_,
            double const& Gamma_, double const& p0_,
            double const& alpha_, double const& beta_, double const& kth_,
            double const& tau_, double const& sigma_, double const& ron_,
            Mesh* const mesh_, Random* random_,
            ForcesType* forces_, VerticesType* const vertices_) :
            VertexForce<ForcesType>("centre",
                {{"K", K_}, {"A0", A0}, {"taur", taur_},
                {"Gamma", Gamma_}, {"p0", p0_},
                {"alpha", alpha_}, {"beta", beta_}, {"kth", kth_},
                {"tau", tau_}, {"sigma", sigma_}, {"ron", ron_}},
                forces_, vertices_),
            mesh(mesh_), random(random_), halfEdges(mesh->getHalfEdges()) {

            for (auto it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == "centre") {   // loop over all cell centres
                    keratin.emplace(it->first,              // zero keratin initial value
                        0);
//                         std::max(
//                             0.,
//                             parameters.at("kth")*(1 + random->gauss())));
                    targetArea.emplace(it->first,           // initial area as initial target area
                        std::max(
                            parameters.at("A0"),
                            mesh->getVertexToNeighboursArea(it->first)));
                }
            }
        }

        std::map<long int, double> const& getKeratin() const
            { return keratin; }
        void setKeratin(std::map<long int, double> const& keratin_)
            { keratin = keratin_; }

        std::map<long int, double> const& getTargetArea() const
            { return targetArea; }
        void setTargetArea(std::map<long int, double> const& targetArea_)
            { targetArea = targetArea_; }

        std::map<long int, double> const& getArea() const
            { return area; }

        std::map<long int, double> const& getPressure() const
            { return pressure; }
        std::map<long int, double> getPressureJunction() const {
            std::map<long int, double> pressure_junction;
            std::map<long int, HalfEdge> const& halfEdges =
                mesh->getHalfEdges();
            std::vector<long int> const halfEdgeIndices =
                mesh->getHalfEdgeIndicesByType("junction");
            for (long int halfEdgeIndex : halfEdgeIndices) {
                pressure_junction.emplace(halfEdgeIndex, 0);
                for (long int index :
                    {halfEdgeIndex,
                        (halfEdges.at(halfEdgeIndex)).getPairIndex()}) {
                    long int const vertexIndex =
                        (halfEdges.at(
                            (halfEdges.at(index)
                                ).getNextIndex())
                            ).getToIndex();
                    if (inMap(pressure, vertexIndex)) {
                        pressure_junction[halfEdgeIndex] +=
                            -pressure.at(vertexIndex)/parameters.at("alpha");
                    }
                }
            }
            return pressure_junction;
        }

        std::map<long int, double> const& getTension() const
            { return tension; }
        std::map<long int, double> getTensionJunction() const {
            std::map<long int, double> tension_junction;
            std::map<long int, HalfEdge> const& halfEdges =
                mesh->getHalfEdges();
            std::vector<long int> const halfEdgeIndices =
                mesh->getHalfEdgeIndicesByType("junction");
            for (long int halfEdgeIndex : halfEdgeIndices) {
                tension_junction.emplace(halfEdgeIndex, 0);
                for (long int index :
                    {halfEdgeIndex,
                        (halfEdges.at(halfEdgeIndex)).getPairIndex()}) {
                    long int const vertexIndex =
                        (halfEdges.at(
                            (halfEdges.at(index)
                                ).getNextIndex())
                            ).getToIndex();
                    if (inMap(tension, vertexIndex)) {
                        tension_junction[halfEdgeIndex] +=
                            tension.at(vertexIndex);
                    }
                }
            }
            return tension_junction;
        }

        void addForce(Vertex const& vertex) override {

            // area elasticity
            double const Kk =                               // keratin-dependent area elastic constant
                parameters.at("K")*keffect(vertex.getIndex());
            area.emplace(vertex.getIndex(),
                mesh->getVertexToNeighboursArea(vertex.getIndex()));
            double const A0 = targetArea.at(vertex.getIndex());
            // perimeter elasticity
            double const Gammak =                           // keratin-dependent perimeter elastic constant
                parameters.at("Gamma")*keffect(vertex.getIndex());
            double const perimeter =
                mesh->getVertexToNeighboursPerimeter(vertex.getIndex());
            double const P0 = parameters.at("p0")*sqrt(A0);

            // area pressure
            pressure.emplace(vertex.getIndex(),
                -Kk*(parameters.at("alpha")/parameters.at("A0"))
                    *(area.at(vertex.getIndex()) - A0)
                    *area.at(vertex.getIndex()));
            // perimeter pressure only for non-boundary cells
            bool isBoundaryCell = false;
            std::vector<long int> const halfEdgesToNeighbours =
                mesh->getNeighbourVertices(vertex.getIndex())[1];
            for (long int index : halfEdgesToNeighbours) {  // loop over half-edges to neighbour vertices
                long int const neighbourVertexIndex =       // neighbour cell index
                    halfEdges.at(
                        halfEdges.at(
                            halfEdges.at(
                                halfEdges.at(index).getNextIndex()
                            ).getPairIndex()
                        ).getNextIndex()
                    ).getToIndex();
                if (vertices->at(neighbourVertexIndex).getBoundary()) {
                    isBoundaryCell = true;
                    break;
                }
            }
            if (!isBoundaryCell) {
                pressure[vertex.getIndex()] +=
                -Gammak*(parameters.at("alpha")/parameters.at("A0"))
                    *(perimeter - P0)
                    *perimeter;
            }
            // tension
            tension.emplace(vertex.getIndex(),
                Gammak*(perimeter - P0));

            // neighbours
            std::vector<long int> const neighbourVerticesIndices =
                mesh->getNeighbourVertices(vertex.getIndex())[0];
            long int const numberNeighbours =
                neighbourVerticesIndices.size();

            std::vector<double> radiusToCorner;
            std::vector<double> toPreviousNeighbour, toNextNeighbour;
            std::vector<double> crossToPreviousNeighbour, crossToNextNeighbour;
            double distToPreviousNeighbour, distToNextNeighbour;
            double force;
            for (int i=0; i < numberNeighbours; i++) {  // loop over neighbours
                assert(vertices->at(neighbourVerticesIndices[i]).getType()
                    != "centre");

                // from cell centre to cell corner
                radiusToCorner = mesh->wrapTo(
                    vertex.getIndex(),
                    neighbourVerticesIndices[i]);
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
                    force = (Kk/2.)*(
                        area.at(vertex.getIndex()) - A0)*(
                            crossToPreviousNeighbour[dim]
                                - crossToNextNeighbour[dim]);
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                    // perimeter force
                    force = tension.at(vertex.getIndex())*(
                        toPreviousNeighbour[dim]    // force from previous neighbour
                            /distToPreviousNeighbour);
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                    force = tension.at(vertex.getIndex())*(
                        toNextNeighbour[dim]        // force from next neighbour
                            /distToNextNeighbour);
                    (*forces)[neighbourVerticesIndices[i]][dim] += force;   // force on vertex i
                }
            }
        }

        void addAllForces() override {
            area.clear(); pressure.clear(); tension.clear();
            VertexForce<ForcesType>::addAllForces();
        }

        void integrate(double const& dt) override {

            // integrate keratin concentration
            double const dt_ = dt/parameters.at("tau");
            double const amp = parameters.at("sigma")*sqrt(2.*dt_);
            for (auto it=pressure.begin(); it != pressure.end(); ++it) {
                // keratin concentration
                double const kon =                                  // on-rate...
                    std::max(0., -pressure.at(it->first))           // ... increases with pressure (set by addAllForces)...
                    + parameters.at("tau")*parameters.at("ron");    // ... and time
                double const koff = keratin.at(it->first);          // keratin-dependent off-rate
                keratin[it->first] +=
                    dt_*(kon - koff)                                // deterministic part
                    + amp*random->gauss();                          // stochastic part
            }
            // integrate target area
            for (auto it=area.begin(); it != area.end(); ++it) {
                double const tauk =         // keratin-dependent target area relaxation time
                    parameters.at("taur")*keffect(it->first);
                double const dtk_ = dt/tauk;
                targetArea[it->first] +=    // relax target area
                    -dtk_*(targetArea.at(it->first) - area.at(it->first));
                targetArea[it->first] =     // minimum to target area
                    std::max(parameters.at("A0"), targetArea.at(it->first));
            }
        }

        pybind11::tuple pybind11_getstate() const override;

};

#endif

