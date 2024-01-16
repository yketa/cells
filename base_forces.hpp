/*
Base types for forces.
*/

#ifndef BASE_FORCES_HPP
#define BASE_FORCES_HPP

#include <map>
#include <string>

#include "mesh.hpp"

typedef std::map<long int, HalfEdge> HalfEdgesType;
typedef std::map<long int, Vertex> VerticesType;
typedef std::map<std::string, double> ParametersType;

// GENERAL BASE

template<class T>
class BaseForce {
/*
All forces.
*/

    protected:

        std::string const type;                         // type on which the force applies (if it is "" then it applies to all types)
        std::map<std::string, double> const parameters; // force parameters
        T* forces;                                      // reference to forces container

    public:

        // CLEAR

        void clear();
        /*
        Clear force values.
        */

        // GETTERS

        std::string const getType() const
            { return type; }
        ParametersType const getParameters() const
            { return parameters; }
        T* getForces()
            { return forces; }

        // CONSTRUCTORS AND DESTRUCTORS

        BaseForce(
            std::string const& type_,
            ParametersType const& parameters_,
            T* forces_) :
            type(type_), parameters(parameters_), forces(forces_) {}

        BaseForce(BaseForce<T> const& bF) : // copy constructor
            BaseForce(bF.getType(), bF.getParameters(), bF.getForces()) {}

        virtual ~BaseForce() = 0;

        // METHODS

        virtual void integrate(double const& dt) {}
        /*
        Integrate internal degrees of freedom.
        */

};

// HALF-EDGE FORCES

template<class T>
class HalfEdgeForce : public BaseForce<T> {
/*
Force which derives from the property of a half-edge.
*/

    protected:

        HalfEdgesType* const halfEdges;  // reference to half-edge container

    public:

        HalfEdgeForce(
            std::string const& type_,
            ParametersType const& parameters_,
            T* forces_,
            HalfEdgesType* halfEdges_) :
            BaseForce<T>(type_, parameters_, forces_), halfEdges(halfEdges_) {}

        HalfEdgeForce(HalfEdgeForce<T> const& hEF) :    // copy constructor
            BaseForce<T>(hEF), halfEdges(hEF.getHalfEdges()) {}

        HalfEdgesType const getHalfEdges() const
            { return *halfEdges; }

        virtual void addForce(HalfEdge const& halfEdge) {}
        /*
        Add to the force container the forces associated to the half-edge.
        This is the customly defined force.
        */

        void addAllForces() {
        /*
        Add to the force container the forces associated to all half-edges.
        */
            HalfEdgesType::const_iterator it;
            for (it=halfEdges->begin(); it != halfEdges->end(); ++it) {
                if ((it->second).getType() == this->type
                    || this->type == "") {
                    addForce(it->second);
                }
            }
        }

};

// VERTEX FORCES

template<class T>
class VertexForce : public BaseForce<T> {
/*
Force which derives from the property of a vertex.
*/

    protected:

        VerticesType* const vertices; // reference to vertex container

    public:

        VertexForce(
            std::string const& type_,
            ParametersType const& parameters_,
            T* forces_,
            VerticesType* vertices_) :
            BaseForce<T>(type_, parameters_, forces_), vertices(vertices_) {}

        VertexForce(VertexForce<T> const& vF) : // copy constructor
            BaseForce<T>(vF), vertices(vF.getVertices()) {}

        VerticesType const getVertices() const
            { return *vertices; }

        virtual void addForce(Vertex const& vertex) {}
        /*
        Add to the force container the forces associated to the vertex.
        This is the customly defined force.
        */

        void addAllForces() {
        /*
        Add to the force container the forces associated to all vertices.
        */
            VerticesType::const_iterator it;
            for (it=vertices->begin(); it != vertices->end(); ++it) {
                if ((it->second).getType() == this->type
                    || this->type == "") {
                    addForce(it->second);
                }
            }
        }

};

#endif

