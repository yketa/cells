/*
Base types for forces.
*/

#ifndef BASE_FORCES_HPP
#define BASE_FORCES_HPP

#include <map>
#include <string>

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mesh.hpp"

typedef std::map<std::string, double> ParametersType;

// GENERAL BASE

template<class T>
class BaseForce {
/*
All forces.
*/

    protected:

        std::string const type;             // type on which the force applies (if it is "" then it applies to all types)
        ParametersType const parameters;    // force parameters
        T* forces;                          // pointer to forces container

    public:

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

        BaseForce(BaseForce<T>& bF) :   // copy constructor
            BaseForce(bF.getType(), bF.getParameters(), bF.getForces()) {}

        virtual ~BaseForce() {}

        // METHODS

        virtual void integrate(double const& dt) {}
        /*
        Integrate internal degrees of freedom.
        */

        virtual void deleteEdge(TopoChangeEdgeInfoType const& del) {}
        /*
        Update internal degrees of freedom from indices of deleted objects
        following egde deletion (see Mesh::deleteEdge).
        */

        virtual void createEdge(TopoChangeEdgeInfoType const& cre) {}
        /*
        Update internal degrees of freedom from indices of created objects
        following egde creation (see Mesh::createEdge).
        */

        virtual pybind11::tuple pybind11_getstate() const
            { return pybind11::make_tuple(); }
        /*
        Tuple with data to regenerate force computation object.
        First element MUST be an unique identifying string for the force object
        type.
        */

};

// HALF-EDGE FORCES

template<class T>
class HalfEdgeForce : public BaseForce<T> {
/*
Force which derives from the property of a half-edge.
*/

    protected:

        HalfEdgesType* const halfEdges; // pointer to half-edge container

    public:

        HalfEdgeForce(
            std::string const& type_,
            ParametersType const& parameters_,
            T* forces_,
            HalfEdgesType* halfEdges_) :
            BaseForce<T>(type_, parameters_, forces_), halfEdges(halfEdges_) {}

        HalfEdgeForce(HalfEdgeForce<T>& hEF) :  // copy constructor
            BaseForce<T>(hEF), halfEdges(hEF.getHalfEdges()) {}

        HalfEdgesType const getHalfEdges() const
            { return *halfEdges; }

        virtual void addForce(HalfEdge const& halfEdge) {}
        /*
        Add to the force container the forces associated to the half-edge.
        This is the customly defined force.
        */

        virtual void addAllForces() {
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

        VerticesType* const vertices;   // pointer to vertex container

    public:

        VertexForce(
            std::string const& type_,
            ParametersType const& parameters_,
            T* forces_,
            VerticesType* vertices_) :
            BaseForce<T>(type_, parameters_, forces_), vertices(vertices_) {}

        VertexForce(VertexForce<T>& vF) :   // copy constructor
            BaseForce<T>(vF), vertices(vF.getVertices()) {}

        VerticesType* const getVertices() const
            { return vertices; }

        virtual void addForce(Vertex const& vertex) {}
        /*
        Add to the force container the forces associated to the vertex.
        This is the customly defined force.
        */

        virtual void addAllForces() {
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

