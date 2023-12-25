/*
Base types for forces.
*/

#ifndef BASE_FORCES_HPP
#define BASE_FORCES_HPP

#include <map>

#include "mesh.hpp"

// GENERAL BASE

template<class T>
class BaseForce {
/*
All forces.
*/

    protected:

        T& forces;  // reference to forces container

    public:

        void reset();
        /*
        Reset force values (and possibly contacts).
        */

        T const& getForces() { return forces; }
        /*
        Return forces.
        */

};

// JUNCTION FORCES

template<class T>
class JunctionForce : public BaseForce<T> {
/*
Force applied on vertices linked by a junction.
*/

    protected:

        std::map<long int, HalfEdge> const& halfEdges;  // reference to half-edge container

    public:

        JunctionForce(
            T& forces_, std::map<long int, HalfEdge> const& halfEdges_)
            : forces(forces_), halfEdges(halfEdges_) {}

        template<typename... Args> void addForce(HalfEdge const& halfEdge,
            Args const&... args);
        /*
        For a halfEdge going from vertex i to vertex j, add the computed force
        to forces[i] and -force to forces[j].
        */

        template<typename... Args> void addAllForces(Args const&... args) {
        /*
        Add forces for all junctions, i.e. half-edges of type "junction".
        */
            std::map<long int, HalfEdge>::const_iterator it;
            for (it=halfEdges.begin(); it != halfEdges.end(); ++it) {
                if ((it->second).getType() == "junction") {
                    addForce(it->second, args);
                }
            }
        }
}

// CELL FORCES

template<class T>
class CellForce : public BaseForce<T> {
/*
Force applied on vertices belonging to a cell.
*/

    protected:

        std::map<long int, Vertex> const& vertices; // reference to vertex container

    public:

        CellForce(
            T& forces_, std::map<long int, Vertex> const& vertices_)
            : forces(forces_), vertices(vertices_) {}

        template<typename... Args> void addForce(Vertex const& vertex,
            Args const&... args);
        /*
        From a cell centre vertex, add force to all neighbours (i.e. cell
        vertices).
        */

        template<typename... Args> void addAllForces(Args const&... args) {
        /*
        Add forces for all cell centres, i.e. vertices of type "centre".
        */
            std::map<long int, Vertex>::const_iterator it;
            for (it=vertices.begin(); it != vertices.end(); ++it) {
                if ((it->second).getType() == "centre") {
                    addForce(it->second, args);
                }
            }
        }
}

// VERTEX FORCES

// PAIRWISE FORCES

template<class T>
class PairwiseForce : public BaseForce<T> {
/*
Force between pair of individuals.
*/

    protected:

        bool const contacts;                // compute number of contacts
        std::vector<int> nContacts(0);      // number of contacts of each individual
        std::vector<int> sumNContacts(0);   // sum of the elements of nContacts

    public:

        PairwiseForce(T& forces_, bool const& contacts_=false) :
            forces(forces_), contacts(contacts_) {}

        template<typename... Args> void addForce(int const& i, int const& j,
            bool const& reciprocal, Args const&... args);
        /*
        Add to forces[i] the force exerted by j on i, and if reciprocal force
        add -forces[i] to forces[j].
        Returns true if the force is non-zero and false otherwise (so that it
        may be used to build contact matrices).
        */

        void computeSumContacts() {
        /*
        Compute sum of the elements of nContacts.
        */
        #ifdef _OPENMP
        #pragma omp master  // this should be performed only once by the master thread
        #endif
            {
                for (int i=0; i < (int) nContacts.size(); i++) {
                    sumNContacts[i + 1] = sumNContacts[i] + nContacts[i];
                }
            } 
        }

        std::vector<int> const& getNContacts() { return nContacts; }
        /* Return number of contacts of each individual. */
        std::vector<int> const& getSumNContacts() { return sumNContacts; }
        /* Return sum of the elements of nContacts. */
        
};

// INDIVIDAL FORCES

template<class T>
class IndividualForce : public BaseForce<T> {
/*
Force exerted on each individual.
*/

    public:

        IndividualForce(T& forces_) : forces(forces_) {}

        template<typename... Args> void addForce(int const& i,
            Args const&... args);
        /*
        Add to forces[i] the force exerted on i.
        */

};

#endif

