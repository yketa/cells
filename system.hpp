/*
Objects for simulation.
*/

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "class_factory.hpp"
#include "forces.hpp"
#include "mesh.hpp"
#include "random.hpp"
#include "tools.hpp"

typedef std::map<long int, std::vector<double>> ForcesType;

class VertexModel : public Mesh {

    private:

        /*
        [mesh.hpp]
        inherited from Mesh
        -------------------
        std::map<long int, Vertex> vertices;    // std::map<long int, Vertex> Mesh::getVertices
        std::map<long int, HalfEdge> halfEdges; // std::map<long int, HalfEdge> Mesh::getHalfEdges
        std::vector<double> systemSize;         // std::vector<double> Mesh::getSystemSize
        */

        ClassFactory<HalfEdgeForce<ForcesType>> halfEdgeForces; // forces deriving from half-edge properties
        ClassFactory<VertexForce<ForcesType>> vertexForces;     // forces deriving from vertex properties
        ForcesType forces;                                      // force applied on each vertex

        long int const seed;    // random number generator seed
        Random random;      	// random number generator

        double time = 0;    // internal time
        long int nT1 = 0;   // number of T1s

    public:

        // CLEAR

        void clear() {
        /*
        Clear all data.
        */
            Mesh::clear();
            time = 0;
            nT1 = 0;
        }

        // GETTERS AND SETTERS

        ClassFactory<HalfEdgeForce<ForcesType>> const& getHalfEdgeForces()
            const
            { return halfEdgeForces; }
        ClassFactory<VertexForce<ForcesType>> const& getVertexForces()
            const
            { return vertexForces; }
        ForcesType const& getForces()
            const
            { return forces; }

        long int getSeed() const { return seed; }
        Random const& getRandom() const { return random; }

        double getTime() const { return time; }
        long int getnT1() const { return nT1; }

        void copyRandom(VertexModel const& vm)
            { random.setGenerator(vm.getRandom().getGenerator()); }

        // CONSTRUCTORS

        VertexModel(long int const& seed_=0) : seed(seed_), random(seed) {}
        /*
        Parameters
        ----------
        seed_ :
            Random number generator seed.
        */

        VertexModel(                            // used to initiate state
            Mesh const& mesh_,
//             ClassFactory<HalfEdgeForce<ForcesType>> const& halfEdgeForces_,
//             ClassFactory<VertexForce<ForcesType>> const& vertexForces_,
            long int const& seed_, Random const& random_,
            double const time_, long int const nT1_) :
            // geometrical objects (Mesh)
            Mesh(mesh_),
//             // forces objects
//             halfEdgeForces(halfEdgeForces_), vertexForces(vertexForces_),
            // integration quantities
            seed(seed_), random(random_), time(time_), nT1(nT1_) {}

        VertexModel(VertexModel const& vM) :   // copy constructor
            // geometrical objects (Mesh)
            Mesh(vM),
            // force objects
            halfEdgeForces(vM.getHalfEdgeForces()),
            vertexForces(vM.getVertexForces()),
            // integration quantities
            seed(vM.getSeed()),
            random(vM.getRandom()),
            time(vM.getTime()),
            nT1(vM.getnT1()) {}

        // FORCE "SETTERS"

        template<class ForceType, typename... Args> void addHalfEdgeForce(
            std::string const& name, Args ...args);
        void removeHalfEdgeForce(
            std::string const& name)
            { halfEdgeForces.remove(name); }

        template<class ForceType, typename... Args> void addVertexForce(
            std::string const& name, Args ...args);
        void removeVertexForce(
            std::string const& name)
            { vertexForces.remove(name); }

        // METHODS

        void integrate(double const& dt=0,
            double const& delta=0.1, double const& epsilon=0.1);
        /*
        Compute forces, integrate one time step, and perform T1s.

        Parameters
        ----------
        dt :
            Integration time step.
        delta :
            Distance between vertices below which these should be merged.
        epsilon :
            Create two vertices at distance `delta' + `epsilon' after T1.
        */

        void computeForces();
        /*
        Compute forces on vertices from forces objects.
        */

        void doT1(double const& delta=0.1, double const& epsilon=0.1);
        /*
        Check all junctions and perform T1s.

        Parameters
        ----------
        delta :
            Distance between vertices below which these should be merge.
        epsilon :
            Create two vertices at distance `delta' + `epsilon'.
        */

        void checkMesh(
            std::vector<std::string> helfEdgeTypes=std::vector<std::string>())
            const;
        /*
        Check mesh with Mesh::checkMesh and that types in `helfEdgeTypes' are
        not defined identically for both half-edges of a single edge.
        */

        void initRegularTriangularLattice(
                long int const& size=6, double const& junctionLength=1);
        /*
        Initialise a regular triangular lattice.

        Parameters
        ----------
        size :
            Number of vertices in both horizontal and vertical directions.
            NOTE: Must be a multiple of 6.
        junctionLength :
            Length of nearest neighbour junctions.
        */

        void initOpenRegularTriangularLattice(
                long int const& size=6, double const& junctionLength=1);
        /*
        Initialise a regular triangular lattice with a cell replaced with a
        hole. (TESTING)

        Parameters
        ----------
        size :
            Number of vertices in both horizontal and vertical directions.
            NOTE: Must be a multiple of 6.
        junctionLength :
            Length of nearest neighbour junctions.
        */

        void initOpenRegularHexagonalLattice(
                long int const& nCells=1, double const& junctionLength=1);
        /*
        Initialise a regular square lattice with open outer bondary.

        Parameters
        ----------
        nCells :
            Number of cells.
            NOTE: Must be the square of an integer.
        junctionLength :
            Length of nearest neighbour junctions.
        */

};

#endif

