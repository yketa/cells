/*
Objects for simulation.
*/

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "class_factory.hpp"
#include "forces.hpp"
#include "integrators.hpp"
#include "mesh.hpp"
#include "random.hpp"
#include "tools.hpp"

class VertexModel : public Mesh {

    private:

        /*
        [mesh.hpp]
        inherited from Mesh
        -------------------
        VerticesType vertices;          // VerticesType Mesh::getVertices
        HalfEdgesType halfEdges;        // HalfEdgesType Mesh::getHalfEdges
        std::vector<double> systemSize; // std::vector<double> Mesh::getSystemSize
        */

        ForcesType forces;                                      // force applied on each vertex
        ClassFactory<HalfEdgeForce<ForcesType>> halfEdgeForces; // forces deriving from half-edge properties
        ClassFactory<VertexForce<ForcesType>> vertexForces;     // forces deriving from vertex properties

        VelocitiesType velocities;                                              // velocities
        std::shared_ptr<BaseIntegrator<ForcesType, VelocitiesType>> integrator; // integrator

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

            forces.clear();
            halfEdgeForces.clear();
            vertexForces.clear();

            velocities.clear();
            integrator = std::make_shared<UnitOverdamped>(
                &forces, &velocities);

            time = 0;
            nT1 = 0;
        }

        // GETTERS AND SETTERS

        ForcesType const& getForces()
            const
            { return forces; }
        ClassFactory<HalfEdgeForce<ForcesType>> const& getHalfEdgeForces()
            const
            { return halfEdgeForces; }
        ClassFactory<VertexForce<ForcesType>> const& getVertexForces()
            const
            { return vertexForces; }

        VelocitiesType const& getVelocities()
            const
            { return velocities; }
        std::shared_ptr<BaseIntegrator<ForcesType, VelocitiesType>> const&
            getIntegrator() const
            { return integrator; }

        long int const& getSeed() const { return seed; }
        Random const& getRandom() const { return random; }
        void setRandom(Random const& random_) { random = random_; }

        double const& getTime() const { return time; }
        long int const& getnT1() const { return nT1; }

        // CONSTRUCTORS

        VertexModel(long int const& seed_=0) :
            // default integrator
            integrator(std::make_shared<UnitOverdamped>(&forces, &velocities)),
            // integration quantities
            seed(seed_), random(seed) {}

        VertexModel(                            // used to initiate state
            Mesh const& mesh_,
            VelocitiesType const& velocities_,
            long int const& seed_, double const time_, long int const nT1_) :
            // geometrical objects (Mesh)
            Mesh(mesh_),
            // velocities
            velocities(velocities_),
            // default integrator
            integrator(std::make_shared<UnitOverdamped>(&forces, &velocities)),
            // integration quantities
            seed(seed_), random(seed), time(time_), nT1(nT1_) {}

        VertexModel(
            VertexModel const& vertexModel_, long int const& seed_) :
            // geometrical objects (Mesh)
            Mesh(vertexModel_),
            // velocities
            velocities(vertexModel_.getVelocities()),
            // default integrator
            integrator(std::make_shared<UnitOverdamped>(&forces, &velocities)),
            // integration quantities
            seed(seed_), random(seed),
            time(vertexModel_.getTime()), nT1(vertexModel_.getnT1()) {}

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

        // INTEGRATOR "SETTERS"

        template<class IntegratorType, typename... Args> void setIntegrator(
            Args ...args);

        // METHODS

        void integrate(double const& dt=0,
            double const& delta=0.1, double const& epsilon=0.1);
        /*
        Compute forces, integrate one time step, and perform T1s.

        WARNING: This explicitly move cell centres to match the centre of mass
                 of cell corners.

        Parameters
        ----------
        dt :
            Integration time step.
        delta :
            Distance between vertices below which these should be merged.
        epsilon :
            Create two vertices at distance `delta' + `epsilon' after T1.
        */

        void integrateVelocities(double const& dt);
        /*
        Compute forces on vertices from forces objects and integrate velocities
        from integrator.

        WARNING: This explicitly remove centre of mass force.
        */

        void doT1(double const& delta=0.1, double const& epsilon=0.1);
        /*
        Check all junctions and perform T1s and update internal degrees of
        freedom.

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
            long int const& size=6, double const& hexagonArea=1);
        /*
        Initialise a regular triangular lattice.

        Parameters
        ----------
        size :
            Number of vertices in both horizontal and vertical directions.
            NOTE: Must be a multiple of 6.
        hexagonArea :
            Area of regular hexagons.
        */

        void initOpenRegularTriangularLattice(
            long int const& size=6, double const& hexagonArea=1);
        /*
        Initialise a regular triangular lattice with a cell replaced with a
        hole. (TESTING)

        Parameters
        ----------
        size :
            Number of vertices in both horizontal and vertical directions.
            NOTE: Must be a multiple of 6.
        hexagonArea :
            Area of regular hexagons.
        */

        void initOpenRegularHexagonalLattice(
            long int const& nCells=1, double const& hexagonArea=1,
            double const& boxLength=3);
        /*
        Initialise a regular square lattice with open outer bondary.

        Parameters
        ----------
        nCells :
            Number of cells.
            NOTE: Must be the square of an integer.
        hexagonArea :
            Area of regular hexagons.
        boxLength :
            Length of the box in units of sqrt(`nCells')*junctionLength where
            junctionLength is the length of a side of a regular hexagon.
        */

};

#endif

