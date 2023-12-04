/*
Objects for simulation.
*/

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "mesh.hpp"
#include "random.hpp"
#include "tools.hpp"

class SPVertex;
class Cell;
class Face;
class Junction;
class VertexModel;

/*
 *  PHYSICS OBJECTS
 *
 */

class SPVertex {
/*
Physical self-propelled vertex.
*/

    private:

        long int const vertexIndex;
        double const v0 = 0;    // self-propulsion velocity
        double const Dr = 0;    // rotational diffusion constant
        double theta = 0;       // angle of self-propulsion

    public:

        SPVertex() : vertexIndex(-1) {}
        SPVertex(long int const& vertexIndex_, double const& theta_=0,
            double const& v0_=1e-1, double const& Dr_=1e-1)
            : vertexIndex(vertexIndex_), v0(v0_), Dr(Dr_), theta(theta_) {}
        /*
        Parameters
        ----------
        vertexIndex_ :
            Unique index for the self-propelled vertex identical to index of
            vertex.
        theta_ :
            Initial self-propulsion direction.
        v0_ :
            Self-propulsion velocity.
        Dr_ :
            Rotational diffusion constant.
        */

        long int getVertexIndex() const { return vertexIndex; }
        double getv0() const { return v0; }
        double getDr() const { return Dr; }
        double gettheta() const { return theta; }

        void settheta(double const& t) { theta = t; }

};

class Cell {
/*
Physical cell enclosed by junctions (HalfEdge) between vertices (Vertex) of the
mesh.
*/

    private:

        long int const vertexIndex;

        double area = -1;       // area of cell
        double const kA = 0;    // area stiffness
        double const A0 = 0;    // target area of cell
        double perimeter = -1;  // perimeter of cell
        double const kP = 0;    // perimeter stiffness
        double const p0 = 0;    // dimensionless target perimeter of cell

    public:

        Cell() : vertexIndex(-1) {}
        Cell(long int const& vertexIndex_,
            double const& A0_, double const& p0_=3.81,
            double const& kA_=1, double const& kP_=1)
            : vertexIndex(vertexIndex_), kA(kA_), A0(A0_), kP(kP_), p0(p0_) {}
        /*
        Parameters
        ----------
        vertexIndex_ :
            Unique index for the cell identical to index of vertex.
        A0_ :
            Target area of cell.
        p0_ :
            Dimensionless target perimeter of cell.
        kA_ :
            Area stiffness.
        kP_ :
            Perimeter stiffness.
        */

        long int getVertexIndex() const { return vertexIndex; }
        double getArea() const { return area; }
        double getkA() const { return kA; }
        double getA0() const { return A0; }
        double getPerimeter() const { return perimeter; }
        double getkP() const { return kP; }
        double getp0() const { return p0; }
        double getP0() const { return p0*std::sqrt(A0); }

        void setArea(double const& a) { area = a; }
        void setPerimeter(double const& p) { perimeter = p; }

};

class Face {
/*
Triangle enclosed by three vertices.
*/

    private:

        long int const halfEdgeIndex;

    public:

        Face() : halfEdgeIndex(-1) {}
        Face(long int const& halfEdgeIndex_) : halfEdgeIndex(halfEdgeIndex_) {}
        /*
        Parameters
        ----------
        halfEdgeIndex_ :
            Unique index of a helf-edge belonging to the face.
        */

        long int getHalfEdgeIndex() const { return halfEdgeIndex; }

};

class Junction {
/*
Physical link between two vertices.
*/

    private:

        long int const halfEdgeIndex;

    public:

        Junction() : halfEdgeIndex(-1) {}
        Junction(long int const& halfEdgeIndex_)
            : halfEdgeIndex(halfEdgeIndex_) {}
        /*
        Parameters
        ----------
        halfEdgeIndex_ :
            Unique index of a helf-edge belonging to the junction.
        */

        long int getHalfEdgeIndex() const { return halfEdgeIndex; }

};

/*
 *  SIMULATION OBJECTS
 *
 */

class VertexModel : public Mesh {

    private:

        /*
        inherited from Mesh
        -------------------
        std::map<long int, Vertex> vertices;    // std::map<long int, Vertex> Mesh::getVertices
        std::map<long int, HalfEdge> halfEdges; // std::map<long int, HalfEdge> Mesh::getHalfEdges
        std::vector<double> systemSize;         // std::vector<double> Mesh::getSystemSize
        */

        MultiIntKeyDict<SPVertex> sPVertices;
        MultiIntKeyDict<Cell> cells;
        MultiIntKeyDict<Face> faces;
        MultiIntKeyDict<Junction> junctions;

        double const v0;    // vertex self-propulsion velocity
        double const Dr;    // vertex propulsion rotational diffusion constant
        double const p0;    // dimensionless target perimeter of cell

        long int const seed;    // random number generator seed
        Random random;      	// random number generator

        double time = 0;
        long int nT1 = 0;   // number of T1s

    public:

        VertexModel(long int const& seed_=0,
            double const& v0_=1e-1, double const& Dr_=1e-1,
            double const& p0_=3.81,
            bool const& boundary_=false)
            : Mesh(boundary_), v0(v0_), Dr(Dr_), p0(p0_),
                seed(seed_), random(seed) {}
        /*
        Parameters
        ----------
        seed :
            Random number generator seed.
        v0_ :
            Vertex self-propulsion velocity.
        Dr_ :
            Vertex propulsion rotational diffusion constant.
        p0_ :
            Dimensionless target perimeter of cell.
        */

        VertexModel(                            // used to load state
            Mesh const& mesh_,
            MultiIntKeyDict<SPVertex> const& sPVertices_,
            MultiIntKeyDict<Cell> const& cells_,
            MultiIntKeyDict<Face> const& faces_,
            MultiIntKeyDict<Junction> const& junctions_,
            double const& v0_, double const& Dr_, double const& p0_,
            long int const& seed_, double const time_, long int const nT1_) :
            // geometrical objects (mesh)
            Mesh(mesh_),
            // physical objects
            sPVertices(sPVertices_), cells(cells_), faces(faces_),
                junctions(junctions_),
            // physical parameters
            v0(v0_), Dr(Dr_), p0(p0_),
            // integration quantities
            seed(seed_), time(time_), nT1(nT1_) {}
        VertexModel(VertexModel const& vm_) :   // copy constructor
            VertexModel(
                // geometrical objects (mesh)
                vm_,
                // physical objects
                vm_.getSPVertices(), vm_.getCells(), vm_.getFaces(),
                    vm_.getJunctions(),
                // physical parameters
                vm_.getv0(), vm_.getDr(), vm_.getp0(),
                // integration quantities
                vm_.getSeed(), vm_.getTime(), vm_.getnT1()) {}

        MultiIntKeyDict<SPVertex> getSPVertices() const { return sPVertices; }
        MultiIntKeyDict<Cell> getCells() const { return cells; }
        MultiIntKeyDict<Face> getFaces() const { return faces; }
        MultiIntKeyDict<Junction> getJunctions() const { return junctions; }

        double getv0() const { return v0; }
        double getDr() const { return Dr; }
        double getp0() const { return p0; }

        long int getSeed() const { return seed; }

        double getTime() const { return time; }
        long int getnT1() const { return nT1; }

        void reset() {
        /*
        Clear all data.
        */
            vertices.clear();
            halfEdges.clear();
            systemSize.clear();
            systemSize.push_back(0); systemSize.push_back(0);
            sPVertices.clear();
            cells.clear();
            faces.clear();
            junctions.clear();
            time = 0;
            nT1 = 0;
        }

        void integrate(double const& dt=0,
            double const& delta=0.1, double const& epsilon=0.1);
        /*
        Integrate one time step of self-propelled vertex model and check for
        and perform T1s.

        Parameters
        ----------
        dt :
            Integration time step.
        delta :
            Distance between vertices below which these should be merged.
        epsilon :
            Create two vertices at distance `delta' + `epsilon' after T1.
        */

        std::map<long int,std::vector<double>> const getForces();
        /*
        Compute forces on vertices from vertex model.

        Returns
        -------
        forces :
            Dictionary which associates vertex indices to force applied on the
            vertex.
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

        long int const mergeVertices(long int const& halfEdgeIndex);
        /*
        Merge two vertices.

        Parameters
        ----------
        halfEdgeIntex :
            Index of half-edge linking two vertices to be merged.

        Returns
        -------
        fromMergeIndex :
            Merged vertex index.
        */

        long int const createJunction(
            long int const& halfEdgeIndex0, long int const& halfEdgeIndex1,
            double const& angle, double const& length=1);
        /*
        Create a new vertex and junction.

        Parameters
        ----------
        halfEdgeIndex0 :
            Index of first half-edge going out of a vertex, and from whose pair
            half-edge it will be separated after the introduction of a new
            junction.
        halfEdgeIndex1 :
            Index of second half-edge going out of the same vertex, and from
            whose pair half-edge it will be separated after the introduction of
            a new junction.
        angle :
            Angle of the new junction with respect to the horizontal axis.
        distance :
            Length to set for the new junction.

        Returns
        -------
        newVertexIndex :
            New vertex index.
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

};

#endif

