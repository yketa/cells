/*
Objects to define the two-dimensional mesh with half-edge procedure.
*/

#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <map>

class Vertex;
class HalfEdge;
class Mesh;

class Vertex {
/*
Individual nodes of the two-dimensional mesh.
*/

    private:

        long int const index;
        double position[2];
        long int halfEdgeIndex;

    public:

        Vertex() : index(-1) {}
        Vertex(long int const& index_, double* const& position_,
            long int const& halfEdgeIndex_=-1)
            : index(index_), position{position_[0], position_[1]},
                halfEdgeIndex(halfEdgeIndex_) {}
        /*
        Parameters
        ----------
        index :
            Unique index for the vertex.
        position :
            Position of the vertex.
        halfEdgeIndex:
            Index of a single halfedge going out of this vertex (the
            association needs not to be bijective).
        */

        long int const& getIndex() { return index; }
        double* getPosition() { return position; }
        long int* getHalfEdgeIndex() { return &halfEdgeIndex; }

};

class HalfEdge {
/*
Directed arrow between vertices (Vertex) of the two-dimensional mesh.
*/

    private:

        long int const index;
        long int fromIndex;
        long int toIndex;
        long int previousIndex;
        long int nextIndex;
        long int pairIndex;

    public:

        HalfEdge() : index(-1) {}
        HalfEdge(long int const& index_,
            long int const& fromIndex_, long int const& toIndex_,
            long int const& previousIndex_=-1,
            long int const& nextIndex_=-1,
            long int const& pairIndex_=-1);
        /*
        Parameters
        ----------
        index :
            Unique index for the half-edge.
        fromIndex :
            Index of the origin vertex.
        toIndex :
            Index of the destination vertex.
        previousIndex :
            Index of the associated `previous' half-edge going towards the
            origin of this half-edge. (default: -1)
        nextIndex :
            Index of the associated `next' half-edge going outwards from the
            destination of this half-edge. (default: -1)
        pairIndex :
            Index of the associated `pair' half-edge going from the destination
            to the origin of this half-edge. (default: -1)
        */

        long int const getIndex() { return index; }
        long int* getFromIndex() { return &fromIndex; }
        long int* getToIndex() { return &toIndex; }
        long int* getPreviousIndex() { return &previousIndex; }
        long int* getNextIndex() { return &nextIndex; }
        long int* getPairIndex() { return &pairIndex; }

};

class Mesh {
/*
Two-dimensional ensembles of vertices and edges.
*/

    private:

        std::map<long int, Vertex> vertices;
        std::map<long int, HalfEdge> halfEdges;
        double systemSize[2] {0, 0};

    public:

        Mesh();

        Vertex* getVertex(long int const& index)
            { return &(vertices[index]); }
        HalfEdge* getHalfEdge(long int const& index)
            { return &(halfEdges[index]); }
        double* getPosition()
            { return systemSize; }

        void wrap(double* position);
        /*
        Wrap position with respect to periodic boundary conditions.

        Parameters
        ----------
        position :
            Pointer to unwrapped position to wrap.
        */

        std::vector<double> wrapDiff(
            double* const& fromPos, double* const& toPos);
        /*
        Wrap difference vector with respect to periodic boundary conditions.

        Parameters
        ----------
        fromPos :
            Pointer to initial point position.
        toPos :
            Pointer to final point position.
        */

        std::vector<double> wrapTo(
            long int const& fromVertexIndex, long int const& toVertexIndex,
            bool const& unit=false);
        /*
        Vector between two vertices, accounting for periodic boundary
        conditions.

        Parameters
        ----------
        fromVertexIndex :
            Index of origin vertex.
        toVertexIndex :
            Index of destination vertex.
        unit :
            Return unitary vector.
        */

        std::vector<double> getHalfEdgeVector(long int const& halfEdgeIndex,
            bool const& unit=false);
        /*
        Vector going from the origin to the destination of a half-edge.

        Parameters
        ----------
        halfEdgeIndex :
            Index of the half-edge.
        unit :
            Return unitary vector.
        */

        double const getEdgeLength(long int const& halfEdgeIndex);
        /*
        Length of edge.

        Parameters
        ----------
        halfEdgeIndex :
            Index of the half-edge.
        */

        std::vector<std::vector<long int>> getNeighbourVertices(
            long int const& vertexIndex);
        /*
        Indices of neighbouring vertices and indices of half-edges towards
        them.

        NOTE: If the half-edge construction is correct, these should be in
              anti-clockwise order.

        Parameters
        ----------
        vertexIndex : int
            Index of the vertex.
        */

        double const getVertexToNeighboursArea(
            long int const& vertexIndex);
        /*
        Area encapsulated by the neighbours of a vertex (shoelace formula).

        Parameters
        ----------
        vertexIndex :
            Index of vertex.
        */

        double const getVertexToNeighboursPerimeter(
            long int const& vertexIndex);
        /*
        Perimeter encapsulated by the neighbours of a vertex.

        Parameters
        ----------
        vertexIndex :
            Index of vertex.
        */

        void checkMesh();
        /*
        Check that the vertices and half-edges define a planar mesh, with
        anticlockwise triangles.
        */

};

#endif

