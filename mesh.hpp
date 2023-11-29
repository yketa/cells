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

        std::vector<double> position;   // position wrapped with respect to periodic boundary conditions
        std::vector<double> uposition;  // unwrapped position
        long int halfEdgeIndex;

    public:

        Vertex() : index(-1) {}
        Vertex(long int const& index_,
            std::vector<double> const& position_,
            long int const& halfEdgeIndex_=-1) :
                index(index_),
                position(position_), uposition(position_),
                halfEdgeIndex(halfEdgeIndex_) {}
        /*
        Parameters
        ----------
        index_ :
            Unique index for the vertex.
        position_ :
            Position of the vertex.
        halfEdgeIndex_ :
            Index of a single halfedge going out of this vertex (the
            association needs not to be bijective).
        */

        Vertex(long int const& index_,
            std::vector<double> const& position_,
            std::vector<double> const& uposition_,
            long int const& halfEdgeIndex_=-1) :
                index(index_),
                position(position_), uposition(uposition_),
                halfEdgeIndex(halfEdgeIndex_) {}

        long int getIndex() const
            { return index; }
        std::vector<double> getPosition() const
            { return position; }
        std::vector<double> getUPosition() const
            { return uposition; }
        long int getHalfEdgeIndex() const
            { return halfEdgeIndex; }

        void setPosition(std::vector<double> const& r)
            { position = r; }
        void setUPosition(std::vector<double> const& r)
            { uposition = r; }
        void setHalfEdgeIndex(long int const& i)
            { halfEdgeIndex = i; }

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
            long int const& previousIndex_=-1, long int const& nextIndex_=-1,
            long int const& pairIndex_=-1)
            : index(index_),
                fromIndex(fromIndex_), toIndex(toIndex_),
                previousIndex(previousIndex_), nextIndex(nextIndex_),
                pairIndex(pairIndex_) {}
        /*
        Parameters
        ----------
        index_ :
            Unique index for the half-edge.
        fromIndex_ :
            Index of the origin vertex.
        toIndex_ :
            Index of the destination vertex.
        previousIndex_ :
            Index of the associated `previous' half-edge going towards the
            origin of this half-edge. (default: -1)
        nextIndex_ :
            Index of the associated `next' half-edge going outwards from the
            destination of this half-edge. (default: -1)
        pairIndex_ :
            Index of the associated `pair' half-edge going from the destination
            to the origin of this half-edge. (default: -1)
        */

        long int getIndex() const { return index; }
        long int getFromIndex() const { return fromIndex; }
        long int getToIndex() const { return toIndex; }
        long int getPreviousIndex() const { return previousIndex; }
        long int getNextIndex() const { return nextIndex; }
        long int getPairIndex() const { return pairIndex; }

        void setFromIndex(long int const& i) { fromIndex = i; }
        void setToIndex(long int const& i) { toIndex = i; }
        void setPreviousIndex(long int const& i) { previousIndex = i; }
        void setNextIndex(long int const& i) { nextIndex = i; }
        void setPairIndex(long int const& i) { pairIndex = i; }

};

class Mesh {
/*
Two-dimensional ensembles of vertices and edges.
*/

    protected:

        std::map<long int, Vertex> vertices;
        std::map<long int, HalfEdge> halfEdges;
        std::vector<double> systemSize{0, 0};

    public:

        Mesh() {}

        Mesh(                   // used to load state
            std::map<long int, Vertex> const& vertices_,
            std::map<long int, HalfEdge> const& halfEdges_,
            std::vector<double> const& systemSize_) {
            // clear maps and vector
            vertices.clear();
            halfEdges.clear();
            systemSize.clear();
            // copy maps
            vertices.insert(vertices_.begin(), vertices_.end());
            halfEdges.insert(halfEdges_.begin(), halfEdges_.end());
            // set vector
            systemSize.push_back(systemSize_[0]);
            systemSize.push_back(systemSize_[1]);
        }
        Mesh(Mesh const& m_) :  // copy constructor
            Mesh(m_.getVertices(), m_.getHalfEdges(), m_.getSystemSize()) {}

        std::map<long int, Vertex> getVertices() const
            { return vertices; }
        std::map<long int, HalfEdge> getHalfEdges() const
            { return halfEdges; }
        std::vector<double> getSystemSize() const
            { return systemSize; }

        std::vector<double> wrap(
            std::vector<double> const& position)
            const;
        /*
        Wrap position to positive values with respect to periodic boundary
        conditions.

        Parameters
        ----------
        position :
            Position to wrap.

        Returns
        -------
        wposition :
            Wraped position.
        */

        std::vector<double> wrapDiff(
            std::vector<double> const& fromPos,
            std::vector<double> const& toPos)
            const;
        /*
        Wrap difference vector with respect to periodic boundary conditions.

        Parameters
        ----------
        fromPos :
            Pointer to initial point position.
        toPos :
            Pointer to final point position.

        Returns
        -------
        disp :
            (2,) difference vector.
        */

        std::vector<double> wrapTo(
            long int const& fromVertexIndex, long int const& toVertexIndex,
            bool const& unit=false)
            const;
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

        Returns
        -------
        fromTo :
            (2,) vector between vertices.
        */

        std::vector<double> getHalfEdgeVector(
            long int const& halfEdgeIndex, bool const& unit=false)
            const;
        /*
        Vector going from the origin to the destination of a half-edge.

        Parameters
        ----------
        halfEdgeIndex :
            Index of the half-edge.
        unit :
            Return unitary vector.

        Returns
        -------
        halfEdgeVector :
            (2,) vector corresponding to half-edge.
        */

        double getEdgeLength(
            long int const& halfEdgeIndex)
            const;
        /*
        Length of edge.

        Parameters
        ----------
        halfEdgeIndex :
            Index of the half-edge.

        Returns
        -------
        edgeLength :
            Length of edge.
        */

        std::vector<std::vector<long int>> getNeighbourVertices(
            long int const& vertexIndex)
            const;
        /*
        Indices of neighbouring vertices and indices of half-edges towards
        them.

        NOTE: If the half-edge construction is correct, these should be in
              anti-clockwise order.

        Parameters
        ----------
        vertexIndex : int
            Index of the vertex.


        Returns
        -------
        neighbourVerticesIndices :
            Indices of vertices neighbouring this vertex.
        halfEdgesToNeighboursIndices :
            Indices of half-edges from this vertex towards neighbour vertices.
        */

        double getVertexToNeighboursArea(
            long int const& vertexIndex)
            const;
        /*
        Area encapsulated by the neighbours of a vertex (shoelace formula).

        Parameters
        ----------
        vertexIndex :
            Index of vertex.

        Returns
        -------
        area :
            Area encapsulated by neighbours.
        */

        double getVertexToNeighboursPerimeter(
            long int const& vertexIndex)
            const;
        /*
        Perimeter encapsulated by the neighbours of a vertex.

        Parameters
        ----------
        vertexIndex :
            Index of vertex.

        Returns
        -------
        perimeter :
            Perimeter encapsulated by neighbours.
        */

        void checkMesh() const;
        /*
        Check that the vertices and half-edges define a planar mesh, with
        anticlockwise triangles.
        */

};

#endif

