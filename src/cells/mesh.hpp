/*
Objects to define the two-dimensional mesh with half-edge procedure.
*/

#ifndef MESH_HPP
#define MESH_HPP

#include <map>
#include <string>
#include <vector>

class Vertex;
class HalfEdge;
class Mesh;

typedef std::map<long int, Vertex> VerticesType;
typedef std::map<long int, HalfEdge> HalfEdgesType;
typedef std::tuple<long int, std::vector<long int>> TopoChangeEdgeInfoType;

class Vertex {
/*
Individual nodes of the two-dimensional mesh.
*/

    private:

        long int const index;           // unique index of vertex
        bool const boundary;            // designates an outer boundary
        std::string const type;         // type of the vertex

        std::vector<double> position;   // position wrapped with respect to periodic boundary conditions
        std::vector<double> uposition;  // unwrapped position
        long int halfEdgeIndex;

    public:

        Vertex() : index(-1), boundary(false), type("") {}

        Vertex(long int const& index_,
            std::vector<double> const& position_,
            long int const& halfEdgeIndex_=-1,
            bool const& boundary_=false,
            std::string const& type_="") :
                index(index_), boundary(boundary_), type(type_),
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
        boundary_ :
            Indicates that vertex is used to designate an outer boundary.
        type_ :
            Name of the type of vertex.
        */

        Vertex(long int const& index_,
            std::vector<double> const& position_,
            std::vector<double> const& uposition_,
            long int const& halfEdgeIndex_,
            bool const& boundary_,
            std::string const& type_) :
                index(index_), boundary(boundary_), type(type_),
                position(position_), uposition(uposition_),
                halfEdgeIndex(halfEdgeIndex_) {}

        long int getIndex() const
            { return index; }
        bool getBoundary() const
            { return boundary; }
        std::string getType() const
            { return type; }
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

        long int const index;   // unique index of the half-edge
        std::string const type; // type of the half-edge

        long int fromIndex;     // index of the origin vertex of the half-edge
        long int toIndex;       // index of the destination vertex of the half-edge
        long int previousIndex; // index of the previous half-edge in the triangle
        long int nextIndex;     // index of the next half-edge in the triangle
        long int pairIndex;     // index of the pair half-edge in the adjacent triangle at this half-edge

    public:

        HalfEdge() : index(-1), type("") {}
        HalfEdge(long int const& index_,
            long int const& fromIndex_, long int const& toIndex_,
            long int const& previousIndex_=-1, long int const& nextIndex_=-1,
            long int const& pairIndex_=-1,
            std::string const& type_="")
            : index(index_), type(type_),
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
            origin of this half-edge.
        nextIndex_ :
            Index of the associated `next' half-edge going outwards from the
            destination of this half-edge.
        pairIndex_ :
            Index of the associated `pair' half-edge going from the destination
            to the origin of this half-edge.
        type_ :
            Name of the type of vertex.
        */

        long int getIndex() const { return index; }
        std::string getType() const { return type; }
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

        VerticesType vertices;
        HalfEdgesType halfEdges;
        std::vector<double> systemSize{0, 0};

    public:

        // CLEAR

        void clear() {
        /*
        Clear all data.
        */
            vertices.clear();
            halfEdges.clear();
            systemSize.clear();
            systemSize.push_back(0); systemSize.push_back(0);
        }

        // GETTERS

        VerticesType const& getVertices() const
            { return vertices; }
        HalfEdgesType const& getHalfEdges() const
            { return halfEdges; }
        std::vector<double> const& getSystemSize() const
            { return systemSize; }
        std::vector<double> getCentre() const
            { return {systemSize[0]/2., systemSize[1]/2.}; }

        std::vector<long int> getVertexIndicesByType(
            std::string const& type, bool const& exclude_boundary=true) const {
            std::vector<long int> vertexIndices;
            for (auto it=vertices.begin(); it != vertices.end(); ++it) {
                if ((it->second).getType() == type || type == "") {
                    if (exclude_boundary && (it->second).getBoundary())
                        { continue; }
                    vertexIndices.push_back(it->first);
                }
            }
            return vertexIndices;
        }

        std::vector<long int> getHalfEdgeIndicesByType(
            std::string const& type) const {
            std::vector<long int> halfEdgeIndices;
            for (auto it=halfEdges.begin(); it != halfEdges.end(); ++it) {
                if ((it->second).getType() == type || type == "") {
                    halfEdgeIndices.push_back(it->first);
                }
            }
            return halfEdgeIndices;
        }

        // CONSTRUCTORS

        Mesh() {}

        Mesh(                   // used to load state
            VerticesType const& vertices_,
            HalfEdgesType const& halfEdges_,
            std::vector<double> const& systemSize_) {
            // clear
            clear();
            // copy maps
            vertices.insert(vertices_.begin(), vertices_.end());
            halfEdges.insert(halfEdges_.begin(), halfEdges_.end());
            // set vector
            systemSize[0] = systemSize_[0];
            systemSize[1] = systemSize_[1];
        }
        Mesh(Mesh const& m_) :  // copy constructor
            Mesh(m_.getVertices(), m_.getHalfEdges(), m_.getSystemSize()) {}

        // METHODS

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
            std::vector<double> const& toPos,
            bool const& unit=false)
            const;
        /*
        Wrap difference vector with respect to periodic boundary conditions.

        Parameters
        ----------
        fromPos :
            Pointer to initial point position.
        toPos :
            Pointer to final point position.
        unit :
            Return unitary vector.

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

        long int getHalfEdgeIndex(
            long int const& fromVertexIndex, long int const& toVertexIndex)
            const;
        /*
        Index of the half-edge between two vertices.
        Throws error if this half-edge does not exist.

        Parameters
        ----------
        fromVertexIndex :
            Index of origin vertex.
        toVertexIndex :
            Index of destination vertex.

        Returns
        -------
        halfEdgeIndex :
            Index of half-edge.
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

        TopoChangeEdgeInfoType deleteEdge(
            long int const& halfEdgeIndex);
        /*
        Delete edge by merging two vertices.

        Parameters
        ----------
        halfEdgeIntex :
            Index of half-edge linking two vertices to be merged.

        Returns
        -------
        deletedVertexIndex :
            Index of deleted vertex.
        deletedHalfEdgeIndices :
            Indices of half-edges deleted by this operation.
        */

        TopoChangeEdgeInfoType createEdge(
            long int const& halfEdgeIndex0, long int const& halfEdgeIndex1,
            double const& angle, double const& length=1,
            std::string type0="", std::string type1="");
        /*
        Create edge by splitting one vertex.

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
        length :
            Length to set for the new junction.
        type0 :
            Name of the type of half-edge for the half-edge going from the
            splitted vertex to the created vertex.
        type1 :
            Name of the type of half-edge for the half-edge going from the
            created vertex to the splitted vertex.
        NOTE: Other created half-edges will inherit the type of their pair
              half-edge, and the created vertex inheretit the type of the
              splitted vertex.

        Returns
        -------
        createdVertexIndex :
            Index of created vertex.
        createdHalfEdgeIndices :
            Indices of half-edges created by this operation.
        */

        void checkMesh() const;
        /*
        Check that the vertices and half-edges define a planar mesh, with
        anticlockwise triangles.
        */

};

void cerrTopoChangeEdgeInfo(TopoChangeEdgeInfoType const& info);
/*
Standard error output of topological changes edge info.
*/

#endif

