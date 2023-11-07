"""
Objects to define the mesh with half-edge procedure.
"""

from cells.tools import MultiIntKeyDict

import numpy as np
import math

class Vertex:
	"""
	Individual nodes of the two-dimensional mesh.
	"""

	def __init__(self, index, position,
		halfEdgeIndex=-1):
		"""
		Parameters
		----------
		index : int
			Unique index for the vertex.
		position : (2,) float array-like
			Position of the vertex.
		halfEdgeIndex: int
			Index of a single halfedge going out of this vertex (the
			association needs not to be bijective). (default: -1)
		"""

		self.index = index
		self.position = np.copy(np.array(position))
		self.halfEdgeIndex = halfEdgeIndex

class HalfEdge:
	"""
	Directed arrow between vertices (Vertex) of the two-dimensional mesh.
	"""

	def __init__(self, index, fromIndex, toIndex,
		previousIndex=-1, nextIndex=-1, pairIndex=-1):
		"""
		Parameters
		----------
		index : int
			Unique index for the half-edge.
		fromIndex : int
			Index of the origin vertex.
		toIndex : int
			Index of the destination vertex.
		previousIndex : int
			Index of the associated `previous' half-edge going towards the
			origin of this half-edge. (default: -1)
		nextIndex : int
			Index of the associated `next' half-edge going outwards from the
			destination of this half-edge. (default: -1)
		pairIndex : int
			Index of the associated `pair' half-edge going from the destination
			to the origin of this half-edge. (default: -1)
		"""

		self.index = index
		self.fromIndex = fromIndex
		self.toIndex = toIndex
		self.previousIndex = previousIndex
		self.nextIndex = nextIndex
		self.pairIndex = pairIndex

class Mesh:
	"""
	Two-dimensional ensembles of vertices and edges.
	"""

	def __init__(self):
		"""
		"""

		self.vertices = dict()
		self.helfEdges = dict()
		self.systemSize = np.array([0, 0])

	def wrap(self, position):
		"""
		Wrap position with respect to periodic boundary conditions.

		Parameters
		----------
		position : (2,) float array-like
			Unwrapped position.

		Returns
		-------
		wrappedPosition : (2,) float numpy array
			Wrapped position.
		"""

		return np.array(
			(position[0]%self.systemSize[0], position[1]%self.systemSize[1]),
			dtype=float)

	def wrapDiff(self, fromPos, toPos):
		"""
		Wrap difference vector with respect to periodic boundary conditions.

		Parameters
		----------
		fromPos : (2,) float array-like
			Initial point.
		toPos : (2,) float array-like
			Final point.

		Returns
		-------
		disp : (2,) float numpy array
			Difference vector.
		"""

		disp = np.array(toPos, dtype=float) - np.array(fromPos, dtype=float)
		assert disp.ndim == 1 and disp.size == 2
		for dim in range(2):
			if np.abs(disp[dim]) >= self.systemSize[dim]/2:
				disp[dim] = -np.sign(disp[dim])*(
					self.systemSize[dim] - np.abs(disp[dim]))

		return disp

	def wrapTo(self, fromVertexIndex, toVertexIndex, unit=False):
		"""
		Vector between two vertices, accounting for periodic boundary
		conditions.

		Parameters
		----------
		fromVertexIndex : int
			Index of origin vertex.
		toVertexIndex : int
			Index of destination vertex.
		unit : bool
			Return unitary vector. (default: False)

		Returns
		-------
		fromTo : (2,) float numpy array
			Vector between vertices.
		"""

		fromTo = self.wrapDiff(
			self.vertices[fromVertexIndex].position,
			self.vertices[toVertexIndex].position)

		if unit: fromTo /= (fromTo**2).sum()
		return fromTo

	def getHalfEdgeVector(self, halfEdgeIndex, unit=False):
		"""
		Vector going from the origin to the destination of a half-edge.

		Parameters
		----------
		halfEdgeIndex : int
			Index of the half-edge.
		unit : bool
			Return unitary vector. (default: False)

		Returns
		-------
		halfEdgeVector : (2,) float numpy array
			Vector corresponding to half-edge.
		"""

		fromVertexIndex = self.halfEdges[halfEdgeIndex].fromIndex
		toVertexIndex = self.halfEdges[halfEdgeIndex].toIndex

		return self.wrapTo(fromVertexIndex, toVertexIndex, unit=unit)

	def getEdgeLength(self, halfEdgeIndex):
		"""
		Length of edge.

		Parameters
		----------
		halfEdgeIndex : int
			Index of the half-edge.

		Returns
		-------
		edgeLength : float
			Length of edge.
		"""

		return np.sqrt((self.getHalfEdgeVector(halfEdgeIndex)**2).sum())

	def getNeighbourVertices(self, vertexIndex, angularSort=False):
		"""
		Indices of neighbouring vertices and half-edges towards them.

		Parameters
		----------
		vertexIndex : int
			Index of the vertex.
		angularSort : bool
			Sort neighbours in anticlockwise order. (default: False)

		Returns
		-------
		neighbourVerticesIndices : (*,) int numpy array
			Indices of vertices neighbouring this vertex.
		halfEdgesToNeighboursIndices : (*,) int numpy array
			Indices of half-edges from this vertex towards neighbour vertices.
		"""

		neighbourVerticesIndices = []
		halfEdgesToNeighboursIndices = []

		# find destination vertex in half-edge construction
		halfEdgeIndex = self.vertices[vertexIndex].halfEdgeIndex
		assert halfEdgeIndex >= 0										# check that the half-edge exists
		assert self.halfEdges[halfEdgeIndex].fromIndex == vertexIndex	# check that the half-edge goes out of this vertex
		firstNeighbourVertexIndex = self.halfEdges[halfEdgeIndex].toIndex

		# loop around neighbours
		toVertexIndex = -1
		while toVertexIndex != firstNeighbourVertexIndex:

			# get next half-edge (next - next - pair)
			nextHalfEdgeIndex = self.halfEdges[halfEdgeIndex].nextIndex
			nextHalfEdgeIndex = self.halfEdges[nextHalfEdgeIndex].nextIndex
			halfEdgeIndex = self.halfEdges[nextHalfEdgeIndex].pairIndex

			# find destination vertex in half-edge construction
			assert halfEdgeIndex >= 0										# check that the half-edge exists
			assert self.halfEdges[halfEdgeIndex].fromIndex == vertexIndex	# check that the half-edge goes out of this vertex
			toVertexIndex = self.halfEdges[halfEdgeIndex].toIndex

			neighbourVerticesIndices += [toVertexIndex]
			halfEdgesToNeighboursIndices += [halfEdgeIndex]

		if angularSort:
			edgeAngles = MultiIntKeyDict()
			for vertexIndex, halfEdgeIndex in zip(
				neighbourVerticesIndices, halfEdgesToNeighboursIndices):
				edgeAngles[vertexIndex, halfEdgeIndex] = (
					math.atan2(
						*self.getHalfEdgeVector(halfEdgeIndex, unit=True)[::-1]
						))
			neighbourVerticesIndices = sorted(
				neighbourVerticesIndices, key=lambda i: edgeAngles[i])			# sorted indices in anticlockwise order
			halfEdgesToNeighboursIndices = sorted(
				halfEdgesToNeighboursIndices, key=lambda i: edgeAngles[i])		# sorted indices in anticlockwise order

		neighbourVerticesIndices = np.array(neighbourVerticesIndices,
			dtype=int)
		halfEdgesToNeighboursIndices = np.array(halfEdgesToNeighboursIndices,
			dtype=int)
		return neighbourVerticesIndices, halfEdgesToNeighboursIndices

	def getVertexToNeighboursArea(self, vertexIndex):
		"""
		Area encapsulated by the neighbours of a vertex (shoelace formula).

		Parameters
		----------
		vertexIndex : int
			Index of vertex.

		Returns
		-------
		area : float
			Area encapsulated by neighbours.
		"""

		_, halfEdgeIndices = self.getNeighbourVertices(vertexIndex,
			angularSort=True)	# neighbours in anticlockwise order
		numberNeighbours = halfEdgeIndices.size
		assert numberNeighbours >= 3

		# compute area
		area = 0
		for i in range(numberNeighbours):
			area += np.cross(	# area of triangle
				self.getHalfEdgeVector(
					halfEdgeIndices[i%numberNeighbours]),
				self.getHalfEdgeVector(
					halfEdgeIndices[(i + 1)%numberNeighbours])
				)/2.

		return area

	def getVertexToNeighboursPerimeter(self, vertexIndex):
		"""
		Perimeter encapsulated by the neighbours of a vertex.

		Parameters
		----------
		vertexIndex : int
			Index of vertex.

		Returns
		-------
		perimeter : float
			Perimeter encapsulated by neighbours.
		"""

		vertexIndices, _ = self.getNeighbourVertices(vertexIndex,
			angularSort=True)	# neighbours in anticlockwise order
		numberNeighbours = vertexIndices.size
		assert numberNeighbours >= 3

		# compute perimeter
		perimeter = 0
		for i in range(numberNeighbours):
			perimeter += np.sqrt(
				(self.wrapTo(	# difference vector between neighbouring vertices
					vertexIndices[i%numberNeighbours],
					vertexIndices[(i + 1)%numberNeighbours],
					unit=False)**2).sum())

		return perimeter

	def checkMesh(self):
		"""
		Check that the vertices and half-edges define a planar mesh, with
		anticlockwise triangles.
		"""

		vertices = list(self.vertices)		# list of vertices
		halfEdges = list(self.halfEdges)	# list of half-edges

		for halfEdgeIndex in self.halfEdges:	# loop over all half-edges

			if not(halfEdgeIndex in halfEdges): continue

			nextHalfEdgeIndex = self.halfEdges[halfEdgeIndex].nextIndex
			previousHalfEdgeIndex = self.halfEdges[halfEdgeIndex].previousIndex

			triangle = (	# three half-edges forming a triangle
				halfEdgeIndex,
				self.halfEdges[halfEdgeIndex].nextIndex,
				self.halfEdges[halfEdgeIndex].previousIndex)

			for i, halfEdge in enumerate(triangle):

				fromVertex = self.halfEdges[halfEdge].fromIndex
				toVertex = self.halfEdges[halfEdge].toIndex

				if fromVertex in vertices:		# looping over all half-edges should encounter all vertices

					halfEdgeBis = self.vertices[fromVertex].halfEdgeIndex
					assert self.halfEdges[halfEdgeBis].fromIndex == fromVertex	# check consistency with identifying half-edge

					vertices.remove(fromVertex)

				assert np.cross(												# check that the triangle has anticlockwise orientation
					self.getHalfEdgeVector(triangle[(i + 0)%3]),
					self.getHalfEdgeVector(triangle[(i + 1)%3])) > 0

				pairHalfEdgeIndex = self.halfEdges[halfEdge].pairIndex
				assert halfEdge == (	# check the indexing of pair
					self.halfEdges[pairHalfEdgeIndex].pairIndex)
				assert toVertex == (
					self.halfEdges[pairHalfEdgeIndex].fromIndex)
				assert fromVertex == (
					self.halfEdges[pairHalfEdgeIndex].toIndex)

				assert halfEdge == (	# check the indexing of previous
					self.halfEdges[triangle[(i + 1)%3]].previousIndex)
				assert toVertex == (
					self.halfEdges[triangle[(i + 1)%3]].fromIndex)

				assert halfEdge == (	# check the indexing of next
					self.halfEdges[triangle[(i - 1)%3]].nextIndex)
				assert fromVertex == (
					self.halfEdges[triangle[(i - 1)%3]].toIndex)

				halfEdges.remove(halfEdge)

		assert len(vertices) == 0
		assert len(halfEdges) == 0
		print('OK mesh construction')

