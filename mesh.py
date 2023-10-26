"""
Objects to define the mesh with half-edge procedure.
"""

from cells.dict import MultiIntKeyDict as MIKD

import numpy as np

import matplotlib.pyplot as plt

# STRUCTURE OBJECTS

class Vertex:
	"""
	Individual nodes of the mesh.
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
		self.position = position
		self.halfEdgeIndex = halfEdgeIndex

class HalfEdge:
	"""
	Directed arrow between vertices (Vertex) of the mesh.
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

# PHYSICS OBJECT

class Cell:
	"""
	Physical cell enclosed by junctions (HalfEdge) between vertices (Vertex) of
	the mesh.
	"""

	def __init__(self, vertexIndex):
		"""
		Parameters
		----------
		vertexIndex : int
			Unique index for the cell identical to index of vertex.
		"""

		self.vertexIndex = vertexIndex

class Face:
	"""
	Triangle enclosed by three vertices.
	"""

	def __init__(self, halfEdgeIndex):
		"""
		Parameters
		----------
		halfEdgeIndex : int
			Unique index of a helf-edge belonging to the face.
		"""

		self.halfEdgeIndex = halfEdgeIndex

class Junction:
	"""
	Physical link between two vertices.
	"""

	def __init__(self, halfEdgeIndex):
		"""
		Parameters
		----------
		halfEdgeIndex : int
			Unique index of a helf-edge belonging to the junction.
		"""

		self.halfEdgeIndex = halfEdgeIndex

		######
		# TEST
		self.w = np.random.randint(1, 4)
		self.phi = np.random.uniform(0, 2*np.pi)
		self.color = 'green'
		######

# SIMULATION OBJECT

class Mesh:
	"""
	Physical system environment.
	"""

	def __init__(self):
		"""
		"""

		self.reset()

	def reset(self):
		"""
		"""

		if hasattr(self, 'vertices'): del self.vertices
		self.vertices = dict()
		if hasattr(self, 'halfEdges'): del self.halfEdges
		self.halfEdges = dict()
		if hasattr(self, 'cells'): del  self.cells
		self.cells = MIKD()
		if hasattr(self, 'faces'): del self.faces
		self.faces = MIKD()
		if hasattr(self, 'junctions'): del self.junctions
		self.junctions = MIKD()
		if hasattr(self, 'systemSize'): del self.systemSize

		######
		# TEST
		self.time = 0
		######

	def getNeighbourVertices(self, vertexIndex, junction=False):
		"""
		Return indices of neighbouring vertices.

		Parameters
		----------
		vertexIndex : int
			Index of the vertex.
		junction : bool
			Restrict to vertices which are linked by a junction.
			(default: False)

		Returns
		-------
		neighboursIndices : (*,) int numpy array
			Indices of vertices neighbouring this vertex.
		junctions : (**,) list of cells.Junction
			List of junctions linked to this vertex.
		"""

		neighboursIndices = []
		junctions = []

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

			# check if the half-edge is a junction
			if not(junction) or halfEdgeIndex in self.junctions:
				neighboursIndices += [toVertexIndex]
			if halfEdgeIndex in self.junctions:
				junctions += [self.junctions[halfEdgeIndex]]

		neighboursIndices = np.array(neighboursIndices, dtype=int)
		return neighboursIndices, junctions

	def initRegularTriangularLattice(self, size=6, junctionLength=1):
		"""
		Initialises a regular triangular lattice.

		Parameters
		----------
		size : int
			Number of vertices in both horizontal and vertical directions.
			(default: 6)
			NOTE: Must be a multiple of 6.
		junctionLength : float
			Distance between nearest neighbour junctions. (default: 1)
		"""

		assert size%6 == 0

		self.reset()		# reset vertices, junctions, and cells
		self.systemSize = (	# length of the periodic box in each direction
			size*junctionLength, size*junctionLength*np.sqrt(3)/2)

		def getIndex(line, column):
			column = (column + (line//size)*(size//2))%size	# the system is periodic along an inclined vertical axis
			line = line%size								# the system is periodic along the horizontal axis
			return line*size + column

		# loop over vertices
		for line in range(size):
			for column in range(size):
				vertexIndex = line*size + column

				# create vertex
				self.vertices[vertexIndex] = Vertex(
					vertexIndex,
					junctionLength*np.array(						# position of the vertex on the regular triangular lattice
						[0.25 + line*0.5 + column, (0.5 + line)*np.sqrt(3)/2]),
					6*vertexIndex + 0)								# index of a single half-edge going out of this vertex
				self.vertices[vertexIndex].position = self.wrap(	# wrap coordinate around periodic boundary condition
					self.vertices[vertexIndex].position)

				# create cell
				if (line - column)%3 == 0:	# condition for vertex to be a cell centre
					self.cells[vertexIndex] = Cell(vertexIndex)

				# vertex indices for half-edge construction
				A = getIndex(line, column)
				assert A == vertexIndex
				B = getIndex(line, column + 1)
				C = getIndex(line + 1, column)
				D = getIndex(line + 1, column - 1)
				E = getIndex(line, column - 1)
				F = getIndex(line - 1, column)
				G = getIndex(line - 1, column + 1)

				# create half-edges
				# (0) from vertex (A) to right (B)
				self.halfEdges[6*A + 0] = HalfEdge(
					6*A + 0,	# halfEdgeIndex
					A,			# fromIndex
					B,			# toIndex
					6*A + 4,	# previousIndex
					6*A + 2,	# nextIndex
					6*A + 1)	# pairIndex
				# (1) from right (B) to vertex (A)
				self.halfEdges[6*A + 1] = HalfEdge(
					6*A + 1,	# halfEdgeIndex
					B,			# fromIndex
					A,			# toIndex
					6*G + 5,	# previousIndex
					6*F + 3,	# nextIndex
					6*A + 0)	# pairIndex
				# (2) from right (B) to top (C)
				self.halfEdges[6*A + 2] = HalfEdge(
					6*A + 2,	# halfEdgeIndex
					B,			# fromIndex
					C,			# toIndex
					6*A + 0,	# previousIndex
					6*A + 4,	# nextIndex
					6*A + 3)	# pairIndex
				# (3) from top (C) to right (B)
				self.halfEdges[6*A + 3] = HalfEdge(
					6*A + 3,	# halfEdgeIndex
					C,			# fromIndex
					B,			# toIndex
					6*C + 1,	# previousIndex
					6*B + 5,	# nextIndex
					6*A + 2)	# pairIndex
				# (4) from top (C) to vertex (A)
				self.halfEdges[6*A + 4] = HalfEdge(
					6*A + 4,	# halfEdgeIndex
					C,			# fromIndex
					A,			# toIndex
					6*A + 2,	# previousIndex
					6*A + 0,	# nextIndex
					6*A + 5)	# pairIndex
				# (5) from vertex (A) to top (C)
				self.halfEdges[6*A + 5] = HalfEdge(
					6*A + 5,	# halfEdgeIndex
					A,			# fromIndex
					C,			# toIndex
					6*E + 3,	# previousIndex
					6*D + 1,	# nextIndex
					6*A + 4)	# pairIndex

				# create faces
				# (0) right above vertex (A)
				self.faces[6*A + 0, 6*A + 2, 6*A + 4] = Face(
					6*A + 0)	# halfEdgeIndex
				# (1) right below vertex (A)
				self.faces[6*A + 1, 6*F + 3, 6*G + 5] = Face(
					6*A + 1)	# halfEdgeIndex

				# create junctions between vertices which are not cell centres
				if ((A//size - A%size)%3 != 0) and ((B//size - B%size)%3 != 0):
					self.junctions[6*A + 0, 6*A + 1] = Junction(
						6*A + 0)	# halfEdgeIndex
				if ((B//size - B%size)%3 != 0) and ((C//size - C%size)%3 != 0):
					self.junctions[6*A + 2, 6*A + 3] = Junction(
						6*A + 2)	# halfEdgeIndex
				if ((C//size - C%size)%3 != 0) and ((A//size - A%size)%3 != 0):
					self.junctions[6*A + 4, 6*A + 5] = Junction(
						6*A + 4)	# halfEdgeIndex

		# for plots
		self.fig, self.ax = plt.subplots()

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
		Wrap difference vector.

		Parameters
		----------
		fromPos : (2,) float array-like
			Initial point.
		toPos : (2,) float array-like
			Final point.

		Parameters
		----------
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

	def plot(self):
		"""
		"""

		assert hasattr(self, 'fig') and hasattr(self, 'ax')

		plt.sca(self.ax)
		self.ax.cla()
		self.ax.set_xlim([0, self.systemSize[0]])
		self.ax.set_ylim([0, self.systemSize[1]])
		self.ax.set_aspect('equal')

		for vertexIndex in self.vertices:
			color = ['blue', 'red'][vertexIndex in self.cells]
			self.ax.scatter(*self.wrap(self.vertices[vertexIndex].position),
				color=color, s=100)

		for halfEdgeIndex in self.halfEdges:

			color, lw = 'black', 1
			if halfEdgeIndex in self.junctions:
				color, lw = self.junctions[halfEdgeIndex].color, 3

			# indices of the vertices linked by the half-edge
			fromVertexIndex = self.halfEdges[halfEdgeIndex].fromIndex
			toVertexIndex = self.halfEdges[halfEdgeIndex].toIndex
			# positions of the vertices
			fromPos = self.wrap(self.vertices[fromVertexIndex].position)
			toPos = self.wrap(self.vertices[toVertexIndex].position)
			disp = self.wrapDiff(fromPos, toPos)

			self.ax.plot(
				[fromPos[0], fromPos[0] + disp[0]],
				[fromPos[1], fromPos[1] + disp[1]],
				color=color, lw=lw)

