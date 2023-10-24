"""
Objects to define the mesh with half-edge procedure.
"""

class Vertex:
	"""
	Individual nodes of the mesh.
	"""

	def __init__(self, index, position,
		halfEdgeIndex=-1, cellIndex=-1):
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
		cellIndex : int
			Cell index if this vertex is the centre of a cell. (default: -1)
		"""

		self.index = index
		self.position = position
		self.halfEdgeIndex = halfEdgeIndex
		self.cellIndex = cellIndex

	def getNeighbouringVertices(self):
		"""
		Returns list of indices of neighbouring vertices.
		"""

		return

class Cell:
	"""
	Physical cell enclosed by junctions (HalfEdge) between nodes (Vertex) of
	the mesh.
	"""

	def __init__(self, index):
		"""
		Parameters
		----------
		index : int
			Unique index for the cell identical to index of vertex.
		"""

		self.index = index

class HalfEdge:
	"""
	Directed junctions between nodes (Vertex) of the mesh.
	"""

	def __init__(self, index, fromIndex, toIndex,
		previousIndex=-1, nextIndex=-1, pairIndex=-1, junction=False):
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
		junction : bool
			Is this half-edge associated to a physical junction between
			vertices? (default: False)
		"""

		self.index = index
		self.fromIndex = fromIndex
		self.toIndex = toIndex
		self.previousIndex = previousIndex
		self.nextIndex = nextIndex
		self.pairIndex = pairIndex
		self.junction = junction

class Face:
	"""
	Triangle enclosed by three vertices.
	"""

	def __init__(self, index, halfEdgeIndex):
		"""
		Parameters
		----------
		index : int
			Unique index for the face.
		halfEdgeIndex : int
			Index of a helf-edge belonging to the face.
		"""

		self.index = index
		self.halfEdgeIndex = halfEdgeIndex

class Mesh:
	"""
	Ensemble of vertices/nodes (Vertex), junctions (HalfEdge), and cells (Cell)
	constituting the system.
	"""

	def __init__(self):
		"""
		"""

		self.vertices = {}
		self.halfedges = {}
		self.cells = {}

	def reset(self):
		"""
		"""

		del self.vertices
		self.vertices = {}
		del self.halfEdges
		self.halfEdges = {}
		del self.cells
		self.cells = {}

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

		self.reset()	# reset vertices, junctions, and cells

		# loop over vertices
		for line in range(size):
			for column in range(size):
				vertexIndex = line*size + column

				# create vertex
				self.vertices[vertexIndex] = Vertex(
					vertexIndex,
					junctionLength*np.array(	# position of the vertex on the regular triangular lattice
						[line*0.5 + column, line*np.sqrt(3)/2]))

				# create cell
				if (line - column)%3 == 0:	# condition for vertex to be a cell centre
					self.cells[vertexIndex] = Cell(vertexIndex)

				# points for half-edge construction
				A = line*size + column
				assert A == vertexIndex
				B = line*size + (column + 1)%size
				C = ((line + 1)%size)*size + column
				D = ((line + 1)%size)*size + (column - 1)%size
				E = line*size + (column - 1)%size
				F = ((line - 1)%size)*size + column
				G = ((line - 1)%size)*size + (column + 1)%size

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
				self.faces[2*A + 0] = Face(
					2*A + 0,	# faceIndex
					6*A + 0)	# halfEdgeIndex
				# (1) right below vertex (A)
				self.faces[2*A + 1] = Face(
					2*A + 1,	# faceIndex
					6*A + 1)	# halfEdgeIndex
