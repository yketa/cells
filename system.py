"""
Objects from simulation and plotting.
"""

from cells.tools import MultiIntKeyDict as MIKD
from cells.tools import Counter
from cells.mesh import Vertex, HalfEdge, Mesh

import numpy as np
import math

import matplotlib.pyplot as plt

# PHYSICS OBJECT

class SPVertex:
    """
    Physical self-propelled vertex.
    """

    def __init__(self, vertexIndex):
        """
        Parameters
        ----------
        vertexIndex : int
            Unique index for the cell identical to index of vertex.
        """

        self.vertexIndex = vertexIndex

        self.v0 = 1e-1							# self-propulsion velocity
        self.Dr = 1e-2							# rotational diffusion constant
        self.theta = 2*np.pi*np.random.random()	# direction of self-propulsion

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

        self.area = -1.			# area of cell
        self.kA	= 1.			# area stiffness
        self.A0 = None          # target area of cell (will be set in VertexModel.initRegularTriangularLattice)
        self.perimeter = -1.	# perimeter of cell
        self.kP = 1.			# perimeter stiffness
        self.p0 = 3.81	        # dimensionless target perimeter of cell

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
        self.color = 'green'
        ######

# SIMULATION OBJECT

class VertexModel(Mesh):
    """
    Physical system environment.
    """

    def __init__(self):
        """
        """

        super().__init__()
        self.reset()

        self.time = 0

    def reset(self):
        """
        """

        if hasattr(self, 'vertices'): del self.vertices
        self.vertices = dict()
        if hasattr(self, 'halfEdges'): del self.halfEdges
        self.halfEdges = dict()
        if hasattr(self, 'sPVertices'): del self.sPVertices
        self.sPVertices = MIKD()
        if hasattr(self, 'cells'): del self.cells
        self.cells = MIKD()
        if hasattr(self, 'faces'): del self.faces
        self.faces = MIKD()
        if hasattr(self, 'junctions'): del self.junctions
        self.junctions = MIKD()
        if hasattr(self, 'systemSize'): del self.systemSize

        self.time = 0

    def integrate(self, dt=0, delta=0.1, epsilon=0.1):
        """
        Integrate one time step of self-propelled vertex model and check for
        and perform T1s.

        Parameters
        ----------
        dt : float
            Integration time step. (default: 0)
        delta : float
            Distance between vertices below which these should be merge.
            (default: 1)
        epsilon : float
            Create two vertices at distance `delta' + `epsilon'. (default: 1)
        """

        # get forces

        forces = self.getForces()

        # integrate

        for vertexIndex in forces:
            self.vertices[vertexIndex].position += forces[vertexIndex]*dt	# Euler integration
            self.vertices[vertexIndex].position = self.wrap(				# wrap w.r.t. periodic boundary conditions
                self.vertices[vertexIndex].position)
        for vertexIndex in self.sPVertices:
            self.sPVertices[vertexIndex].theta += np.sqrt(
                2*self.sPVertices[vertexIndex].Dr*dt)*np.random.normal(0, 1)

        # move cell centres

        for vertex in self.cells:
            fromPos = self.vertices[vertex].position.copy()

            neighbours, _ = self.getNeighbours(vertex, junction=False)
            for neighbourVertexIndex in neighbours:
                toPos = self.vertices[neighbourVertexIndex].position

                disp = self.wrapDiff(fromPos, toPos)
                self.vertices[vertex].position += disp/len(neighbours)

            self.vertices[vertex].position = self.wrap(
                self.vertices[vertex].position)

        # perform T1s

        self.doT1(delta=delta, epsilon=epsilon)

        # update time

        self.time += dt

    def getForces(self):
        """
        Compute forces on vertices from vertex model.

        Returns
        -------
        forces : {int: (2,) float numpy array}
            Forces at each vertex.
        """

        forces = {index: np.array([0, 0], dtype=float)	# nitialise forces at each vertex
            for index in self.vertices}

        for cell in self.cells.values():	# loop over cells

            cell.area = self.getVertexToNeighboursArea(				# update area
                cell.vertexIndex)
            cell.perimeter = self.getVertexToNeighboursPerimeter(	# update perimeter
                cell.vertexIndex)

            vertexIndices, _ = self.getNeighbourVertices(
                cell.vertexIndex, angularSort=True)
            numberVertices = len(vertexIndices)
            for i, vertexIndex in enumerate(vertexIndices):	# loop over vertices of the cell
                assert not(vertexIndex in self.cells)

                previousVertexIndex = vertexIndices[(i - 1)%numberVertices]
                nextVertexIndex = vertexIndices[(i + 1)%numberVertices]

                # area term

                areaTerm = (
                    -cell.kA*(cell.area - cell.A0)*(1./2.)*(
                        np.cross(
                            self.wrapTo(	# vector from cell to previous vertex
                                cell.vertexIndex, previousVertexIndex,
                                unit=False),
                            [0, 0, 1])
                        - np.cross(
                            self.wrapTo(	# vector from cell to next vertex
                                cell.vertexIndex, nextVertexIndex,
                                unit=False),
                            [0, 0, 1]))
                    )[:2]
                forces[vertexIndex] += areaTerm
                forces[cell.vertexIndex] -= areaTerm	# enforce Newton's third law with cell centres

                # perimeter term

                forces[vertexIndex] += (
                    -cell.kP*(cell.perimeter - cell.p0*np.sqrt(cell.A0))*(
                        self.wrapTo(	# unit vector from next vertex to vertex
                            nextVertexIndex, vertexIndex, unit=True)
                        + self.wrapTo(	# unit vector from previous vertex to vertex
                            previousVertexIndex, vertexIndex, unit=True)))

        for sPVertex in self.sPVertices.values():	# loop over self-propelled vertices

            forces[sPVertex.vertexIndex] += sPVertex.v0*np.array(	# self-propulsion force
                [np.cos(sPVertex.theta), np.sin(sPVertex.theta)])

        return forces

    def doT1(self, delta=1, epsilon=1):
        """
        Check all junctions and perform T1s.

        Parameters
        ----------
        delta : float
            Distance between vertices below which these should be merge.
            (default: 1)
        epsilon : float
            Create two vertices at distance `delta' + `epsilon'. (default: 1)
        """

        # identify small junctions

        halfEdgeIndices = []
        for junction in self.junctions.values():
            if self.getEdgeLength(junction.halfEdgeIndex) < delta:
                halfEdgeIndices += [junction.halfEdgeIndex]	# this way each junction appears through an unique half-edge

        # perform T1s

        np.random.shuffle(halfEdgeIndices)	# do T1s in a random order
        for mergeHalfEdgeIndex in halfEdgeIndices:
            print('merge:', mergeHalfEdgeIndex)

            # identify half-edge to split to create new junction

            fromMergeIndex = self.halfEdges[mergeHalfEdgeIndex].fromIndex	# (first) vertex to be merge into neighbour
            toMergeIndex = self.halfEdges[mergeHalfEdgeIndex].toIndex		# (second) vertex towards which neighbour is merged

            neighboursFromMerge, halfEdgesNeighboursFromMerge = (	# neighbours of first vertex and half-edges towards them
                self.getNeighbourVertices(fromMergeIndex))
            neighboursToMerge, halfEdgesNeighboursToMerge = (		# neighbours of second vertex and half-edges towards them
                self.getNeighbourVertices(toMergeIndex))

            createHalfEdgeIndex0 = np.random.choice(	# randomly pick a half-edge towards a cell centre neighbour of first vertex
                [halfEdgeIndex
                    for vertexIndex, halfEdgeIndex in
                        zip(neighboursFromMerge, halfEdgesNeighboursFromMerge)
                    if vertexIndex in self.cells
                    and not(vertexIndex in neighboursToMerge)])
            createHalfEdgeIndex1 = np.random.choice(	# randomly pick a half-edge towards a cell centre neighbour of second vertex
                [halfEdgeIndex
                    for vertexIndex, halfEdgeIndex in
                        zip(neighboursToMerge, halfEdgesNeighboursToMerge)
                    if vertexIndex in self.cells
                    and not(vertexIndex in neighboursFromMerge)])

            angle = np.pi/2. + (	# create new junction orthogonal to the deleted junction
                math.atan2(*self.getHalfEdgeVector(mergeHalfEdgeIndex)[::-1]))

            # merge vertices

            self.mergeVertices(mergeHalfEdgeIndex)

            # create new vertex

            self.createJunction(createHalfEdgeIndex0, createHalfEdgeIndex1,
                length=delta + epsilon, angle=angle)

        if len(halfEdgeIndices) > 0: self.checkMesh()

    def mergeVertices(self, halfEdgeIndex):
        """
        Merge two vertices.

        Parameters
        ----------
        halfEdgeIntex : int
            Index of half-edge linking two vertices to be merged.

        Returns
        -------
        fromMergeIndex : int
            Merged vertex index.
        """

        # identify vertices to merge

        fromMergeIndex = self.halfEdges[halfEdgeIndex].fromIndex	# vertex to be merge into neighbour
        toMergeIndex = self.halfEdges[halfEdgeIndex].toIndex		# vertex towards which neighbour is merged

        # relabel half-edges origins and destinations

        previousHalfEdgeIndex = self.halfEdges[halfEdgeIndex].previousIndex	# half-edge pointing to the vertex to be removed in the first face to be removed
        endPreviousHalfEdgeIndex = self.halfEdges[halfEdgeIndex].pairIndex	# half-edge pointing to the vertex to be removed in the second face to be removed

        while True:													# loop through half-edge construction starting from first face to be removed until second face to removed

            pairHalfEdgeIndex = (									# half-edge pointing from the vertex to be removed and to be relabelled
                self.halfEdges[previousHalfEdgeIndex].pairIndex)

            # half-edge coming from vertex to be removed
            assert (
                self.halfEdges[pairHalfEdgeIndex].fromIndex
                    == fromMergeIndex)
            self.halfEdges[pairHalfEdgeIndex].fromIndex = toMergeIndex

            # half-edge going to vertex to be removed
            assert (
                self.halfEdges[previousHalfEdgeIndex].toIndex
                    == fromMergeIndex)
            self.halfEdges[previousHalfEdgeIndex].toIndex = toMergeIndex

            if previousHalfEdgeIndex == endPreviousHalfEdgeIndex:	# all faces which the vertex to be removed belongs to have been explored and half-edges relabelled
                assert pairHalfEdgeIndex == halfEdgeIndex
                break

            previousHalfEdgeIndex = (								# half-edge pointing to the vertex to be removed and to be relabelled
                self.halfEdges[pairHalfEdgeIndex].previousIndex)

        # reassign identifying half-edges of points belonging to faces to be removed

        previousHalfEdgeIndex = self.halfEdges[halfEdgeIndex].previousIndex	# half-edge pointing to the vertex to be removed in the first face to be removed
        pairHalfEdgeIndex = self.halfEdges[previousHalfEdgeIndex].pairIndex # half-edge pointing from the vertex to be removed

        assert self.halfEdges[pairHalfEdgeIndex].fromIndex == toMergeIndex
        self.vertices[toMergeIndex].halfEdgeIndex = pairHalfEdgeIndex

        nextHalfEdgeIndex = self.halfEdges[halfEdgeIndex].nextIndex		# half-edge pointing from the merge vertex in the first face to be removed
        pairHalfEdgeIndex = self.halfEdges[nextHalfEdgeIndex].pairIndex # half-edge pointing to the merge vertex
        thirdVertexIndex = self.halfEdges[pairHalfEdgeIndex].fromIndex  # index of third vertex in the first face to be removed

        self.vertices[thirdVertexIndex].halfEdgeIndex = pairHalfEdgeIndex

        pairHalfEdgeIndex = self.halfEdges[halfEdgeIndex].pairIndex		# half-edge pointing to the vertex to be removed in the second face to be removed
        nextHalfEdgeIndex = self.halfEdges[pairHalfEdgeIndex].nextIndex # half-edge pointing from the vertex to be removed in the second face to be removed
        pairHalfEdgeIndex = self.halfEdges[nextHalfEdgeIndex].pairIndex # half-edge pointing to the vertex to be removed
        thirdVertexIndex = self.halfEdges[pairHalfEdgeIndex].fromIndex  # index of third vertex in the second face to be removed

        self.vertices[thirdVertexIndex].halfEdgeIndex = pairHalfEdgeIndex

        # relabel half-edges pairs and delete/merge junctions

        fromHalfEdgeIndex = self.halfEdges[halfEdgeIndex].previousIndex
        fromHalfEdgeIndex = self.halfEdges[fromHalfEdgeIndex].pairIndex

        toHalfEdgeIndex = self.halfEdges[halfEdgeIndex].nextIndex
        toHalfEdgeIndex = self.halfEdges[toHalfEdgeIndex].pairIndex

        self.halfEdges[fromHalfEdgeIndex].pairIndex = toHalfEdgeIndex
        self.halfEdges[toHalfEdgeIndex].pairIndex = fromHalfEdgeIndex

        fromHalfEdgeIndex = self.halfEdges[halfEdgeIndex].pairIndex
        fromHalfEdgeIndex = self.halfEdges[fromHalfEdgeIndex].previousIndex
        fromHalfEdgeIndex = self.halfEdges[fromHalfEdgeIndex].pairIndex

        toHalfEdgeIndex = self.halfEdges[halfEdgeIndex].pairIndex
        toHalfEdgeIndex = self.halfEdges[toHalfEdgeIndex].nextIndex
        toHalfEdgeIndex = self.halfEdges[toHalfEdgeIndex].pairIndex

        self.halfEdges[fromHalfEdgeIndex].pairIndex = toHalfEdgeIndex
        self.halfEdges[toHalfEdgeIndex].pairIndex = fromHalfEdgeIndex

        # move merge vertex

        self.vertices[toMergeIndex].position += (	# move vertex to half point between two merged vertices
            self.wrapTo(toMergeIndex, fromMergeIndex, unit=False)/2.)

        # delete faces, junction, half-edges, and vertex

        assert halfEdgeIndex == self.halfEdges[halfEdgeIndex].index	# index of half-edge from deleted vertex
        pairHalfEdgeIndex = self.halfEdges[halfEdgeIndex].pairIndex	# index of half-edge to deleted vertex

        # first face
        previousIndex = self.halfEdges[halfEdgeIndex].previousIndex
        nextIndex = self.halfEdges[halfEdgeIndex].nextIndex
        del (																# delete face
            self.faces[halfEdgeIndex],
            self.faces[previousIndex],
            self.faces[nextIndex])
        del self.halfEdges[previousIndex], self.halfEdges[nextIndex]		# delete all half-edges but the one from the erased junction
        print('del1', previousIndex, nextIndex)

        # second face
        previousIndex = self.halfEdges[pairHalfEdgeIndex].previousIndex
        nextIndex = self.halfEdges[pairHalfEdgeIndex].nextIndex
        del (																# deleta face
            self.faces[pairHalfEdgeIndex],
            self.faces[previousIndex],
            self.faces[nextIndex])
        del self.halfEdges[previousIndex], self.halfEdges[nextIndex]		# delete all half-edges but the one from the erased junction
        print('del2:', previousIndex, nextIndex)

        del (	# delete junction
            self.junctions[halfEdgeIndex],
            self.junctions[pairHalfEdgeIndex])

        del (	# delete half-edges in the erased junction
            self.halfEdges[halfEdgeIndex],
            self.halfEdges[pairHalfEdgeIndex])
        print('del3:', halfEdgeIndex, pairHalfEdgeIndex)

        del self.vertices[fromMergeIndex]	# delete vertex
        del self.sPVertices[fromMergeIndex]	# delete self-propelled vertex

        return fromMergeIndex

    def createJunction(self, halfEdgeIndex0, halfEdgeIndex1, length=1,
        angle=None):
        """
        Create a new vertex and junction.

        Parameters
        ----------
        halfEdgeIndex0 : int
            Index of first half-edge going out of a vertex, and from whose pair
            half-edge it will be separated after the introduction of a new
            junction.
        halfEdgeIndex1 : int
            Index of second half-edge going out of the same vertex, and from
            whose pair half-edge it will be separated after the introduction of
            a new junction.
        distance : float
            Length to set for the new junction. (default: 1)
        angle : float or None
            Angle of the new junction with respect to the horizontal axis.
            (default: None)
            NOTE: if angle == None then the new vertices are displaced from the
                  original vertex position in an average direction of
                  neighbours.

        Returns
        -------
        newVertexIndex : int
            New vertex index.
        """

        # junctions cannot be split

        assert not(halfEdgeIndex0 in self.junctions)
        assert not(halfEdgeIndex1 in self.junctions)

        # create new vertex

        vertexIndex = self.halfEdges[halfEdgeIndex0].fromIndex
        assert vertexIndex == self.halfEdges[halfEdgeIndex1].fromIndex	# check that both half-edges go out of the same vertex
        self.vertices[vertexIndex].halfEdgeIndex = halfEdgeIndex0		# re-assign identifying half-edge to `halfEdgeIndex0' which will remain attached to `vertexIndex'

        newVertexIndex = max(self.vertices) + 1
        self.vertices[newVertexIndex] = Vertex(
            newVertexIndex,							# vertexIndex
            self.vertices[vertexIndex].position,	# position
            halfEdgeIndex1)							# halfEdgeIndex
        self.sPVertices[newVertexIndex] = SPVertex(
            newVertexIndex)							# vertexIndex

        # relabel origins and destinations of half of the half-edges

        halfEdgeIndex = halfEdgeIndex1
        while halfEdgeIndex != halfEdgeIndex0:

            assert self.halfEdges[halfEdgeIndex].fromIndex == vertexIndex
            self.halfEdges[halfEdgeIndex].fromIndex = newVertexIndex

            previousHalfEdgeIndex = self.halfEdges[halfEdgeIndex].previousIndex
            assert self.halfEdges[previousHalfEdgeIndex].toIndex == vertexIndex
            self.halfEdges[previousHalfEdgeIndex].toIndex = newVertexIndex

            halfEdgeIndex = self.halfEdges[previousHalfEdgeIndex].pairIndex

        # create new half-edges, faces, and junctions

        c = Counter(max(self.halfEdges) + 1)

        newHalfEdgeIndex0, previousNewHalfEdgeIndex0, nextNewHalfEdgeIndex0 = (
            c(), c(), c())
        newHalfEdgeIndex1, previousNewHalfEdgeIndex1, nextNewHalfEdgeIndex1 = (
            c(), c(), c())

        # first face
        self.halfEdges[newHalfEdgeIndex0] = HalfEdge(
            newHalfEdgeIndex0,							# halfEdgeIndex
            self.halfEdges[halfEdgeIndex0].toIndex,		# fromIndex
            vertexIndex,								# toIndex
            previousNewHalfEdgeIndex0,					# previousIndex
            nextNewHalfEdgeIndex0,						# nextIndex
            halfEdgeIndex0)								# pairIndex
        self.halfEdges[previousNewHalfEdgeIndex0] = HalfEdge(
            previousNewHalfEdgeIndex0,					# halfEdgeIndex
            newVertexIndex,								# fromIndex
            self.halfEdges[halfEdgeIndex0].toIndex,		# toIndex
            nextNewHalfEdgeIndex0,						# previousIndex
            newHalfEdgeIndex0,							# nextIndex
            self.halfEdges[halfEdgeIndex0].pairIndex)	# pairIndex
        self.halfEdges[nextNewHalfEdgeIndex0] = HalfEdge(
            nextNewHalfEdgeIndex0,						# halfEdgeIndex
            vertexIndex,								# fromIndex
            newVertexIndex,								# toIndex
            newHalfEdgeIndex0,							# previousIndex
            previousNewHalfEdgeIndex0,					# nextIndex
            nextNewHalfEdgeIndex1)						# pairIndex
        self.faces[
            newHalfEdgeIndex0, previousNewHalfEdgeIndex0, nextNewHalfEdgeIndex0
            ] = Face(newHalfEdgeIndex0)					# create new face
        self.halfEdges[self.halfEdges[halfEdgeIndex0].pairIndex].pairIndex = (
            previousNewHalfEdgeIndex0)					# relabel pair of old pair of halfEdgeIndex0
        self.halfEdges[halfEdgeIndex0].pairIndex = (
            newHalfEdgeIndex0)							# relabel pair of halfEdgeIndex0

        # second face
        self.halfEdges[newHalfEdgeIndex1] = HalfEdge(
            newHalfEdgeIndex1,							# halfEdgeIndex
            self.halfEdges[halfEdgeIndex1].toIndex,		# fromIndex
            newVertexIndex,								# toIndex
            previousNewHalfEdgeIndex1,					# previousIndex
            nextNewHalfEdgeIndex1,						# nextIndex
            halfEdgeIndex1)								# pairIndex
        self.halfEdges[previousNewHalfEdgeIndex1] = HalfEdge(
            previousNewHalfEdgeIndex1,					# halfEdgeIndex
            vertexIndex,								# fromIndex
            self.halfEdges[halfEdgeIndex1].toIndex,		# toIndex
            nextNewHalfEdgeIndex1,						# previousIndex
            newHalfEdgeIndex1,							# nextIndex
            self.halfEdges[halfEdgeIndex1].pairIndex)	# pairIndex
        self.halfEdges[nextNewHalfEdgeIndex1] = HalfEdge(
            nextNewHalfEdgeIndex1,						# halfEdgeIndex
            newVertexIndex,								# fromIndex
            vertexIndex,								# toIndex
            newHalfEdgeIndex1,							# previousIndex
            previousNewHalfEdgeIndex1,					# nextIndex
            nextNewHalfEdgeIndex0)						# pairIndex
        self.faces[
            newHalfEdgeIndex1, previousNewHalfEdgeIndex1, nextNewHalfEdgeIndex1
            ] = Face(newHalfEdgeIndex1)					# create new face
        self.halfEdges[self.halfEdges[halfEdgeIndex1].pairIndex].pairIndex = (
            previousNewHalfEdgeIndex1)					# relabel pair of old pair of halfEdgeIndex1
        self.halfEdges[halfEdgeIndex1].pairIndex = (
            newHalfEdgeIndex1)							# relabel pair of halfEdgeIndex1

        # create junction

        self.junctions[nextNewHalfEdgeIndex0, nextNewHalfEdgeIndex1] = (
            Junction(nextNewHalfEdgeIndex0))	# halfEdgeIndex

        # move vertices apart in the direction of neighbours' barycentre


        if type(angle) == type(None):

            neighboursVertex, _ = self.getNeighbourVertices(vertexIndex)
            toNeighbours = np.array([0, 0], dtype=float)			# vector in the direction of neighbours' barycentre
            for neighbour in neighboursVertex:
                toNeighbours += self.wrapTo(
                    vertexIndex, neighbour, unit=False)
            toNeighbours /= np.sqrt((toNeighbours**2).sum())		# normalise vector

            neighboursNewVertex, _ = self.getNeighbourVertices(newVertexIndex)
            toNewNeighbours = np.array([0, 0], dtype=float)			# vector in the direction of neighbours' barycentre
            for neighbour in neighboursNewVertex:
                toNewNeighbours += self.wrapTo(
                    newVertexIndex, neighbour, unit=False)
            toNewNeighbours /= np.sqrt((toNewNeighbours**2).sum())	# normalise vector

            assert (toNeighbours != toNewNeighbours).all()	# vectors to neighbours' barycentres have to be different
            assert np.cross(								# check that the created face has anticlockwise orientation
                self.getHalfEdgeVector(newHalfEdgeIndex0),
                toNewNeighbours - toNeighbours) > 0

            diff = np.sqrt(((toNeighbours - toNewNeighbours)**2).sum())	# norm of the difference between unitary vectors in the directions of neighbours' barycentre
            self.vertices[vertexIndex].position += (
                toNeighbours*length/diff)
            self.vertices[newVertexIndex].position += (
                toNewNeighbours*length/diff)

        else:

            direction = lambda _: np.array([np.cos(_), np.sin(_)])	# unit vector of angle `_'

            if np.cross(
                self.getHalfEdgeVector(newHalfEdgeIndex0),
                direction(angle)) < 0:
                angle += np.pi	# enforces that the newly created face has anticlockwise orientation

            self.vertices[vertexIndex].position -= (
                direction(angle)*length/2.)
            self.vertices[newVertexIndex].position += (
                direction(angle)*length/2.)

        return newVertexIndex

    def getNeighbours(self, vertexIndex, junction=False):
        """
        Return indices of neighbouring vertices and junctions.

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

        neighbourVerticesIndices, halfEdgesToNeighboursIndices = (
            self.getNeighbourVertices(vertexIndex))

        neighboursIndices, junctions = [], []
        for vertexIndex, halfEdgeIndex in zip(
            neighbourVerticesIndices, halfEdgesToNeighboursIndices):
            # check if the half-edge is a junction
            if not(junction) or (halfEdgeIndex in self.junctions):
                neighboursIndices += [vertexIndex]
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

                # create cell or self-propelled vertex
                if (line - column)%3 == 0:	                            # condition for vertex to be a cell centre
                    self.cells[vertexIndex] = Cell(vertexIndex)
                    self.cells[vertexIndex].A0 = (junctionLength**2)*(  # target area is the area of the regular hexagon
                        (3./2.)/np.tan(2*np.pi/12))
                else:						                            # then it is a self-propelled vertex
                    self.sPVertices[vertexIndex] = SPVertex(vertexIndex)

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

        self.checkMesh()	# check mesh construction

        # for plots
        if not(hasattr(self, 'fig')) or not(hasattr(self, 'ax')):
            self.fig, self.ax = plt.subplots()

    def plot(self):
        """
        Update plot of system.
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

