from cells.mesh import Mesh

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize as ColorsNormalise
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

m = Mesh()
m.initRegularTriangularLattice(size=6)

cax = make_axes_locatable(m.ax).append_axes('right', size='5%', pad=0.05)
cmap = plt.cm.PiYG
norm = ColorsNormalise(-1, 1)
colormap = mpl.colorbar.ColorbarBase(cax, cmap, norm, orientation='vertical')
scalarMap = ScalarMappable(norm, cmap)

def update(iterations=100, dt=5e-3, A=5e-2):

	for iteration in range(iterations):

		# compute forces

		forces = {vertex: np.array([0., 0.]) for vertex in m.vertices}

		for halfEdgeIndex in m.halfEdges:
			if halfEdgeIndex in m.junctions:	# loop over junctions

				fromIndex = m.halfEdges[halfEdgeIndex].fromIndex
				toIndex = m.halfEdges[halfEdgeIndex].toIndex
				disp = (
					m.vertices[toIndex].position
					- m.vertices[fromIndex].position)
				for dim in range(2):                                            
					if np.abs(disp[dim]) >= m.systemSize[dim]/2:             
						disp[dim] = -np.sign(disp[dim])*(                       
							m.systemSize[dim] - np.abs(disp[dim])) 
				disp /= np.sqrt((disp**2).sum())	# normalised vector out of vertex

				amp = np.cos(
					m.junctions[halfEdgeIndex].w*m.time
					+ m.junctions[halfEdgeIndex].phi)
				forces[fromIndex] += A*disp*amp
				m.junctions[halfEdgeIndex].color = scalarMap.to_rgba(amp)

		# integrate positions

		for vertex in m.vertices:

			m.vertices[vertex].position += forces[vertex]*dt
			m.vertices[vertex].position = m.wrap(m.vertices[vertex].position)

		# move cell centres

		for vertex in m.cells:
			fromPos = m.vertices[vertex].position.copy()

			neighbours, _ = m.getNeighbourVertices(vertex, junction=False)
			for neighbourVertexIndex in neighbours:
				toPos = m.vertices[neighbourVertexIndex].position

				disp = m.wrapDiff(fromPos, toPos)
				m.vertices[vertex].position += disp/len(neighbours)

			m.vertices[vertex].position = m.wrap(m.vertices[vertex].position)

		# update time

		m.time += dt

	# plot

	m.plot()

anim = animation.FuncAnimation(m.fig, update, repeat=True, interval=0)
plt.show()
