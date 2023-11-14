"""
Objects from simulation and plotting.
"""

from cells.bind import VertexModel as VM

import numpy as np
import math

import matplotlib.pyplot as plt

# SIMULATION OBJECT

class VertexModel(VM):
    """
    Physical system environment.
    """

    def __init__(self):
        """
        """

        super().__init__()
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
                color, lw = 'green', 3

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

        plt.pause(0.001)
        self.fig.canvas.draw()

