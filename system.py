"""
Objects from simulation and plotting.
"""

from cells.bind import VertexModel as VM
from cells.bind import getLinesHalfEdge, getLinesJunction

import numpy as np
import math

import matplotlib.pyplot as plt

# SIMULATION OBJECT

class VertexModel(VM):
    """
    Physical system environment.
    """

    def __init__(self, seed=0, v0=1e-1, Dr=1e-1, p0=3.81):
        """
        Parameters
        ----------
        seed : int
            Random number generator seed. (default: 0)
        v0 : float
            Vertex self-propulsion velocity. (default: 1e-1)
        Dr : float
            Vertex propulsion rotational diffusion constant. (default: 1e-1)
        p0 : float
            Dimensionless target perimeter of cell. (default: 3.81)
        """

        super().__init__(seed, v0, Dr, p0)
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

        self.ax.plot(*getLinesHalfEdge(self), color='blue', lw=1)   # all half-edges
        self.ax.plot(*getLinesJunction(self), color='red', lw=3)    # all junctions

        self.fig.canvas.draw_idle()
        self.fig.canvas.start_event_loop(0.001)

