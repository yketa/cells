"""
Objects from simulation and plotting.
"""

from cells.bind import VertexModel
from cells.bind import getLinesHalfEdge, getLinesJunction, getPolygonsCell

import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection

# SIMULATION OBJECT

class ModelSystem(VertexModel):
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

    def plot(self):
        """
        Update plot of system.
        """

        # initialise figure
        if not(hasattr(self, 'fig')) or not(hasattr(self, 'ax')):
            self.fig, self.ax = plt.subplots()
            self.cmap = plt.cm.PiYG
            self.norm = Normalize(-0.05, 0.05)
            self.cax = make_axes_locatable(self.ax).append_axes(
                'right', size='5%', pad=0.05)
            self.colormap = mpl.colorbar.ColorbarBase(self.cax,
                cmap=self.cmap, norm=self.norm, orientation='vertical')
            self.scalarMap = ScalarMappable(self.norm, self.cmap)

        plt.sca(self.ax)
        self.ax.cla()
        self.ax.set_xlim([0, self.systemSize[0]])
        self.ax.set_ylim([0, self.systemSize[1]])
        self.ax.set_aspect('equal')

#         self.ax.plot(*getLinesHalfEdge(self), color='blue', lw=1)   # all half-edges
        self.ax.plot(*getLinesJunction(self), color='red', lw=3)    # all junctions

        polygons = PatchCollection(list(map(
                lambda vertices: plt.Polygon(vertices, closed=True),
                getPolygonsCell(self))))
        polygons.set_color(list(map(
            lambda i: self.scalarMap.to_rgba(
                self.cells[i].area/self.cells[i].A0 - 1),
            self.cells)))
        self.ax.add_collection(polygons)

        self.ax.set_title(
            r'$t=%.3f, N_{\mathrm{T}_1}=%.3e, N_{\mathrm{cells}}=%i,$'
                % (self.time, self.nT1, len(self.cells))
            + r'$v_0=%1.e, D_r=%.1e, p_0=%.2f$'
                % (self.v0, self.Dr, self.p0))

        self.fig.canvas.draw_idle()
        self.fig.canvas.start_event_loop(0.001)

