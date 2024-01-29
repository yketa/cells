"""
Define objects and functions to read data file and plot vertex model object.
When executed this checks the consistency of a simulation file.
"""

from cells.bind import VertexModel
from cells.bind import getLinesHalfEdge, getLinesJunction, getPolygonsCell

import sys

from math import ceil
import numpy as np

import pickle

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

# TOOLS

# progress bar (https://stackoverflow.com/a/3160819/7385044)
_toolbar_width = 40
def _progressbar(p):
    out = "[%s] (%2i%%)" % (
        "="*ceil(p*_toolbar_width)
            + " "*(_toolbar_width - ceil(p*_toolbar_width)),
        ceil(100*p))
    sys.stderr.write(out)
    sys.stderr.flush()
    sys.stderr.write("\b"*len(out))
    if p == 1: sys.stderr.write("\n")   # end the progress bar

# READ SIMULATION

class Read:

    def __init__(self, fname, check=True):
        """
        Extract metadata and check consistency of file.

        Parameters
        ----------
        fname : str
            Path to input file.
        check : bool
            Check file consistency. This is mandatory in order to access data,
            otherwise only metadata is read and accessible. (default: True)
        """

        self.filename = fname

        with open(self.filename, "rb") as dump:
            self.metadata = pickle.load(dump)
            assert type(self.metadata) == dict
            current_pointer = dump.tell()   # current position of the read pointer

            # default metadata
            self.t0 = self.metadata["t0"]           # initial times
            self.t = self.metadata["t"]             # lag times
            self.frames = self.metadata["frames"]   # array of computed frames
            self.dt = self.metadata["dt"]           # integration time step

            # check file consistency and build skip directory
            if not(check): return
            self.skip = np.array([], dtype=int)     # position of each object in file
            _max_diff_t = 0                         # check consistency with metadata in computed times
            for i, time in enumerate(self.frames*self.dt):
                _progressbar(i/self.frames.size)    # display progress bar
                self.skip = np.append(self.skip, current_pointer)
                vm = pickle.load(dump)              # load vertex model object
                assert(type(vm) == VertexModel)     # check object has correct type
                if time == 0:                       # absolute difference in time
                    _max_diff_t = max(_max_diff_t, np.abs(time - vm.time))
                    assert(_max_diff_t == 0)
                else:                               # relative difference in time
                    _max_diff_t = max(_max_diff_t, np.abs(time - vm.time)/time)
                current_pointer = dump.tell()       # current position of the read pointer
            _progressbar(1)
            try:
                pickle.load(dump)   # this should raise an EOFError if the file was read completely
                raise ValueError("File size is not consistent with metadata.")
            except EOFError:
                pass

        assert _max_diff_t < 1e-8   # relative difference between python and C++ times

        self.fig, self.ax = None, None  # used for plotting

    def plot(self, frame, rainbow=None):
        """
        Plot vertex model state corresponding to frame.

        Parameters
        ----------
        frame : int
            Index of frame.
        rainbow : int or None
            Index of previous frame of the system with respect to which colour
            cells. (default: None)
            NOTE: if rainbow != None then this overrides all other cell
                  colouring.
        """

        self.fig, self.ax = plot(self[frame], fig=self.fig, ax=self.ax,
            rainbow=rainbow if rainbow is None else self[rainbow])

    def play(self, frames=None, rainbow=None):
        """
        Plot frames.

        Parameters
        ----------
        frames : int array-like or None
            Frames to plot.
            NOTE: if frames == None then all frames in self.frames are plotted.
        rainbow : int or None
            Index of previous frame of the system with respect to which colour
            cells. (default: None)
            NOTE: if rainbow != None then this overrides all other cell
                  colouring.
        """

        if frames is None:
            frames = self.frames
        frames = list(frames)

        _progressbar(0)
        self.plot(frames[0])
        plt.ion()
        plt.show()
        for frame in frames[1:]:
            _progressbar(frames.index(frame)/len(frames))
            self.plot(frame, rainbow=rainbow)
        _progressbar(1)

    def __getitem__(self, frame):
        """
        Returns saved vertex model state corresponding to frame.

        Parameters
        ----------
        frame : int
            Index of frame.

        Returns
        -------
        vm : cells.bind.VertexModel
            Vertex model state.
        """

        assert frame in self
        with open(self.filename, "rb") as dump:
            dump.seek(self.skip[self.frames.tolist().index(frame)]) # skip until object (https://stackoverflow.com/questions/76252112)
            vm = pickle.load(dump)
            assert type(vm) == VertexModel

        return vm

    def __contains__(self, frame):  # is frame in saved frames?
        return self.frames.__contains__(frame)

    def __iter__(self):             # iterator over frames
        return self.frames.__iter__()

    def __next__(self):             # iterator over frames
        return self.frames.__next__()

    def __len__(self):              # number of frames
        return self.frames.__len__()

# PLOT VERTEX MODEL OBJECT

# area colourbar
cmap_area = plt.cm.bwr                                          # colourmap
norm_area = Normalize(-0.05, 0.05)                              # interval of value represented by colourmap
scalarMap_area = ScalarMappable(norm_area, cmap_area)           # conversion from scalar value to colour
# tension colourbar
cmap_tension = plt.cm.Spectral                                  # colourmap
norm_tension = Normalize(-1.5, 1.5)                             # interval of value represented by colourmap
scalarMap_tension = ScalarMappable(norm_tension, cmap_tension)  # conversion from scalar value to colour

def plot(vm, fig=None, ax=None, rainbow=None):
    """
    Plot vertex model.

    Parameters
    ----------
    vm : cells.bind.VertexModel
        State of the system to plot.
    fig : matplotlib.figure.Figure or None
        Figure on which to plot. (default: None)
        NOTE: if fig == None then a new figure and axes subplot is created.
    ax : matplotlib.axes._subplots.AxesSubplot or None
        Axes subplot on which to plot. (default: None)
        NOTE: if ax == None then a new figure and axes subplot is created.
    rainbow : cells.bind.VertexModel or None
        Previous state of the system with respect to which colour cells.
        (default: None)
        NOTE: if rainbow != None then this overrides all other cell colouring.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    # forces parameters

    if "area" in vm.vertexForces:
        A0 = vm.vertexForces["area"].parameters["A0"]
    if "perimeter" in vm.vertexForces:
        if "area" in vm.vertexForces:
            p0 = vm.vertexForces["perimeter"].parameters["P0"]/np.sqrt(A0)
        else:
            P0 = vm.vertexForces["perimeter"].parameters["P0"]
    if "abp" in vm.vertexForces:
        v0 = vm.vertexForces["abp"].parameters["v0"]
        taup = vm.vertexForces["abp"].parameters["taup"]
    if "out" in vm.halfEdgeForces:
        t0 = vm.halfEdgeForces["out"].parameters["t0"]
        taup = vm.halfEdgeForces["out"].parameters["taup"]
    for model in ("model0", "model1"):
        if model in vm.halfEdgeForces:
            if "area" in vm.vertexForces:
                p0 = vm.halfEdgeForces[model].parameters["P0"]/np.sqrt(A0)
            else:
                P0 = vm.halfEdgeForces[model].parameters["P0"]
    for model in ("model0", "model1", "model2", "model3", "model4"):
        if model in vm.halfEdgeForces:
            sT0 = vm.halfEdgeForces[model].parameters["sigma"]
            taup = vm.halfEdgeForces[model].parameters["taup"]
    for model in ("model2", "model4"):
        if model in vm.halfEdgeForces:
            taur = vm.halfEdgeForces[model].parameters["taur"]

    # initialise figure

    if fig is None or ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)                                         # set figure size
        try: fig.canvas.window().setFixedSize(fig.canvas.window().size())   # set window size
        except AttributeError: pass
        if "area" in vm.vertexForces:
            cbar_area = plt.colorbar(
                mappable=scalarMap_area, ax=ax, shrink=0.5)
            cbar_area.set_label(r"$A_i/A_0 - 1$", rotation=270)
        if "out" in vm.halfEdgeForces:
            cbar_tension = plt.colorbar(
                mappable=scalarMap_tension, ax=ax, shrink=0.5)
            cbar_tension.set_label(r"$t_i/t_0 - 1$", rotation=270)
        for model in ("model0", "model1", "model2", "model3", "model4"):
            if model in vm.halfEdgeForces:
                cbar_tension = plt.colorbar(
                    mappable=scalarMap_tension, ax=ax, shrink=0.5)
                cbar_tension.set_label(r"$t_i/\sigma$", rotation=270)

    # plot

    plt.sca(ax)
    ax.cla()
    ax.set_xlim([0, vm.systemSize[0]])
    ax.set_ylim([0, vm.systemSize[1]])
    ax.set_aspect("equal")

    # junctions and half-edges
    lines = LineCollection(getLinesJunction(vm), colors="red", linewidths=1.5)  # all junctions
    if "t0" in locals() or "sT0" in locals():
        junctions = [i for i in sorted(vm.halfEdges)
            if vm.halfEdges[i].type == "junction"]
        if "t0" in locals():
            tensions = np.concatenate(list(map(
                lambda i: [vm.halfEdgeForces["out"].tension[i]]*2,
                junctions)))
            lines.set_color(list(map(
                lambda tension: scalarMap_tension.to_rgba(tension/t0 - 1),
                tensions)))
        elif "sT0" in locals():
            for model in ("model0", "model1", "model2", "model3", "model4"):
                try:
                    tensions = np.concatenate(list(map(
                        lambda i: [vm.halfEdgeForces[model].tension[i]]*2,
                        junctions)))
                except:
                    continue
            lines.set_color(list(map(
                lambda tension: scalarMap_tension.to_rgba(tension/sT0),
                tensions)))
    ax.add_collection(lines)
    #ax.plot(*getLinesHalfEdge(vm), color="blue", lw=1)                          # all half-edges

    # cells
    polygons = PatchCollection(
        list(map(                           # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none")
    cells = [i for i in sorted(vm.vertices)
        if vm.vertices[i].type == "centre"]
    if type(rainbow) is VertexModel:
        scalarMap_rainbow = ScalarMappable(
            Normalize(0, rainbow.systemSize[0]), plt.cm.hsv)
        positions0 = np.array(list(map(
            lambda i: rainbow.vertices[i].position,
            cells)))
        polygons.set_color(list(map(        # colour according to previous position
            lambda position0: scalarMap_rainbow.to_rgba(position0[0]),
            positions0)))
    elif "A0" in locals():
        areas = np.array(list(map(
            lambda i: vm.getVertexToNeighboursArea(i),
            cells)))
        polygons.set_color(list(map(        # colour according to area
            lambda area: scalarMap_area.to_rgba(area/A0 - 1),
            areas)))
    ax.add_collection(polygons)

    title = (r"$t=%.3f, N_{\mathrm{T}_1}=%.3e, N_{\mathrm{cells}}=%i$"
        % (vm.time, vm.nT1, len(cells)))
    if "p0" in locals():
        title += r"$, p_0=%.2f$" % p0
    if "P0" in locals():
        title += r"$, P_0=%.2f$" % P0
    if "v0" in locals():
        title += r"$, v_0=%.1e, \tau_p=%.1e$" % (v0, taup)
    if "t0" in locals():
        title += r"$, t_0=%1.e, \tau_p=%.1e$" % (t0, taup)
    if "sT0" in locals():
        title += r"$, \sigma=%.1e, \tau_p=%.1e$" % (sT0, taup)
    if "taur" in locals():
        title += r"$, \tau_r=%.1e$" % taur
    ax.set_title(title)

    fig.canvas.draw_idle()
    fig.canvas.start_event_loop(0.001)

    return fig, ax

# CHECK

if __name__ == "__main__":

    if len(sys.argv) == 2:
        try:
            Read(sys.argv[1], check=True)
            print("true")   # print lowercase to be treated as bash boolean
        except:
            print("false")  # print lowercase to be treated as bash boolean

