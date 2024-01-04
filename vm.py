"""
Run simulation of vertex model. Define objects and functions to read data file
and plot vertex model object.
"""

from cells.bind import VertexModel
from cells.bind import getLinesHalfEdge, getLinesJunction, getPolygonsCell
from cells.exponents import float_to_letters

import sys
from math import ceil

from datetime import datetime
import atexit, signal

import argparse

import numpy as np
from collections import OrderedDict

import pickle

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection

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

            # progress bar (https://stackoverflow.com/a/3160819/7385044)
            toolbar_width = 40
            def progressbar(p):
                out = "[%s] (%2i%%)" % (
                    "="*ceil(p*toolbar_width)
                        + " "*(toolbar_width - ceil(p*toolbar_width)),
                    ceil(100*p))
                sys.stdout.write(out)
                sys.stdout.flush()
                sys.stdout.write("\b"*len(out))

            # check file consistency and build skip directory
            if not(check): return
            self.skip = np.array([], dtype=int) # position of each object in file
            _max_diff_t = 0                     # check consistency with metadata in computed times
            for i, time in enumerate(self.frames*self.dt):
                progressbar(i/self.frames.size) # display progress bar
                self.skip = np.append(self.skip, current_pointer)
                vm = pickle.load(dump)          # load vertex model object
                assert(type(vm) == VertexModel) # check object has correct type
                if time == 0:                   # absolute difference in time
                    _max_diff_t = max(_max_diff_t, np.abs(time - vm.time))
                    assert(_max_diff_t == 0)
                else:                           # relative difference in time
                    _max_diff_t = max(_max_diff_t, np.abs(time - vm.time)/time)
                current_pointer = dump.tell()   # current position of the read pointer
            progressbar(1)
            sys.stdout.write("\n")              # end the progress bar
            try:
                pickle.load(dump)   # this should raise an EOFError if the file was read completely
                raise ValueError("File size is not consistent with metadata.")
            except EOFError:
                pass

        assert _max_diff_t < 1e-8   # relative difference between python and C++ times

        self.fig, self.ax = None, None  # used for plotting

    def plot(self, frame):
        """
        Plot vertex model state corresponding to frame.

        Parameters
        ----------
        frame : int
            Index of frame.
        """

        self.fig, self.ax = plot(self[frame], fig=self.fig, ax=self.ax)

    def play(self):
        """
        Plot all frames.
        """

        self.plot(self.frames[0])
        plt.ion()
        plt.show()
        for frame in self.frames[1:]:
            self.plot(frame)

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

def filename(N, v0, Dr, p0, identifier):
    """
    Standard filename for simulation file.

    Parameters
    ----------
    N : int
        Number of vertices.
    v0 : float
        Vertex self-propulsion velocity.
    Dr : float
        Vertex propulsion rotational diffusion constant.
    p0 : float
        Dimensionless target perimeter of cell.

    Returns
    -------
    name : str
        File name.
    """

    return ("N%s_v0%s_Dr%s_p0%s_%%%s.p"
        % tuple(map(float_to_letters, (N, v0, Dr, p0, identifier))))

# PLOT VERTEX MODEL OBJECT

cmap = plt.cm.jet                       # colourmap (colourmap without white for TESTING)
norm = Normalize(-0.05, 0.05)           # interval of value represented by colourmap
scalarMap = ScalarMappable(norm, cmap)  # conversion from scalar value to colour

def plot(vm, fig=None, ax=None):
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

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    # forces parameters

    hasForces = (
        "perimeter" in vm.vertexForces
        and "area" in vm.vertexForces
        and "abp" in vm.vertexForces)
    if hasForces:
        A0 = vm.vertexForces["area"].parameters["A0"]
        p0 = vm.vertexForces["perimeter"].parameters["P0"]/np.sqrt(A0)
        v0 = vm.vertexForces["abp"].parameters["v0"]
        Dr = 1./vm.vertexForces["abp"].parameters["taup"]

    # initialise figure

    if type(fig) == type(None) or type(ax) == type(None):
        fig, ax = plt.subplots()
        cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.05)
        colormap = mpl.colorbar.ColorbarBase(cax,
            cmap=cmap, norm=norm, orientation="vertical")

    # plot

    plt.sca(ax)
    ax.cla()
    ax.set_xlim([0, vm.systemSize[0]])
    ax.set_ylim([0, vm.systemSize[1]])
    ax.set_aspect("equal")

    #ax.plot(*getLinesHalfEdge(vm), color="blue", lw=1) # all half-edges
    ax.plot(*getLinesJunction(vm), color="red", lw=3)   # all junctions

    cells = [i for i in vm.vertices if vm.vertices[i].type == "centre"]
    areas = np.array(list(map(
        lambda i: vm.getVertexToNeighboursArea(i),
        cells)))

    polygons = PatchCollection(list(map(    # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))))
    if hasForces:
        polygons.set_color(list(map(        # colour according to area
            lambda area: scalarMap.to_rgba(area/A0 - 1),
            areas)))
        ax.add_collection(polygons)

    ax.set_title(
        r"$t=%.3f, N_{\mathrm{T}_1}=%.3e, N_{\mathrm{cells}}=%i$"
            % (vm.time, vm.nT1, len(cells))
        + ("" if not(hasForces) else
            r"$, v_0=%1.e, D_r=%.1e, p_0=%.2f$" % (v0, Dr, p0)))

    fig.canvas.draw_idle()
    fig.canvas.start_event_loop(0.001)

    return fig, ax

# SIMULATION

if __name__ == "__main__":

    # COMPUTATION TIME

    start_t = datetime.now()
    print("Started on %s." % start_t)

    def exit_handler():
        end_t = datetime.now()
        print("Stopped on %s (elapsed: %s)." % (end_t, end_t - start_t))
    # print elapsed time on exit
    signal.signal(signal.SIGINT, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    atexit.register(exit_handler)

    # PARAMETERS
    # `python -m cells.vertex_model -h` to display arguments and default values

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # physics
    parser.add_argument("-N", type=int, default=9,
        help="[open] total number of cells (! square number)")
    parser.add_argument("-n", type=int, default=6,
        help="[close] number of vertices in each direction (! multiple of 6)")
    parser.add_argument("-v0", type=float, default=1e-1,
        help="vertex self-propulsion velocity")
    parser.add_argument("-Dr", type=float, default=1e-1,
        help="vertex propulsion rotational diffusion constant")
    parser.add_argument("-p0", type=float, default=3.81,
        help="dimensionless target perimeter of cell")
    parser.add_argument("-open", "-o",
        action=argparse.BooleanOptionalAction,
        help="turn on specific checks for boundary vertices")
    # algorithm
    parser.add_argument("-seed", type=int, default=0,
        help="random number generator seed")
    parser.add_argument("-delta", type=float, default=0.1,
        help="length below which to perform T1")
    parser.add_argument("-epsilon", type=float, default=0.1,
        help="create junction with length epsilon above threshold after T1")
    # integration
    parser.add_argument("-dt", type=float, default=1e-2,
        help="intergation time step")
    parser.add_argument("-init", type=int, default=0,
        help="number of initial iterations")
    parser.add_argument("-niter", type=int, default=1000,
        help="number of production iterations")
    # saving
    parser.add_argument("-dtmin", type=int, default=1,
        help="lag time between each frame or [-log-frames] minimum lag time")
    parser.add_argument("-dtmax", type=int, default=500,
        help="[-log-frames] maximum lag time")
    parser.add_argument("-nmax", type=int, default=50,
        help="[-log-frames] (maximum) number of lag times")
    parser.add_argument("-intmax", type=int, default=20,
        help="[-log-frames] (maximum) number of initial times")
    parser.add_argument("-log-frames", "-log",
        action=argparse.BooleanOptionalAction,
        help="compute logarithmically spaced frames")
    parser.add_argument("-id", type=int, default=0,
        help="numerical identifier for simulation file")

    args = parser.parse_args()

    metadata = {        # metadata for simulation file
        "filename": filename(args.n**2, args.v0, args.Dr, args.p0, args.id),
        "dt": args.dt,  # used to check computed times
        "args": args,   # save all arguments
    }

    # CHOOSE FRAMES

    if not(args.log_frames):    # LINEARLY SPACED FRAMES
        if args.niter%args.dtmin:
            raise ValueError(
                "No integer multiple of %i frames in a total of %i frames."
                    % (args.dtmin, args.niter))
        else:
            metadata["t0"] = np.array(                                      # initial times
                [args.init                                                  # initialisation frames
                    + i*args.dtmin for i in range(args.niter//args.dtmin)], # frames spaced with dtmin
                dtype=int)
            metadata["t"] = np.array([0, args.dtmin])                       # lag times
    else:                       # LOGARITHMICALLY SPACED FRAMES
        dtmax = np.min([args.dtmax, args.niter - (args.intmax - 1)])        # maximum lag time
        metadata["t0"] = np.array(                                          # initial times
            [args.init                                                      # initialisation frames
                + i*(args.niter - dtmax)/(args.intmax - 1)
                    for i in range(args.intmax)],
            dtype=int)
        metadata["t"] = np.array(                                           # lag times
            [args.dtmin]                                                    # minimum lag time
                + [args.dtmin
                    + np.exp(i*np.log(dtmax - args.dtmin)/(args.nmax - 1))  # frames logarithmically spaced
                        for i in range(1, args.nmax - 1)]
                + [dtmax],                                                  # maximum lag time
            dtype=int)

    metadata["frames"] = np.array(list(OrderedDict().fromkeys(sorted(       # ensemble of unique frames
        [0, *[t0 + t
            for t0 in metadata["t0"]
            for t in [0, *metadata["t"]]]]))))

    # DUMP METADATA

    with open(metadata["filename"], "wb") as dump:
        pickle.dump(metadata, dump)
    print("Writing to \"%s\"." % metadata["filename"])

    # SIMULATION

    # initialisation of mesh
    vm = VertexModel(args.seed, args.v0, args.Dr, args.p0, args.open)
    if args.open:
        vm.initOpenRegularHexagonalLattice(nCells=args.N)
#         vm.initOpenRegularTriangularLattice(size=args.n)
    else:
        vm.initRegularTriangularLattice(size=args.n)

    # simulation of vertex model
    for t in np.diff(metadata["frames"], prepend=0):
        vm.nintegrate(t, metadata["dt"], args.delta, args.epsilon)
        vm.checkMesh()
        with open(metadata["filename"], "ab") as dump:
            pickle.dump(vm, dump)

