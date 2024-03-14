"""
Routine to run and plot in real time a simulation of the keratin vertex model.
This does not save data.
"""

from cells.init import init_vm, A0
from cells.plot import plot, _measure_fig, _resize_fig, _update_canvas
from cells.bind import getLinesHalfEdge, getPolygonsCell
from cells import __path__

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

import subprocess
import os
import signal
from tempfile import TemporaryDirectory
import atexit
import traceback

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np

# keratin colourbar
cmap_keratin = plt.cm.PRGn                                      # colourmap
norm_keratin = Normalize(-1.5, 1.5)                             # interval of value represented by colourmap
scalarMap_keratin = ScalarMappable(norm_keratin, cmap_keratin)  # conversion from scalar value to colour

# tension colourbar
cmap_tension = plt.cm.bwr                                       # colourmap
norm_tension = Normalize(-0.3, 0.3)                             # interval of value represented by colourmap
scalarMap_tension = ScalarMappable(norm_tension, cmap_tension)  # conversion from scalar value to colour

def plot_keratin(vm, fig=None, ax=None):

    assert("keratin" in vm.vertexForces)
    set_call = (fig is None or ax is None)
    fig, ax = plot(vm, fig=fig, ax=ax, only_set=True)   # initialise figure and axis

    # colorbars

    if set_call:

        ax_size, fig_width, fig_height = _measure_fig(ax)   # measure

        cbar_keratin = fig.colorbar(
            mappable=scalarMap_keratin, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_keratin.set_label(r"$[\mathrm{ker}]_i $",
            rotation=270, labelpad=20)
        ax_size, fig_width, fig_height = (                  # resize
            _resize_fig(ax, ax_size, fig_width, fig_height))

        cbar_tension = fig.colorbar(
            mappable=scalarMap_tension, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_tension.set_label(r"$t_i$",
            rotation=270, labelpad=20)
        ax_size, fig_width, fig_height = (                  # resize
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # cells

    cells = vm.getVertexIndicesByType("centre")
    polygons = PatchCollection(
        list(map(                               # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none")
    polygons.set_color(list(map(                # colour according to keratin
        lambda i: scalarMap_keratin.to_rgba(
            vm.vertexForces["keratin"].keratin[i]),
        cells)))
    ax.add_collection(polygons)

    # bonds

    halfEdges = list(vm.vertexForces["keratin"].tension.keys())
    if True:    # remove inner cell half-edges (which go from or to a centre vertex)
        halfEdges = list(set(halfEdges) - set(vm.getCentreHalfEdges()))
    if len(halfEdges) > 0:
        lines = LineCollection(
            getLinesHalfEdge(vm, halfEdges),    # all half-edges with tension
            linewidth=1.5)
        lines.set_color(list(map(               # colour according to tension
            lambda i: scalarMap_tension.to_rgba(
                vm.vertexForces["keratin"].tension[i]),
            halfEdges)))
        ax.add_collection(lines)

    # title

    param = vm.vertexForces["keratin"].parameters
    title = (r"$t=%.3f, N_{\mathrm{T}_1}=%.3e, N_{\mathrm{cells}}=%i$"
        % (vm.time, vm.nT1, len(cells)))
    title += r"$, p_0=%.2f, l_0=%.2f, \alpha=%.1e$" % (
        param["P0"]/np.sqrt(param["A0"]), param["l0"], param["alpha"])
    title += r"$, [\mathrm{ker}]_{\mathrm{th.}}=%.1e$" % param["kth"]
    title += "\n"
    title += r"$\tau=%.1e, \sigma=%.1e$" % (
        param["tau"], param["sigma"])
    title += r"$, \tau_{\mathrm{on}}=$" + (
        r"$\infty$" if param["ron"] == 0 else r"$%.1e$" % (1./param["ron"]))
    title += r"$, k_0=%.1e, T_0=%.1e$" % (
        param["k0"], param["p0"])

    ax.set_title(title)

    # update canvas
    _update_canvas(fig)

    return fig, ax

if __name__ == "__main__":

    def exit_handler(*_args, **_kwargs):
        # make movie on exit
        if "args" in globals() and args.movie:
            subprocess.call([os.path.join(__path__[0], "movie.sh"),
                "-d", tmpdir.name, "-p", sys.executable, # "-F", args.ffmpeg,
                "-y"])
            tmpdir.cleanup()
        # exit
        os._exit(0)
    signal.signal(signal.SIGINT, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    atexit.register(exit_handler)

    # INITIALISATION

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    # K is 1 by default
    # A0 is (3./2.)/np.tan(np.pi/6.) by default
    # Gamma is 1 by default
    parser.add_argument("-vmp0", type=float, default=3.72,
        help="dimensionless target perimeter")
    parser.add_argument("-l0", type=float, default=1,
        help="bond rest length")
    parser.add_argument("-alpha", type=float, default=0.1,
        help="bond elasticity per keratin concentration above threshold")
    parser.add_argument("-kth", type=float, default=0,
        help="bond keratin concentration threshold")
    parser.add_argument("-tau", type=float, default=100,
        help="keratin concentration evolution time scale")
    # sigma is defined by {sigma}
    parser.add_argument("-ron", type=float, default=0,
        help="keratin concentration on-rate evolution time rate (= 1/tauon)")
    parser.add_argument("-k0", type=float, default=1,
        help="keratin concentration off-rate inverse pressure constant")
    parser.add_argument("-pr0", type=float, default=-0.3,
        help="keratin concentration off-rate inflection pressure")

    args, vm = init_vm(parser)

    if args.movie: tmpdir = TemporaryDirectory()

    # KERATIN

#     assert len(vm.vertexForces) == 0
#     assert len(vm.halfEdgeForces) == 0
    vm.addKeratinModel("keratin",
        1, A0, 1, args.vmp0*np.sqrt(A0), args.l0, args.alpha, args.kth,
        args.tau, args.sigma, args.ron, args.k0, args.pr0)

    fig, ax = plot_keratin(vm)

    # RUN

    plt.ion()
    plt.show()
    while True:
        # save frame
        if args.movie:
            try: count += 1
            except NameError: count = 0
            fig.savefig(os.path.join(tmpdir.name, "%05d.png" % count))
        # integrate
        try:
            vm.nintegrate(args.iterations,
                dt=args.dt, delta=args.delta, epsilon=args.epsilon)
#         vm.checkMesh(["junction"])
        except:
            print(traceback.format_exc(), file=sys.stderr)  # print traceback
            exit_handler()                                  # exit with handler
        # plot
        plot_keratin(vm, fig=fig, ax=ax)

