"""
Routine to run and plot in real time a simulation of the keratin vertex model.
This does not save data.
"""

import sys
sys.argv[0] = "run.py"  # appear as "run.py"

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
cmap_keratin = plt.cm.PRGn                                        # colourmap
norm_keratin = (lambda x: Normalize(-x, x))(1.5)                                 # interval of value represented by colourmap
# norm_keratin = Normalize(-2, 2)                                 # interval of value represented by colourmap
scalarMap_keratin = ScalarMappable(norm_keratin, cmap_keratin)  # conversion from scalar value to colour

# tension colourbar
cmap_tension = plt.cm.bwr                                     # colourmap
norm_tension = (lambda x: Normalize(-x, x))(0.3)                                 # interval of value represented by colourmap
# norm_tension = Normalize(-2, 2)                                 # interval of value represented by colourmap
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
        cbar_keratin.set_label(
            r"$([\mathrm{ker}]_i - \left<[\mathrm{ker}]_i\right>)$"
                + r"$/\mathrm{std}([\mathrm{ker}]_i)$",
            rotation=270, labelpad=20)
        ax_size, fig_width, fig_height = (                  #resize
            _resize_fig(ax, ax_size, fig_width, fig_height))

        cbar_tension = fig.colorbar(
            mappable=scalarMap_tension, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_tension.set_label(
            r"$(t_i - \left<t_i\right>)/\mathrm{std}(t_i)$",
            rotation=270, labelpad=20)
        ax_size, fig_width, fig_height = (                  #resize
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # cells

    cells = vm.getVertexIndicesByType("centre")
    polygons = PatchCollection(
        list(map(                           # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none")
    ker_mean, ker_std = 0, 1
#     ker_mean = np.mean(list(vm.vertexForces["keratin"].keratin.values()))
#     ker_std = np.std(list(vm.vertexForces["keratin"].keratin.values()))
    if ker_std != 0:
        polygons.set_color(list(map(        # colour according to keratin
            lambda i: (lambda s_ker: scalarMap_keratin.to_rgba(s_ker))(
                (vm.vertexForces["keratin"].keratin[i] - ker_mean)/ker_std),
            cells)))
    ax.add_collection(polygons)

    # bonds

    halfEdges = list(vm.vertexForces["keratin"].tension.keys())
    if len(halfEdges) > 0:
        lines = LineCollection(
            getLinesHalfEdge(vm, halfEdges),    # all half-edges with tension
            linewidth=0.5)
        ten_mean, ten_std = 0, 1
#         ten_mean = np.mean(list(vm.vertexForces["keratin"].tension.values()))
#         ten_std = np.std(list(vm.vertexForces["keratin"].tension.values()))
        if ten_std != 0:
            lines.set_color(list(map(           # colour according to tension
                lambda i: (lambda s_ten: scalarMap_tension.to_rgba(s_ten))(
                    (vm.vertexForces["keratin"].tension[i] - ten_mean
                        )/ten_std),
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
    title += r"$\tau=%.1e, \sigma=%.1e, \tau_{\mathrm{on}}=%.1e$" % (
        param["tau"], param["sigma"], param["tauon"])
    title += r"$, k_0=%.1e, T_0=%.1e$" % (
        param["k0"], param["p0"])

#     for param, tex in (
#         ("K", "K"), ("A0", "A_0"),
#         ("Gamma", "\\Gamma"), ("P0", "P_0"),
#         ("l0", "l_0"), ("alpha", "\\alpha"),
#             ("kth", "[\\mathrm{ker}]_ {\\mathrm{th}.}"),
#         ("tau", "\\tau"), ("sigma", "\\sigma"),
#         ("tauon", "\\tau_{\\mathrm{on}}"), ("k0", "k_0"), ("p0", "p_0")):
#         title += r"$, %s=%.2e$" % (
#             tex, vm.vertexForces["keratin"].parameters[param])
#         count += 1
#         if count % 5 == 0:
#             title += "\n"
    ax.set_title(title)

    # update canvas
    _update_canvas(fig)

    return fig, ax

if __name__ == "__main__":

    def exit_handler(*_args, **_kwargs):
        # make movie on exit
        if "args" in globals() and args.movie:
            subprocess.call([os.path.join(__path__[0], "movie.sh"),
                "-d", tmpdir.name, "-p", sys.executable, # "-f", args.ffmpeg,
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
    parser.add_argument("-tauon", type=float, default=100,
        help="keratin concentration on-rate evolution time scale")
    parser.add_argument("-k0", type=float, default=1,
        help="keratin concentration off-rate inverse pressure constant")
    parser.add_argument("-pr0", type=float, default=-0.3,
        help="keratin concentration off-rate inflection pressure")

    args, vm = init_vm(parser)

    if args.movie: tmpdir = TemporaryDirectory()

    # KERATIN

    assert len(vm.vertexForces) == 0
    assert len(vm.halfEdgeForces) == 0
    vm.addKeratinModel("keratin",
        1, A0, 1, args.vmp0*np.sqrt(A0), args.l0, args.alpha, args.kth,
        args.tau, args.sigma, args.tauon, args.k0, args.pr0)

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
        #if vm.time > 10: exit()
        # plot
        plot_keratin(vm, fig=fig, ax=ax)

