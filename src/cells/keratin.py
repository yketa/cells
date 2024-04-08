"""
Routine to run and plot in real time a simulation of the keratin vertex model.
This does not save data.
"""

from cells.init import init_vm, K, A0
from cells.plot import plot, _measure_fig, _resize_fig, _update_canvas
from cells.bind import getLinesHalfEdge, getPolygonsCell
from cells.run import run

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

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

def plot_keratin(vm, fig=None, ax=None, **kwargs):

    assert("keratin" in vm.vertexForces)
    set_call = (fig is None or ax is None)
    fig, ax = plot(vm, fig=fig, ax=ax, only_set=True)   # initialise figure and axis

#     print(vm.vertexForces["keratin"].pressure)

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

    halfEdges = list(
        set(list(vm.vertexForces["keratin"].tension.keys())).intersection(
            set(list(vm.halfEdges))))
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
    title += r"$, \tau_r=%.2f, p_0=%.2f, l_0=%.2f, \alpha=%.1e$" % (
        param["taur"], param["p0"], param["l0"], param["alpha"])
    title += r"$, [\mathrm{ker}]_{\mathrm{th.}}=%.1e$" % param["kth"]
    title += "\n"
    title += r"$\tau=%.1e, \sigma=%.1e$" % (
        param["tau"], param["sigma"])
    title += r"$, \tau_{\mathrm{on}}=$" + (
        r"$\infty$" if param["ron"] == 0 else r"$%.1e$" % (1./param["ron"]))
    title += r"$, k_0=%.1e, T_0=%.1e$" % (
        param["k0"], param["pr0"])
    if "pull" in vm.vertexForces:
        title += r"$, F_{\mathrm{pull}}=%.1e$" % (
            vm.vertexForces["pull"].parameters["Fpull"])

    ax.set_title(title)

    # update canvas
    _update_canvas(fig)

    return fig, ax

if __name__ == "__main__":

    # INITIALISATION

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    # K is set by default in cells.init
    # A0 is set by default in cells.init
    # taur is defined by {taur}
    # Gamma is defined by {Gamma}
    # p0 is defined by {p0}
    parser.add_argument("-l0", type=float, default=1,
        help="bond rest length")
    parser.add_argument("-alpha", type=float, default=0.2,
        help="bond elasticity per keratin concentration above threshold")
    parser.add_argument("-kth", type=float, default=0,
        help="bond keratin concentration threshold")
    parser.add_argument("-tau", type=float, default=10,
        help="keratin concentration evolution time scale")
    # sigma is defined by {sigma}
    parser.add_argument("-ron", type=float, default=0,
        help="keratin concentration on-rate evolution time rate (= 1/tauon)")
    parser.add_argument("-k0", type=float, default=0.5,
        help="keratin concentration off-rate inverse pressure constant")
    parser.add_argument("-pr0", type=float, default=-0.3,
        help="keratin concentration off-rate inflection pressure")
    parser.add_argument("-fpull", type=float, default=1,
        help="outer vertices pulling force")

    args, vm = init_vm(parser=parser)

    # KERATIN

    for _ in vm.vertexForces: vm.removeVertexForce(_)
    for _ in vm.halfEdgeForces: vm.removeHalfEdgeForce(_)
    vm.addKeratinModel("keratin",
        K, args.taur, args.Gamma, args.p0, args.l0, args.alpha, args.kth,
        args.tau, args.sigma, args.ron, args.k0, args.pr0)
    vm.addEdgePullForce("pull",
        args.fpull)

    # RUN

    run(args, vm, plot_function=plot_keratin)

