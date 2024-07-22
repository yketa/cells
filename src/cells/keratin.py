"""
Routine to run and plot in real time a simulation of the keratin vertex model.
This does not save data.
"""

from cells.init import init_vm, K, A0
from cells.plot import plot, _measure_fig, _resize_fig, _update_canvas,\
    _cbar_labelpad
from cells.bind import getLinesHalfEdge, getPolygonsCell, getLinesJunction
from cells.run import run

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np

# keratin colourbar
cmap_keratin = plt.cm.Greens                                    # colourmap
norm_keratin = Normalize(0, 25)                                 # interval of value represented by colourmap
scalarMap_keratin = ScalarMappable(norm_keratin, cmap_keratin)  # conversion from scalar value to colour

def plot_keratin(vm, time0=0, fig=None, ax=None, update=True, **kwargs):

    assert("keratin" in vm.vertexForces)
    param = vm.vertexForces["keratin"].parameters
    set_call = (fig is None or ax is None)
    fig, ax = plot(vm, fig=fig, ax=ax, only_set=True)   # initialise figure and axis

    # colourbars

    if set_call:

        ax_size, fig_width, fig_height = _measure_fig(ax)       # measure

        cbar_keratin = fig.colorbar(
            mappable=scalarMap_keratin, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_keratin.set_label(r"$[\mathrm{ker}]_i $",
            rotation=270, labelpad=_cbar_labelpad)
        cbar_keratin.ax.fill_between(                           # show keratin threshold...
            cbar_keratin.ax.get_xlim(), 0, param["kth"], color="white")
        cbar_keratin.ax.axhline(y=param["kth"], color="red")    # ... on the colourbar
        ax_size, fig_width, fig_height = (                      # resize
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # cells

    cells = vm.getVertexIndicesByType("centre")
    polygons = PatchCollection(
        list(map(                                               # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none")
    polygons.set_color(                                         # colour according to keratin
        (lambda keratin: list(map(
            lambda i:
                (lambda ki:
                    scalarMap_keratin.to_rgba(
                        ki, alpha=1 if ki > param["kth"] else 0))(
                keratin[i]),
            cells)))(
        vm.vertexForces["keratin"].keratin))
    ax.add_collection(polygons)

    # junctions

    lines = LineCollection(getLinesJunction(vm), colors="pink", # all junctions
        linewidths=2.5/max(1, np.sqrt(len(vm.vertices))/12))    # scale junction width with linear system size
    ax.add_collection(lines)

    # title

    title = (r"$t=%.3f, N_{\mathrm{T}_1}=%.3e, N_{\mathrm{cells}}=%i$"
        % (vm.time - time0, vm.nT1, len(cells)))
    title += r"$, \tau_r=%.2f, p_0=%.2f, T=%.2f, \tilde{p}_0=%.2f$" % (
        param["taur"], param["p0"], param["T"],
        param["p0"] - param["T"]/(2*param["Gamma"]*np.sqrt(param["A0"])))
    title += "\n"
    title += r"$\alpha=%.1e, \beta=%.1e$" % (
        param["alpha"], param["beta"])
    title += r"$, [\mathrm{ker}]_{\mathrm{th.}}=%.1e$" % param["kth"]
    title += r"$, \tau=%.1e, \sigma=%.1e$" % (
        param["tau"], param["sigma"])
    title += r"$, \tau_{\mathrm{on}}=$" + (
        r"$\infty$" if param["ron"] == 0 else r"$%.1e$" % (1./param["ron"]))
    if "pull" in vm.vertexForces:
        title += r"$, F_{\mathrm{pull}}=%.1e$" % (
            vm.vertexForces["pull"].parameters["F"])

    ax.set_title(title)

    # update canvas
    if update: _update_canvas(fig)

    return fig, ax

if __name__ == "__main__":

    # INITIALISATION

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    # K is set by default in cells.init
    # A0 is set by default in cells.init
    # taur is defined by {taur}
    # Gamma is defined by {Gamma}
    # p0 is defined by {p0}
    parser.add_argument("-T", type=float, default=0.1,
        help="junction tension")
    parser.add_argument("-alpha", type=float, default=1,
        help="keratin on-rate pressure-dependence parameter")
    parser.add_argument("-beta", type=float, default=1,
        help="keratin to area elasticity constant parameter")
    parser.add_argument("-kth", type=float, default=0.1,
        help="keratin concentration threshold")
    parser.add_argument("-tau", type=float, default=1e-1,
        help="keratin concentration evolution time scale")
    # sigma is defined by {sigma}
    parser.add_argument("-ron", type=float, default=0,
        help="keratin concentration on-rate evolution time rate (= 1/tauon)")
    parser.add_argument("-fpull", type=float, default=2,
        help="outer vertices pulling force scale")

    args, vm = init_vm(parser=parser)
    time0 = vm.time

    # KERATIN

    for _ in vm.vertexForces: vm.removeVertexForce(_)
    for _ in vm.halfEdgeForces: vm.removeHalfEdgeForce(_)
    vm.addKeratinModel("keratin",
        K, A0, args.taur,
        args.Gamma, args.p0, args.T,
        args.alpha, args.beta, args.kth,
        args.tau, args.sigma, args.ron)
    vm.addPressureForce("pull",
        args.fpull, True)

    # RUN

    run(args, vm, plot_function=plot_keratin, time0=time0)

