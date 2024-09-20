"""
Routine to run and plot in real time a simulation of the keratin vertex model.
This does not save data.
"""

from cells.init import init_vm, K, A0, get_areas
from cells.plot import plot, _measure_fig, _resize_fig, _update_canvas,\
    _cbar_labelpad, norm_tension
from cells.bind import getLinesHalfEdge, getPolygonsCell, getLinesJunction
from cells.run import run

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import copy
import pickle

# tension colourbar
scalarMap_tension = ScalarMappable(norm_tension, plt.cm.bwr)    # conversion from scalar value to colour

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

        if "kmax" in kwargs:                                    # colourmap from user-defined max
            global scalarMap_keratin
            norm_keratin = Normalize(0, kwargs["kmax"])
            scalarMap_keratin = ScalarMappable(norm_keratin, cmap_keratin)

        # keratin colourbar
        cbar_keratin = fig.colorbar(
            mappable=scalarMap_keratin, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_keratin.set_label(r"$[\mathrm{ker}]_i $",
            rotation=270, labelpad=_cbar_labelpad)
#         cbar_keratin.ax.fill_between(                           # ignore colourbar values below threshold
#             cbar_keratin.ax.get_xlim(), 0, param["kth"], color="white")
        cbar_keratin.ax.axhline(y=param["kth"], color="red")    # display threshold on colourbar
        ax_size, fig_width, fig_height = (                      # resize
            _resize_fig(ax, ax_size, fig_width, fig_height))

        # tension colourbar
        cbar_tension = fig.colorbar(
            mappable=scalarMap_tension, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_tension.set_label(
            r"$(t_i - \left<t_i\right>)/\mathrm{std}(t_i)$",
            rotation=270, labelpad=_cbar_labelpad)
        ax_size, fig_width, fig_height = (
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # displayed quantities

    keratin = vm.vertexForces["keratin"].keratin    # keratin concentration
    tension = (                                     # scaled cell tension
        lambda t: (
            lambda m, s: {i: (t[i] - m)/s for i in t})(
            np.mean(list(t.values())), np.std(list(t.values()))))(
        vm.vertexForces["keratin"].tension)

    # cells

    cells = vm.getVertexIndicesByType("centre")
    polygons = PatchCollection(
        list(map(                                               # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none")
    polygons.set_linewidth(                                     # outer line width
        2.5/max(1, np.sqrt(len(vm.vertices))/12))
    polygons.set_edgecolor(                                     # edge colour according to keratin
        list(map(
            lambda i: scalarMap_tension.to_rgba(tension[i]),
            cells)))
    polygons.set_facecolor(                                     # face colour according to keratin
        list(map(
            lambda i: (
                lambda ki: scalarMap_keratin.to_rgba(ki))(
                keratin[i]),
#                 keratin[i] if keratin[i] > param["kth"] else 0),
            cells)))
    ax.add_collection(polygons)

#     # junctions
#
#     lines = LineCollection(getLinesJunction(vm), colors="pink", # all junctions
#         linewidths=2.5/max(1, np.sqrt(len(vm.vertices))/12))    # scale junction width with linear system size
#     ax.add_collection(lines)

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
    parser.add_argument("-xi", type=float, default=1,
        help="vertex drag coefficient")
    parser.add_argument("-K", type=float, default=K,
        help="area elasticity")
    parser.add_argument("-A0", type=float, default=A0,
        help="minimum target area")
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
    parser.add_argument("-keffmax", type=float, default=np.inf,
        help="maximum keratin effect")
    parser.add_argument("-tau", type=float, default=1e-1,
        help="keratin concentration evolution time scale")
    # sigma is defined by {sigma}
    parser.add_argument("-ron", type=float, default=0,
        help="keratin concentration on-rate evolution time rate (= 1/tauon)")
    parser.add_argument("-fpull", type=float, default=2,
        help="outer vertices pulling force scale")

    args, vm = init_vm(parser=parser, boxLength=5)
    time0 = vm.time

    # FORCES

    # remove forces
    for _ in vm.vertexForces: vm.removeVertexForce(_)
    for _ in vm.halfEdgeForces: vm.removeHalfEdgeForce(_)

    # rescale distances such that mean force is low (dichotomic search)
    vm0 = copy.deepcopy(vm)
    scalemin, scalemax = 0, 100
    while scalemax - scalemin > 1e-6:

        vm = copy.deepcopy(vm0)

        # --- set forces
        vm.setOverdampedIntegrator(
            args.xi)
        vm.addKeratinModel("keratin",
            args.K, args.A0, args.tau,
            args.Gamma, args.p0, args.T,
            args.alpha, args.beta,
            args.kth, agrs.keffmax,
            args.tau, args.sigma, args.ron)

        # --- compute forces
        scale = (scalemin + scalemax)/2
        vm.scale(scale*np.sqrt(args.A0/get_areas(vm).mean()))   # scale distances
        vm.nintegrate(1, 0)                                     # integrate with dt=0 to get forces
        meanForce = 0
        vertices = vm.getVertexIndicesByType("vertex")
        positions = (
            lambda p: np.array(list(map(
                lambda i: p[i],
                vertices))))(
            vm.getPositions())
        posCM = positions.mean(axis=0)
        forces = vm.forces
        for i, index in enumerate(vertices):
            meanForce += np.dot(positions[i] - posCM, forces[index]
                )/len(vertices)
        if np.abs(meanForce) < 1e-2*args.fpull: break           # initial radial force is small enough
        if meanForce > 0: scalemin = scale
        if meanForce < 0: scalemax = scale

    # pullin force
    vm.addPressureForce("pull",
        args.fpull, True)

    # RUN

    with open("keratin.init_vm.p", "wb") as dump: pickle.dump(vm, dump)
    run(args, vm, plot_function=plot_keratin, time0=time0)

