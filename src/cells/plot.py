"""
Define objects functions to plot vertex model object.
"""

from cells.bind import VertexModel, angle2
from cells.bind import getLinesHalfEdge, getLinesJunction, getPolygonsCell

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

# PLOT VERTEX MODEL OBJECT

def _update_canvas(fig):
    fig.canvas.draw_idle()
    fig.canvas.start_event_loop(0.001)
    assert(plt.fignum_exists(fig.number))

def plot(vm, fig=None, ax=None, update=True,
    rainbow=None, clear=False, only_set=False):
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
    update : bool
        Update figure canvas. (default: True)
    rainbow : cells.bind.VertexModel or None
        Previous state of the system with respect to which colour cells.
        (default: None)
        NOTE: if rainbow != None then this overrides all other cell colouring.
    clear : bool
        Clear the plot of all cell colouring. (default: False)
        NOTE: `rainbow' and `clear' options are exclusionary.
    only_set : bool
        Only set the figure and do not plot anything. (default: False)
        NOTE: `only_set' implies `clear'.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    rainbow_plot = (type(rainbow) == VertexModel)   # is it a rainbow plot?
    if only_set: clear = True
    if rainbow_plot and clear:
        raise ValueError("`rainbow' and `clear' options are exclusionary.")

    # forces parameters

    if not(clear):
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

        plt.ioff()
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)                                         # set figure size
        try: fig.canvas.window().setFixedSize(fig.canvas.window().size())   # set window size
        except AttributeError: pass

        # all force-related colorbars
        if not(rainbow_plot):

            ax_size, fig_width, fig_height = _measure_fig(ax)

            if "area" in vm.vertexForces:
                cbar_area = fig.colorbar(
                    mappable=scalarMap_area, ax=ax,
                    shrink=0.75, pad=0.01)
                cbar_area.set_label(
                    r"$(A_i - \left<A_i\right>)/\mathrm{std}(A_i)$",
                    rotation=270, labelpad=20)
                # resize
                ax_size, fig_width, fig_height = (
                    _resize_fig(ax, ax_size, fig_width, fig_height))

            if "out" in vm.halfEdgeForces:
                cbar_tension = fig.colorbar(
                    mappable=scalarMap_tension, ax=ax,
                    shrink=0.75, pad=0.01)
                cbar_tension.set_label(
                    r"$(t_i - \left<t_i\right>)/\mathrm{std}(t_i)$",
                    rotation=270, labelpad=20)
                # resize
                ax_size, fig_width, fig_height = (
                    _resize_fig(ax, ax_size, fig_width, fig_height))

            for model in ["model%i" % i for i in range(5)]:
                if model in vm.halfEdgeForces:
                    cbar_tension = fig.colorbar(
                        mappable=scalarMap_tension, ax=ax,
                        shrink=0.75, pad=0.01)
                    cbar_tension.set_label(
                        r"$(t_i - \left<t_i\right>)/\mathrm{std}(t_i)$",
                        rotation=270, labelpad=20)
                    # resize
                    ax_size, fig_width, fig_height = (
                        _resize_fig(ax, ax_size, fig_width, fig_height))

    # plot

    plt.sca(ax)
    ax.cla()
    ax.set_xlim([0, vm.systemSize[0]])
    ax.set_ylim([0, vm.systemSize[1]])
    ax.set_aspect("equal")
    if only_set: return fig, ax

    # junctions and half-edges
    lines = LineCollection(getLinesJunction(vm), colors="red", linewidths=2.5)  # all junctions
    if not(rainbow_plot):
        if ("t0" in locals() or "sT0" in locals()):
            junctions = vm.getHalfEdgeIndicesByType("junction")
            if "t0" in locals():
                tensions = np.concatenate(list(map(
                    lambda i: [vm.halfEdgeForces["out"].tension[i]]*2,
                    junctions)))
            elif "sT0" in locals():
                for model in ["model%i" % i for i in range(5)]:
                    try:
                        tensions = np.concatenate(list(map(
                            lambda i: [vm.halfEdgeForces[model].tension[i]]*2,
                            junctions)))
                    except:
                        continue
            tensions_std = tensions.std()
            if tensions_std != 0:
                lines.set_color(list(map(
                    lambda s_tension: scalarMap_tension.to_rgba(s_tension),
                    (tensions - tensions.mean())/tensions_std)))
    ax.add_collection(lines)
    #ax.plot(*getLinesHalfEdge(vm), color="blue", lw=1)                          # all half-edges

    # cells
    polygons = PatchCollection(
        list(map(                           # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none")
    cells = vm.getVertexIndicesByType("centre")
    if rainbow_plot:
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
        areas_std = areas.std()
        if areas_std != 0:
            polygons.set_color(list(map(    # colour according to area
                lambda s_area: scalarMap_area.to_rgba(s_area),
                (areas - areas.mean())/areas_std)))
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

    # update canvas
    if update: _update_canvas(fig)

    return fig, ax

def plot_forces(vm, fig=None, ax=None, zero=1e-10, av_norm=0.2, centres=False,
    **kwargs):
    """
    Plot vertex model with forces on vertices.

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
    zero : float
        Discard arrows with infinite norm below threshold. (default: 1e-10)
    av_norm : float
        Average norm of arrows. (default: 0.2)
    centres : bool
        Only display average forces of cell cornes on centre vertices.
        (default: False)

    Additional keywords arguments are passed to cells.plot.plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    fig, ax = plot(vm, fig=fig, ax=ax, update=False, **kwargs)

    positions = vm.getPositions(wrapped=True)

    if centres: forces = vm.getCentreForces()
    else: forces = vm.getForces()
    forces = (
        lambda f: {i: f[i]
            for i in f if np.abs(f[i]).max() > zero and i in positions})(
        forces)
    av_norm = (1./av_norm)*np.sqrt(
        (np.array(list(forces.values()), dtype=float)**2).sum(axis=-1).mean())

    if av_norm != 0:
        scalarMap = ScalarMappable(
            norm=Normalize(vmin=0, vmax=2*np.pi), cmap=plt.cm.hsv)
        arrows = PatchCollection(
            list(map(
                lambda i:
                    plt.arrow(*positions[i], *np.array(forces[i])/av_norm,
                        width=5e-2, head_width=4e-1,
                        length_includes_head=False,
                        color=scalarMap.to_rgba(angle2(*forces[i])%(2*np.pi)),
                        zorder=10),
                forces)))
        ax.add_collection(arrows)

    _update_canvas(fig)

    return fig, ax

# COLOURBARS

# area colourbar
cmap_area = plt.cm.bwr                                                      # colourmap
norm_area = Normalize(-2, 2)                                                # interval of value represented by colourmap
scalarMap_area = ScalarMappable(norm_area, cmap_area)                       # conversion from scalar value to colour

# tension colourbar
cmap_tension = plt.cm.PRGn                                                  # colourmap
norm_tension = Normalize(-2, 2)                                             # interval of value represented by colourmap
scalarMap_tension = ScalarMappable(norm_tension, cmap_tension)              # conversion from scalar value to colour

# orientation colourbar
cmap_orientation = plt.cm.hsv                                               # colourmap
norm_orientation = Normalize(0, 2*np.pi)                                    # interval of value represented by colourmap
scalarMap_orientation = ScalarMappable(norm_orientation, cmap_orientation)  # conversion from scalar value to colour

# resize figure with colorbars
# (https://github.com/matplotlib/matplotlib/issues/15010#issuecomment-524438047)
def _measure_fig(ax):
    ax_size = ax.get_position().size.copy()
    try:
        fig_width, fig_height = (
            lambda s: (s.width(), s.height()))(
            ax.figure.canvas.window().size())
    except AttributeError:
        fig_width, fig_height = None, None
    return ax_size, fig_width, fig_height
def _resize_fig(ax, ax_size, fig_width, fig_height):
    r = ax_size/ax.get_position().size  # rescaling factors
    # rescale
    ax.figure.set_size_inches(ax.figure.get_size_inches()*r)
    try:
        ax.figure.canvas.window().setFixedSize(
            int(fig_width*r.max()), fig_height)
    except AttributeError: pass
    # measure again
    return _measure_fig(ax)

