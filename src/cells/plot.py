"""
Define objects functions to plot vertex model object.
"""

from cells import init
from cells.bind import VertexModel, angle2, getPercentageKeptNeighbours,\
    getLinesHalfEdge, getLinesJunction, getPolygonsCell, hexagonEdgeLength

import numpy as np
from operator import itemgetter

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm,\
    LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection, LineCollection

# PLOT VERTEX MODEL OBJECT

class WindowClosedException(Exception): pass

def _update_canvas(fig):
    fig.canvas.draw_idle()
    fig.canvas.start_event_loop(0.001)
    if not(plt.fignum_exists(fig.number)):  # throw error when window is closed
        raise WindowClosedException

def plot(vm, fig=None, ax=None, update=True,
    rainbow=None, clear=False, vertex_indices=False, only_set=False):
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
    vertex_indices : bool
        Write indices alongside vertices. (default: False)
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
    rainbow_means_relaxation = True                 # is rainbow plot a relaxation plot?
    if only_set: clear = True
    if rainbow_plot and clear:
        raise ValueError("`rainbow' and `clear' options are exclusionary.")

    # forces parameters

    if "area" in vm.vertexForces:
        A0 = vm.vertexForces["area"].parameters["A0"]
    if "perimeter" in vm.vertexForces:
        if "area" in vm.vertexForces:
            p0 = vm.vertexForces["perimeter"].parameters["P0"]/np.sqrt(A0)
        else:
            P0 = vm.vertexForces["perimeter"].parameters["P0"]
    if "volume" in vm.vertexForces or "linear_volume" in vm.vertexForces:
        if "volume" in vm.vertexForces:
            volume_force = vm.vertexForces["volume"]
        else:
            volume_force = vm.vertexForces["linear_volume"]
        h0 = volume_force.parameters["H0"]/volume_force.parameters["A0"]
    if "boundary_tension" in vm.vertexForces:
        gamma = vm.vertexForces["boundary_tension"].parameters["gamma"]
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

    # integrators parameters

    if "eta" in vm.integrator.parameters:
        eta = vm.integrator.parameters["eta"]

    # initialise figure

    def _set_lim(ax_, vm_):
        if True:
            ax_.set_xlim([0, vm_.systemSize[0]])
            ax_.set_ylim([0, vm_.systemSize[1]])
        else:
            # zoom
            n = 4
            ax_.set_xlim([vm_.systemSize[0]/n, (n - 1)*vm_.systemSize[0]/n])
            ax_.set_ylim([vm_.systemSize[0]/n, (n - 1)*vm_.systemSize[0]/n])
        ax_.set_aspect("equal")

    if fig is None or ax is None:

        plt.ioff()
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)                                         # set figure size
        try: fig.canvas.window().setFixedSize(fig.canvas.window().size())   # set window size
        except AttributeError: pass

        # all force-related colourbars
        if not(rainbow_plot) and not(clear):

            ax_size, fig_width, fig_height = _measure_fig(ax)

            if "volume" in vm.vertexForces:
                cbar_volume = fig.colorbar(
                    mappable=scalarMap_area, ax=ax,
                    shrink=0.75, pad=0.01)
                cbar_volume.set_label(
                    r"$(h_i - \left<h_i\right>)/\mathrm{std}(h_i)$",
                    rotation=270, labelpad=_cbar_labelpad)
                # resize
                ax_size, fig_width, fig_height = (
                    _resize_fig(ax, ax_size, fig_width, fig_height))

            if "area" in vm.vertexForces and not("volume" in vm.vertexForces):
                cbar_area = fig.colorbar(
                    mappable=scalarMap_area, ax=ax,
                    shrink=0.75, pad=0.01)
                cbar_area.set_label(
                    r"$(A_i - \left<A_i\right>)/\mathrm{std}(A_i)$",
                    rotation=270, labelpad=_cbar_labelpad)
                # resize
                ax_size, fig_width, fig_height = (
                    _resize_fig(ax, ax_size, fig_width, fig_height))

            if "out" in vm.halfEdgeForces:
                cbar_tension = fig.colorbar(
                    mappable=scalarMap_tension, ax=ax,
                    shrink=0.75, pad=0.01)
                cbar_tension.set_label(
                    r"$(t_i - \left<t_i\right>)/\mathrm{std}(t_i)$",
                    rotation=270, labelpad=_cbar_labelpad)
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
                        rotation=270, labelpad=_cbar_labelpad)
                    # resize
                    ax_size, fig_width, fig_height = (
                        _resize_fig(ax, ax_size, fig_width, fig_height))

        # relaxation-related colourbar
        elif rainbow_plot and rainbow_means_relaxation:

            ax_size, fig_width, fig_height = _measure_fig(ax)

            cbar_relaxation = fig.colorbar(
                mappable=scalarMap_relaxation, ax=ax,
                shrink=0.75, pad=0.01)
            cbar_relaxation.set_label(
                r"\% lost neighbours",
                rotation=270, labelpad=_cbar_labelpad)
            cbar_relaxation.set_ticks([0/3, 1/3, 2/3, 3/3])
            cbar_relaxation.set_ticklabels(
                [r"$0$", r"$\frac{1}{3}$", r"$\frac{2}{3}$", r"$1$"])
            # resize
            ax_size, fig_width, fig_height = (
                _resize_fig(ax, ax_size, fig_width, fig_height))

        # set figure limits
        _set_lim(ax, vm)
        fig.canvas.mpl_connect("button_press_event",    # reset figure limits on double click
            lambda event: event.dblclick and _set_lim(ax, vm))

    # plot

    plt.sca(ax)
    try:
        # make zoom persistent
        # https://discourse.matplotlib.org/t/how-to-make-zoom-persistent/22663/3
        fig.canvas.toolbar.push_current()
        ax.cla()
        fig.canvas.toolbar.back()
    except AttributeError:
        ax.cla()
        _set_lim(ax, vm)
    if only_set: return fig, ax

    # junctions and half-edges
    lines = LineCollection(getLinesJunction(vm), colors="pink",                 # all junctions
        linewidths=2.5/max(1, np.sqrt(len(vm.vertices))/12))                    # scale junction width with linear system size
    if not(rainbow_plot) and not(clear):
        if ("t0" in locals() or "sT0" in locals()):
            junctions = vm.getHalfEdgeIndicesByType("junction")
            if "t0" in locals():
                tensions = (
                    lambda tension:
                        np.concatenate(list(map(
                            lambda i: [tension[i]]*2,
                            junctions))))(
                    vm.halfEdgeForces["out"].tension)
            elif "sT0" in locals():
                for model in ["model%i" % i for i in range(5)]:
                    try:
                        tensions = (
                            lambda tension: np.concatenate(list(map(
                                lambda i: [tension[i]]*2,
                                junctions))))(
                            vm.halfEdgeForces[model].tension)
                    except:
                        continue
            tensions_std = tensions.std()
            if tensions_std != 0:
                lines.set_color(list(map(
                    lambda s_tension: scalarMap_tension.to_rgba(s_tension),
                    (tensions - tensions.mean())/tensions_std)))
    ax.add_collection(lines)
#     ax.add_collection(LineCollection(                                           # all half-edges
#         getLinesHalfEdge(vm, list(vm.halfEdges)), colors="blue", lw=1))

    # cells
    polygons = PatchCollection(
        list(map(                                   # all cells
            lambda vertices: plt.Polygon(vertices, closed=True),
            getPolygonsCell(vm))),
        facecolors="none", edgecolors="none")
    cells = vm.getVertexIndicesByType("centre")
    if not(clear):
        if rainbow_plot:
            if not(rainbow_means_relaxation):
                positions0 = (
                    lambda vertices: np.array(list(map(
                        lambda i: vertices[i].position[0],
                        cells))))(
                    rainbow.vertices)
                scalarMap_rainbow = ScalarMappable(Normalize(
                        positions0.min(), positions0.max(), cmap_orientation))
                polygons.set_facecolor(list(map(    # colour according to previous position
                    lambda position0: scalarMap_rainbow.to_rgba(position0),
                    positions0)))
            else:
                pct = getPercentageKeptNeighbours(rainbow, vm)
                polygons.set_facecolor(list(map(    # colour according to percentage of lost neighbours
                    lambda i: scalarMap_relaxation.to_rgba(1 - pct[i]),
                    cells)))
        elif "volume_force" in locals():
            heights = (
                lambda height: np.array(list(map(
                    lambda i: height[i],
                    cells))))(
                volume_force.height)
            heights_mean, heights_std = heights.mean(), heights.std()
            areas = np.array(list(map(lambda i: vm.getVertexToNeighboursArea(i), cells)))
            volumes = areas*heights
            print(volumes.min(), volumes.mean(), volumes.max())
            if heights_std != 0:
                polygons.set_facecolor(list(map(    # colour according to height
                    lambda s_height: scalarMap_area.to_rgba(s_height),
                    (heights - heights_mean)/heights_std)))
        elif "A0" in locals():
            areas = np.array(list(map(
                lambda i: vm.getVertexToNeighboursArea(i),
                cells)))
            areas_mean, areas_std = areas.mean(), areas.std()
            if areas_std != 0:
                polygons.set_facecolor(list(map(    # colour according to area
                    lambda s_area: scalarMap_area.to_rgba(s_area),
                    (areas - areas_mean)/areas_std)))
    ax.add_collection(polygons)

    # vertex indices
    if vertex_indices:
        (lambda vertices: list(map(
            lambda i: ax.text(*vertices[i].position, i),
            vertices)))(
        vm.vertices)

    title = (r"$t=%.3f, N_{\mathrm{T}_1}=%.3e, N_{\mathrm{cells}}=%i$"
        % (vm.time, vm.nT1, len(cells)))
    if "eta" in locals():
        title += r"$, \eta=%.1e$" % eta
    if "gamma" in locals():
        title += r"$, \gamma=%.1e$" % gamma
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

def plot_neighbours(vm, fig=None, ax=None, update=True,
    **kwargs):
    """
    Plot vertex model with number of neighbours per cell.

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

    Additional keyword arguments are passed to cells.plot.plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    assert(not("clear" in kwargs) or not(kwargs["clear"]))  # not compatible with clear plot

    set_call = (fig is None or ax is None)
    fig, ax = plot(vm, fig=fig, ax=ax, clear=True, update=False, **kwargs)  # initialise figure and axis

    # compute number of neighbours
    cells = vm.getVertexIndicesByType("centre")
    nNeigh = np.array(list(map(
        lambda i: vm.getNeighbourVertices(i)[0].size,
        cells)))

    # colourbars
    if set_call:
        # measure
        ax_size, fig_width, fig_height = _measure_fig(ax)
        # colourbar
        cbar_neigh = fig.colorbar(
            mappable=scalarMap_neigh, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_neigh.set_label(
            r"number of neighbours",
            rotation=270, labelpad=_cbar_labelpad)
        cbar_neigh.set_ticks([
            min(n_neigh)
                + (0.5 + i)*(max(n_neigh) - min(n_neigh))/(len(n_neigh) - 1.)
            for i in range(len(n_neigh) - 1)])
        cbar_neigh.set_ticklabels(
            [r'$%i$' % i for i in n_neigh[:-1]])
        # resize
        ax_size, fig_width, fig_height = (
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # set colours
    ax.collections[1].set_facecolor(scalarMap_neigh.to_rgba(nNeigh))

    # update canvas
    if update: _update_canvas(fig)

    return fig, ax

def plot_hexatic(vm, fig=None, ax=None, update=True, transparency=True,
    **kwargs):
    """
    Plot vertex model with hexatic bond orientational order parameter per cell.

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
    transparency : bool
        Set transparency to match norm of hexatic bond orientional order
        parameter. (default: True)

    Additional keyword arguments are passed to cells.plot.plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    assert(not("clear" in kwargs) or not(kwargs["clear"]))  # not compatible with clear plot

    set_call = (fig is None or ax is None)
    fig, ax = plot(vm, fig=fig, ax=ax, clear=True, update=False, **kwargs)  # initialise figure and axis

    # compute hexatic bond orientational order parameter
    cells = vm.getVertexIndicesByType("centre")
    psi6 = np.array(itemgetter(*cells)(vm.getPAticOrderParameters(p=6)))

    # colourbars
    if set_call:
        # measure
        ax_size, fig_width, fig_height = _measure_fig(ax)
        # colourbar
        cbar_psi6 = fig.colorbar(
            mappable=scalarMap_orientation, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_psi6.set_label(
            r"$\mathrm{arg}(\psi_{6,i})$",
            rotation=270, labelpad=_cbar_labelpad)
        cbar_psi6.set_ticks(
            [-np.pi, -2*np.pi/3, -np.pi/3, 0,
                np.pi/3, 2*np.pi/3, np.pi])
        cbar_psi6.set_ticklabels(
            [r"$-\pi$", r"$-\frac{2\pi}{3}$", r"$\frac{\pi}{3}$", r"$0$",
                r"$\frac{\pi}{3}$", r"$\frac{2\pi}{3}$", r"$\pi$"])
        # resize
        ax_size, fig_width, fig_height = (
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # set colours
    ax.collections[1].set_facecolor(scalarMap_orientation.to_rgba(  # argument
        angle2(psi6.real, psi6.imag)))
    if transparency: ax.collections[1].set_alpha(np.abs(psi6))      # norm

    # update canvas
    if update: _update_canvas(fig)

    return fig, ax

def plot_translational(vm, fig=None, ax=None, update=True, **kwargs):
    """
    Plot vertex model with orientation of translational order parameter per
    cell.

    This involves a pair of reciprocal vectors of the lattice, whose
    orientation is determined by the argument of the average hexatic order
    parameter and whose norm is 4\\pi/a where a is the distance between two
    centroids in a perfect hexagonal lattice, such that for points within a
    periodic hexagonal lattice the dot product of the reciprocal vectors with
    their position is an integer multiple of 2\\pi.

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

    Additional keyword arguments are passed to cells.plot.plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    assert(not("clear" in kwargs) or not(kwargs["clear"]))  # not compatible with clear plot

    set_call = (fig is None or ax is None)
    fig, ax = plot(vm, fig=fig, ax=ax, clear=True, update=False, **kwargs)  # initialise figure and axis

    # compute reciprocal vector
    a = np.sqrt(3)*hexagonEdgeLength(init.A0)   # centres distance in initial periodic hexagonal lattice
    cells = vm.getVertexIndicesByType("centre")
    vectorsToNeighbours = vm.getVectorsToNeighbouringCells()
    psi = np.mean(list(map(                     # average hexatic order parameter
        lambda i: np.mean(list(map(
            lambda v: np.exp(1j*6*angle2(*v)),
            vectorsToNeighbours[i]))),
        cells)))
    theta = angle2(psi.real, psi.imag)          # orientation of average hexatic order parameter
    q0 = (2*np.pi/((1/2)*a))*np.array(          # reciprocal vector of lattice oriented by hexatic order parameter (note there is a pi/3 rotation in the initial matrix)
        [np.sin(theta + np.pi/3), -np.cos(theta + np.pi/3)])
    q1 = np.array([                             # other reciprocal vector
        -(1/2)*q0[0] + (np.sqrt(3)/2)*q0[1],    # rotation by -2\pi/3
        -(np.sqrt(3)/2)*q0[0] - (1/2)*q0[1]])

    # compute translational order parameter
    pos = vm.getPositions(wrapped=True)
    thetar = np.array(list(map(
        lambda i: (
            lambda psi: angle2(psi.real, psi.imag)%(2*np.pi))(
            (np.exp(1j*np.dot(q0, pos[i])) + np.exp(1j*np.dot(q1, pos[i])))/2),
        cells)))

    # colourbars
    if set_call:
        # measure
        ax_size, fig_width, fig_height = _measure_fig(ax)
        # colourbar
        cbar_psir = fig.colorbar(
            mappable=scalarMap_orientation, ax=ax,
            shrink=0.75, pad=0.01)
        cbar_psir.set_label(
            r"$\mathrm{arg}"
            r"\left(\frac{1}{2}\left["
            r"e^{\mathrm{i}\boldsymbol{q}_0\cdot\boldsymbol{r}_i}"
            r"+ e^{\mathrm{i}\boldsymbol{q}_1\cdot\boldsymbol{r}_i}"
            r"\right]\right)$",
            rotation=270, labelpad=_cbar_labelpad)
        cbar_psir.set_ticks(
            [-np.pi, -2*np.pi/3, -np.pi/3, 0,
                np.pi/3, 2*np.pi/3, np.pi])
        cbar_psir.set_ticklabels(
            [r"$-\pi$", r"$-\frac{2\pi}{3}$", r"$\frac{\pi}{3}$", r"$0$",
                r"$\frac{\pi}{3}$", r"$\frac{2\pi}{3}$", r"$\pi$"])
        # resize
        ax_size, fig_width, fig_height = (
            _resize_fig(ax, ax_size, fig_width, fig_height))

    # set colours
    ax.collections[1].set_facecolor(scalarMap_orientation.to_rgba(thetar))

    # update canvas
    if update: _update_canvas(fig)

    return fig, ax

def plot_velocities(vm, fig=None, ax=None, update=True,
    av_norm=0.2, hide_centres=True, hide_vertices=False, override=None,
    **kwargs):
    """
    Plot vertex model with velocities on vertices.

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
    av_norm : float
        Average norm of arrows. (default: 0.2)
    hide_centres : bool
        Do not display velocities of cell centres. (default: True)
    hide_vertices : bool
        Do not display velocities of vertices. (default: False)
    override : function or None
        Override plotting function. (default: None)
        NOTE: This is expected to return a figure and an axes subplot, and
              to support keyword argument `update'. (see cells.plot.plot)

    Additional keyword arguments are passed to the plotting function.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes subplot.
    """

    assert(not("clear" in kwargs) or not(kwargs["clear"]))  # not compatible with clear plot

    fig, ax = (
        # plotting function
        override if not(override is None) else plot)(
        # plotting argument
        vm, fig=fig, ax=ax, update=False, **kwargs)

    positions = vm.getPositions(wrapped=True)
    velocities = vm.velocities

    if len(velocities) > 0:
        indices = list(
            set(list(velocities)).intersection(set(list(positions))))   # vertices which have both a velocity and a position

        if hide_centres:                                                # remove cell centres
            centres = vm.getVertexIndicesByType("centre")
            indices = list(set(indices) - set(centres))
        else:
            velocities.update(vm.getCentreVelocities())

        if hide_vertices:                                               # remove vertices
            vertices = vm.getVertexIndicesByType("vertex")
            indices = list(set(indices) - set(vertices))

        if len(indices) > 0:
            av_norm = (1./av_norm)*np.sqrt(
                (np.array(list(map(lambda i: velocities[i], indices)))**2)
                    .sum(axis=-1)
                        .mean())

            if av_norm != 0:
                arrows = PatchCollection(
                    list(map(
                        lambda i:
                            plt.arrow(*positions[i],
                                *np.array(velocities[i])/av_norm,
                                width=4e-2, head_width=1.5e-1,
                                length_includes_head=False,
                                color=scalarMap_orientation.to_rgba(
                                    angle2(*velocities[i])),
                                zorder=10),
                        indices)))
                ax.add_collection(arrows)

    # update canvas
    if update: _update_canvas(fig)

    return fig, ax

# COLOURBARS

# resize figure with colourbars
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

_cbar_labelpad = 40 # label spacing from axis

def _make_cycle(colours):
    """
    Cyclic colourmap from linear interpolation between colours.

    https://wildsilicon.com/blog/2018/cyclical-colormaps/

    Parameters
    ----------
    colours : (*, 3) float array-like
        Array of RGB colours.

    Returns
    -------
    cmap : matplotlib.colors.ListedColormap
        Cyclic colourmap.
    """

    colours = np.array(colours)
    assert(colours.ndim == 2 and (colours.shape[1] in (3, 4)))
    N = colours.shape[0]
    C = 256

    def _mod(a, b): return min(abs(a)%abs(b), abs(a)%-abs(b), key=abs)
    def _dist(x, y, L): return np.abs(_mod(x - y, L))
    def _kernel(i, n, N): return max(0, 1 - _dist(i/C, (1/2 + n)/N, 1)*N)

    cycle = []
    for i in range(C):
        weights = np.array(list(map(
            lambda n: _kernel(i, n, N),
            range(N))))
        cycle += [np.sum(
            list(map(
                lambda c, w: w*c/weights.sum(),
                *(colours, weights))),
            axis=0)]

    return ListedColormap(cycle)

# area colourbar
cmap_area = (                                                               # colourmap
    LinearSegmentedColormap.from_list("aroace", (                           # https://en.wikipedia.org/wiki/Pride_flag#/media/File:Aroace_flag.svg
        (0/4, (0.125, 0.220, 0.337)),
        (1/4, (0.384, 0.682, 0.863)),
        (2/4, (1.000, 1.000, 1.000)),
        (3/4, (0.925, 0.804, 0.000)),
        (4/4, (0.886, 0.549, 0.000))))
    or plt.cm.bwr)
norm_area = Normalize(-2, 2)                                                # interval of value represented by colourmap
scalarMap_area = ScalarMappable(norm_area, cmap_area)                       # conversion from scalar value to colour

# tension colourbar
cmap_tension = (                                                            # colourmap
    LinearSegmentedColormap.from_list("abrosexual", (                       # https://en.wikipedia.org/wiki/Pride_flag#/media/File:Abrosexual_flag.svg
        (0/4, (0.851, 0.267, 0.431)),
        (1/4, (0.906, 0.588, 0.718)),
        (2/4, (1.000, 1.000, 1.000)),
        (3/4, (0.706, 0.894, 0.800)),
        (4/4, (0.396, 0.761, 0.525))))
    or LinearSegmentedColormap.from_list("new_gay_men", (                   # https://en.wikipedia.org/wiki/Gay_men%27s_flags#/media/File:New_Gay_Pride_Flag.svg
        (0/6, (0.239, 0.102, 0.471)),
        (1/6, (0.314, 0.286, 0.800)),
        (2/6, (0.482, 0.678, 0.886)),
        (3/6, (1.000, 1.000, 1.000)),
        (4/6, (0.596, 0.910, 0.757)),
        (5/6, (0.149, 0.808, 0.667)),
        (6/6, (0.027, 0.553, 0.439))))
    or LinearSegmentedColormap.from_list("genderqueer", (                   # https://en.wikipedia.org/wiki/Non-binary_gender#/media/File:Genderqueer_Pride_Flag.svg
        (0/2, (0.290, 0.506, 0.130)),
        (1/2, (1.000, 1.000, 1.000)),
        (2/2, (0.710, 0.494, 0.863))))
    or plt.cm.PRGn)
norm_tension = Normalize(-2, 2)                                             # interval of value represented by colourmap
scalarMap_tension = ScalarMappable(norm_tension, cmap_tension)              # conversion from scalar value to colour

# orientation colourbar
cmap_orientation = (                                                        # colourmap
    _make_cycle((                                                           # https://en.wikipedia.org/wiki/Rainbow_flag_(LGBT)#/media/File:Gay_Pride_Flag.svg
        (0.467, 0.000, 0.533),
        (0.000, 0.298, 1.000),
        (0.008, 0.506, 0.129),
        (1.000, 0.933, 0.000),
        (1.000, 0.553, 0.000),
        (0.898, 0.000, 0.000)))
    or plt.cm.hsv)
norm_orientation = Normalize(-np.pi, np.pi)                                 # interval of value represented by colourmap
scalarMap_orientation = ScalarMappable(norm_orientation, cmap_orientation)  # conversion from scalar value to colour

# neighbours colourbar
n_neigh = (lambda n: range(6 - n, 6 + n + 2))(3)                            # interval of connectivity
cmap_neigh = ListedColormap([                                               # colourmap
    ScalarMappable(
        norm=Normalize(vmin=min(n_neigh), vmax=max(n_neigh) - 1),
            cmap=(LinearSegmentedColormap.from_list("lesbian_2018", (       # https://en.wikipedia.org/wiki/Lesbian_flags#/media/File:Lesbian_pride_flag_2018.svg
                (0/6, (0.639, 0.008, 0.384)),
#                 (1/6, (0.710, 0.337, 0.565)),
#                 (2/6, (0.820, 0.384, 0.643)),
                (3/6, (1.000, 1.000, 1.000)),
#                 (4/6, (1.000, 0.604, 0.337)),
#                 (5/6, (0.937, 0.463, 0.153)),
                (6/6, (0.835, 0.176, 0.000))))
            or plt.cm.PiYG)
        ).to_rgba(x)
    for x in n_neigh])
norm_neigh = BoundaryNorm(n_neigh, len(n_neigh))                            # interval of value represented by colourmap
scalarMap_neigh = ScalarMappable(cmap=cmap_neigh, norm=norm_neigh)          # conversion from scalar value to colour

# relaxation colourbar
cmap_relaxation = ListedColormap([                                          # colourmap
    ScalarMappable(
        norm=Normalize(vmin=0, vmax=1),
            cmap=(LinearSegmentedColormap.from_list("polysexual", (         # https://en.wikipedia.org/wiki/Pride_flag#/media/File:Polysexuality_Pride_Flag.svg
                (0/2, (0.110, 0.573, 0.965)),
                (1/2, (0.027, 0.835, 0.412)),
                (2/2, (0.965, 0.110, 0.725))))
            or plt.cm.cool),
        ).to_rgba(x)
    for x in (0/2, 1/2, 2/2)])
norm_relaxation = Normalize(0, 1)                                           # interval of value represented by colourmap
scalarMap_relaxation = ScalarMappable(norm_relaxation, cmap_relaxation)     # conversion from scalar value to colour

