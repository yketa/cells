"""
Module plot provides objects and functions to be used in matplotlib plots.
"""

import numpy as np

from collections import OrderedDict

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
cmap = plt.cm.jet

import seaborn as sns

# DEFAULT VARIABLES

_markers = (                        # default markers list
    # mpl.markers.MarkerStyle.filled_markers)
    'o', '^', 's', '*', 'X', 'D', '8', 'v', '<', '>', 'h', 'H', 'p', 'd', 'P')
if mpl.__version__ >= '3':
    _linestyles = (
        (0, ()),                    # solid
        (0, (5, 1)),                # densely dashed
        (0, (1, 1)),                # densely dotted
        (0, (3, 1, 1, 1)),          # densely dashdotted
        (0, (3, 1, 1, 1, 1, 1)),    # densely dashdotdotted
        (0, (5, 5)),                # dashed
        (0, (1, 5)),                # dotted
        (0, (3, 5, 1, 5)),          # dashdotted
        (0, (3, 5, 1, 5, 1, 5)),    # dashdotdotted
        (0, (5, 10)),               # loosely dashed
        (0, (1, 10)),               # loosely dotted
        (0, (3, 10, 1, 10)),        # loosely dashdotted
        (0, (3, 10, 1, 10, 1, 10))  # loosely dashdotdotted
        )
else:
    _linestyles = (
        '-',    # solid
        '--',   # dashed
        '-.',   # dash-dotted
        ':'     # dotted
        )

# LINE STYLE

def list_colormap(value_list, colormap='colorblind', sort=True):
    """
    Creates hash table of colors from colormap, defined according to value_list
    index, with value_list elements as keys.

    Parameters
    ----------
    value_list : list
        List of values.
    colormap : matplotlib colormap or seaborn color palette
        Colormap or color palette to use. (default: 'colorblind')
    sort : bool
        Sort list of values before assigning colors. (default: True)

    Returns
    -------
    colors : hash table
        Hash table of colors.
    """

    value_list = list(OrderedDict.fromkeys(value_list))
    if sort: value_list = sorted(value_list)

    try:    # matplotlib colormap

        cmap = plt.get_cmap(colormap)                               # colormap
        norm = colors.Normalize(vmin=0, vmax=len(value_list) + 1)   # normalise colormap according to list index
        scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)        # associates scalar to color

        return {value_list[index]: scalarMap.to_rgba(index + 1)
            for index in range(len(value_list))}

    except ValueError:  # seaborn palette

        return {value: color
            for value, color in zip(
                value_list,
                sns.color_palette(colormap, len(value_list)))}

def list_markers(value_list, marker_list=_markers, sort=True):
    """
    Creates hash table of markers from markers_list, defined according to
    value_list index, with value_list elements as keys.

    Parameters
    ----------
    value_list : list
        List of values.
    marker_list : list of matplotlib markers
        List of markers to use. (default: cells.plot._markers)
    sort : bool
        Sort list of values before assigning markers. (default: True)

    Returns
    -------
    markers : hash table
        Hash table of markers.
    """

    value_list = list(OrderedDict.fromkeys(value_list))
    if sort: value_list = sorted(value_list)

    return {value_list[index]: marker_list[index]
        for index in range(len(value_list))}

def list_linestyles(value_list, linestyle_list=_linestyles, sort=True):
    """
    Creates hash table of line styles from linestyle_list, defined according to
    value_list index, with value_list elements as keys.

    Parameters
    ----------
    value_list : list
        List of values.
    linestyle_list : list of matplotlib line styles
        List of line styles to use.
        (default: cells.plot._linestyles)
    sort : bool
        Sort list of values before assigning line styles. (default: True)

    Returns
    -------
    linestyles : hash table
        Hash table of line styles.
    """

    value_list = list(OrderedDict.fromkeys(value_list))
    if sort: value_list = sorted(value_list)

    return {value_list[index]: linestyle_list[index]
        for index in range(len(value_list))}

# CONTOUR PLOT

def contours(x, y, z, vmin=None, vmax=None, contours=20, cmap=plt.cm.jet,
        colorbar_position='right', colorbar_orientation='vertical'):
    """
    Plot contours from 3D data.

    Parameters
    ----------
    x : (*,) float array-like
        x-axis data.
    y : (**,) or (*, **) float array-like
        y-axis data.
    z : (*, **) float array-like
        z-axis data to represent with color map.
    vmin : float or None
        Minimum value for the colorbar. (default: None)
        NOTE: if vmin == None then min(z) is taken.
    vmax : float or None
        Maximum value for the colorbar. (default: None)
        NOTE: if vmax == None then max(z) is taken.
    contours : int
        Number of contour lines. (default: 20)
        (see matplotlib.pyplot.tricontourf)
    cmap : matplotlib colorbar
        Matplotlib colorbar to be used. (default: matplotlib.pyplot.cm.jet)
    colorbar_position : string
        Position of colorbar relative to axis. (default: 'right')
    colorbar_orientation : string
        Orientation of colorbar. (default: 'vertical')

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes.Axes
        Axes object.
    colorbar : matplotlib.colorbar
        Colorbar.
    """

    x = np.array(x)
    assert x.ndim == 1
    y = np.array(y)
    if y.ndim == 1: y = np.full((x.size, y.size), fill_value=y)
    assert y.ndim == 2 and y.shape[0] == x.size
    z = np.array(z)
    assert z.shape == y.shape

    vmin = vmin if vmin != None else z.min()
    vmax = vmax if vmax != None else z.max()
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    fig, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(colorbar_position, size='5%', pad=0.05)
    colorbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
        orientation=colorbar_orientation)

    ax.tricontourf(
        *np.transpose([[x[i], y[i, j], z[i, j]]
            for i in range(x.size) for j in range(y[i].size)]),
        contours, cmap=cmap, norm=norm)

    return fig, ax, colorbar

# DISPLAY FITTING LINE

class FittingLine:
    """
    Provided a matplotlib.axes.Axes object, this object:
    > draws a staight line on the corresponding figure, either in a log-log
    (powerlaw fit) or in a lin-log (exponential fit) plot,
    > displays underneath the figure a slider which controls the slope of the
    line, the slider can be hidden / shown by scrolling,
    > enables switching between powerlaw and exponential fit at double click,
    > shows fitting line expression in legend.

    Clicking on the figure updates the position of the line such that it passes
    through the clicked point.

    Instances
    ---------

    FittingLine.ax : matplotlib.axes.Axes object
        Plot Axes object.

    FittingLine.x_fit : string
        x data name in legend.
    FittingLine.y_fit : string
        y data name in legend.
    FittingLine.color : any matplotlib color
        Color of fitting line.
    FittingLine.linestyle : any matplotlib linestyle
        linestyle of fitting line

    FittingLine.x0 : float
        x coordinate of clicked point
    FittingLine.y0 : float
        y coordinate of clicked point
    FittingLine.slope : float
        Slope of fitting line.

    FittingLine.line : matplotlib.lines.Line2D object
        Line2D representing fitting line.

    FittingLine.slider : matplotlib Slider widget
        Slope slider.

    FittingLine.law : string (either 'Powerlaw' or 'Exponential')
        Fitting line law.
    FittingLine.func : function (either Powerlaw or Exponential)
        Fitting line function.
    """

    def __init__(self, ax, slope, slope_min=None, slope_max=None,
        color='black', linestyle='--', slider=True,
        legend=True, exp_format='{:.2e}', font_size=None,
        legend_frame=True, handlelength=None, **kwargs):
        """
        Parameters
        ----------
        ax : matplotlib.axes.Axes object
            Axes object on which to draw fitting line.
        slope : float
            Initial slope of fitting line in log-log plot.
        slope_min : float
            Minimum slope of fitting line for slider.
            NOTE: if slope_min=None, then slope_min is taken to be slope.
            DEFAULT: None
        slope_max : float
            Maximum slope of fitting line for slider.
            NOTE: if slope_max=None, then slope_max is taken to be slope.
            DEFAULT: None
        color : any matplotlib color
            Color of fitting line.
            DEFAULT: black
        linestyle : any matplotlib line style
            Line style of fitting line.
            DEFAULT: --
        slider : bool
            Display slider for slope.
            DEFAULT: True
        legend : bool
            Display legend.
            DEFAULT: True
        exp_format : string
            Exponent string format in legend.
            NOTE: Only if legend=True.
            DEFAULT: {:.2e}
        font_size : float
            Legend font size.
            NOTE: if font_size=None, the font size is not imposed.
            DEFAULT: None
        legend_frame : bool
            Display legend frame.
            DEFAULT: True
        handlelength : float
            Horizontal line length in legend.
            DEFAULT: None

        Optional keyword arguments
        --------------------------
        x_fit : string
            Custom name of x data for fitting line expression in legend.
        y_fit : string
            Custom name of y data for fitting line expression in legend.
        """

        self.ax = ax                # Axes object
        plt.sca(self.ax)            # set current axis
        self.ax.set_yscale('log')   # setting y-axis on logarithmic scale

        self.x_fit = (kwargs['x_fit'] if 'x_fit' in kwargs
            else self.ax.get_xlabel()).replace('$', '') # x data name in legend
        self.y_fit = (kwargs['y_fit'] if 'y_fit' in kwargs
            else self.ax.get_ylabel()).replace('$', '') # y data name in legend
        self.color = color                              # color of fitting line
        self.linestyle = linestyle                      # linestyle of fitting line

        self.x0 = np.exp(np.ma.log(self.ax.get_xlim()).mean())  # x coordinate of clicked point
        self.y0 = np.exp(np.ma.log(self.ax.get_ylim()).mean())  # y coordinate of clicked point
        self.slope = slope                                      # slope of fitting line

        self.line, = self.ax.plot([], [], label=' ',
            color=self.color, linestyle=self.linestyle) # Line2D representing fitting line

        self.display_legend = legend                                # display legend
        if self.display_legend:
            self.x_legend = np.mean(self.ax.get_xlim())             # x coordinate of fitting line legend
            self.y_legend = np.mean(self.ax.get_ylim())             # y coordinate of fitting line legend
            self.legend = plt.legend(handles=[self.line], loc=10,
                bbox_to_anchor=(self.x_legend, self.y_legend),
                bbox_transform=self.ax.transData,
                frameon=legend_frame, handlelength=handlelength)    # fitting line legend
            self.set_fontsize(font_size)                            # set legend font size
            self.legend_artist = self.ax.add_artist(self.legend)    # fitting line legend artist object
            self.legend_artist.set_picker(10)                       # epsilon tolerance in points to fire pick event
        self.on_legend = False                                      # has the mouse been clicked on fitting line legend
        self.exp_format = exp_format                                # exponent string format in legend

        self.display_slider = slider                    # display slider
        if self.display_slider:
            self.slider_ax = make_axes_locatable(self.ax).append_axes(
                'bottom', size='5%', pad=0.6)           # slider Axes
            self.slider = Slider(self.slider_ax, 'slope',
                slope_min if slope_min != None else slope,
                slope_max if slope_max != None else slope,
                valinit=slope)                          # slider
            self.slider.on_changed(self.update_slope)   # call self.update_slope when slider value is changed

        self.law = 'powerlaw'       # fitting line law
        self.update_law()           # initialises fitting line function, updates figure and sets legend

        self.cid_click = self.line.figure.canvas.mpl_connect(
            'button_press_event', self.on_click)        # call on click on figure
        self.cid_pick = self.line.figure.canvas.mpl_connect(
            'pick_event', self.on_pick)                 # call on artist pick on figure
        self.cid_release = self.line.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)    # call on release on figure
        self.cid_scroll = self.line.figure.canvas.mpl_connect(
            'scroll_event', self.on_scroll)             # call on scroll

    def set_fontsize(self, font_size):
        """
        Set legend font size.

        Parameters
        ----------
        font_size : float
            Legend font size.
            NOTE: if font_size=None, the font size is not changed.
        """

        self.font_size = font_size
        if self.font_size != None:
            self.legend.get_texts()[0].set_fontsize(self.font_size) # set legend font size

    def on_click(self, event):
        """
        Executes on click.

        Double click switches between powerlaw and exponential laws and updates
        figure.
        Simple click makes fitting line pass through clicked point and updates
        figure.
        """

        if event.inaxes != self.ax:   # if Axes instance mouse is over is different than figure Axes
            return

        elif self.on_legend:    # if fitting line legend is being dragged
            return

        elif event.dblclick:    # if event is a double click
            self.update_law()   # update fitting line law (and update figure)

        else:
            self.x0 = event.xdata   # x coordinate of clicked point
            self.y0 = event.ydata   # y coordinate of clicked point
            self.draw()             # update figure

    def on_pick(self, event):
        """
        Executes on picking.

        Fitting line legend can be moved if dragged.
        """

        if self.display_legend == False: return

        if event.artist == self.legend_artist:  # if fitting line legend is clicked
            self.on_legend = True               # fitting line legend has been clicked

    def on_release(self, event):
        """
        Executes on release.

        Moves fitting line legend to release position.
        """

        if self.display_legend == False: return

        if not(self.on_legend): return      # if fitting line legend has not been clicked
        self.x_legend = event.xdata         # x coordinate of fitting line legend
        self.y_legend = event.ydata         # y coordinate of fitting line legend
        self.legend.set_bbox_to_anchor(bbox=(self.x_legend, self.y_legend),
            transform=self.ax.transData)    # move legend to release point
        self.line.figure.canvas.draw()      # updates legend
        self.on_legend = False              # fitting line legend has been released

    def on_scroll(self, event):
        """
        Executes on scroll.

        Hide slider if slider is visible, and vice versa.
        """

        if not(self.display_slider): return

        self.slider_ax.set_visible(self.slider_ax.get_visible() == False)   # hide or show slider Axes
        self.line.figure.canvas.draw()                                      # updates figure

    def update_slope(self, val):
        """
        Set fitting line slope according to slider value and updates figure.
        """

        self.slope = self.slider.val    # updates slope of fitting line
        self.update_legend()            # updates legend and figure

    def update_law(self):
        """
        Switches between powerlaw and exponential laws and updates figure.
        """

        self.law = ['powerlaw', 'exponential'][self.law == 'powerlaw']  # switches between powerlaw and exponential
        self.func = {'powerlaw': _powerlaw,
            'exponential': _exponential}[self.law]                      # fitting line function
        self.ax.set_xscale(['linear', 'log'][self.law == 'powerlaw'])   # set x-axis scale according to fitting law
        self.update_legend()                                            # updates legend and figure

    def update_legend(self):
        """
        Updates fitting line legend.
        """

        if self.law == 'powerlaw':
            self.line.set_label(r'$%s \propto %s^{%s}$' % (self.y_fit,
                self.x_fit, self.exp_format.format(self.slope)))    # fitting line label
        elif self.law == 'exponential':
            self.line.set_label(r'$%s \propto e^{%s%s}$' % (self.y_fit,
                self.exp_format.format(self.slope), self.x_fit))    # fitting line label

        if self.display_legend == True:
            self.legend.get_texts()[0].set_text(self.line.get_label())  # updates fitting line legend
        self.draw()                                                     # updates figure

    def draw(self):
        """
        Updates figure with desired fitting line.
        """

        self.line.set_data(self.ax.get_xlim(), list(map(
            lambda x: self.func(self.x0, self.y0, self.slope, x),
            self.ax.get_xlim()
            )))                           # line passes through clicked point according to law
        self.line.figure.canvas.draw()    # updates figure

def _powerlaw(x0, y0, slope, x):
    """
    From point (x0, y0) and parameter slope, returns y = f(x) such that:
    > f(x) = a * (x ** slope)
    > f(x0) = y0

    Parameters
    ----------
    x0, y0, slope, x : float

    Returns
    -------
    y = f(x) : float
    """

    return y0 * ((x/x0) ** slope)

def _exponential(x0, y0, slope, x):
    """
    From point (x0, y0) and parameter slope, returns y = f(x) such that:
    > f(x) = a * exp(x * slope)
    > f(x0) = y0

    Parameters
    ----------
    x0, y0, slope, x : float

    Returns
    -------
    y = f(x) : float
    """

    return y0 * np.exp((x - x0) * slope)
