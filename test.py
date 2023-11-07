from cells.system import VertexModel

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize as ColorsNormalise
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

m = VertexModel()
m.initRegularTriangularLattice(size=6)

cax = make_axes_locatable(m.ax).append_axes('right', size='5%', pad=0.05)
cmap = plt.cm.PiYG
norm = ColorsNormalise(-1, 1)
colormap = mpl.colorbar.ColorbarBase(cax, cmap, norm, orientation='vertical')
scalarMap = ScalarMappable(norm, cmap)

def update(iterations=100, dt=1e-3, A=5e-2):

    for iteration in range(iterations):
        m.integrate(dt=dt, delta=0.5, epsilon=0.1)

    # plot

    m.plot()
    m.fig.suptitle('t=%s' % m.time)

anim = animation.FuncAnimation(m.fig, update, repeat=True, interval=0)
plt.show()
