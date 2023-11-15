from cells.system import VertexModel

import numpy as np

from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize as ColorsNormalise
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

# PARAMETERS

parser = ArgumentParser()
# physics
parser.add_argument("-n", type=int, default=6,
    help="Number of vertices in each direction (must be multiple of 6).")
parser.add_argument("-v0", type=float, default=1e-1,
    help="Vertex self-propulsion velocity.")
parser.add_argument("-Dr", type=float, default=1e-1,
    help="Vertex propulsion rotational diffusion constant.")
parser.add_argument("-p0", type=float, default=3.81,
    help="Dimensionless target perimeter of cell.")
# algorithm
parser.add_argument("-seed", type=int, default=0,
    help="Random number generator seed.")
# integration
parser.add_argument("-dt", type=float, default=1e-3,
    help="Intergation time step.")
parser.add_argument("-delta", type=float, default=0.5,
    help="Length below which to perform T1.")
parser.add_argument("-epsilon", type=float, default=0.1,
    help="Create junction with length epsilon above threshold after T1.")
parser.add_argument("-iterations", type=int, default=1000,
    help="Number of iterations between plots.")

args = parser.parse_args()

# INITIALISATION

m = VertexModel(args.seed, args.v0, args.Dr, args.p0)
m.initRegularTriangularLattice(size=args.n)

# cax = make_axes_locatable(m.ax).append_axes('right', size='5%', pad=0.05)
# cmap = plt.cm.PiYG
# norm = ColorsNormalise(-1, 1)
# colormap = mpl.colorbar.ColorbarBase(cax,
#     cmap=cmap, norm=norm, orientation='vertical')
# scalarMap = ScalarMappable(norm, cmap)

# RUN

plt.ion()
plt.show()
while True:
    # integrate
    m.integrate(
        args.iterations, dt=args.dt, delta=args.delta, epsilon=args.epsilon)
    # plot
    m.plot()
    m.fig.suptitle('t=%.3f' % m.time)

