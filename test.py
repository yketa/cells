from cells.system import ModelSystem

import numpy as np

from argparse import ArgumentParser

import matplotlib.pyplot as plt

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

m = ModelSystem(args.seed, args.v0, args.Dr, args.p0)
m.initRegularTriangularLattice(size=args.n)

# RUN

plt.ion()
plt.show()
while True:
    # integrate
    m.integrate(
        args.iterations, dt=args.dt, delta=args.delta, epsilon=args.epsilon)
#     m.checkMesh()
    # plot
    m.plot()

