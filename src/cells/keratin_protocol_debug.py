"""
Define functions to initialise a vertex model.
"""

from cells.bind import VertexModel
from cells.read import Read
from cells.plot import plot, plot_forces
from cells import __path__

import pickle

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,\
    BooleanOptionalAction

import numpy as np
import matplotlib.pyplot as plt

import os
import sys

A0 = (3./2.)/np.tan(np.pi/6.)           # area of a regular hexagon with edge length 1

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-N", type=int, default=400)
parser.add_argument("-p0", type=float, default=3.81,help="dimensionless target perimeter of cell")
 # active forces
parser.add_argument("-taup", type=float, default=1e0,help="persistence time of active force")
# active Brownian force
parser.add_argument("-v0", type=float, default=0.1,help="vertex self-propulsion velocity")
parser.add_argument("-Fpull", type=float, default=0.1,help="outer pulling force")


args = parser.parse_args()
vm = VertexModel(0)
vm.initOpenRegularHexagonalLattice(nCells=args.N)
vm.addPerimeterForce("perimeter",1, args.p0*np.sqrt(A0))
vm.addAreaForce("area",1, A0)
vm.addActiveBrownianForce("abp",args.v0, args.taup)
vm.addEdgePullForce("pull",-0.1)

fig, ax = plot(vm)
plt.ion()
plt.show()
for k in range(100):
    vm.nintegrate(100, 0.01, 0.1, 0.11)
    plot(vm, fig=fig, ax=ax)

vm.removeVertexForce("abp")
vm.addEdgePullForce("pull",args.Fpull)

for k in range(100):
    vm.nintegrate(100, 0.01, 0.1, 0.11)
    plot(vm, fig=fig, ax=ax)
