"""
Define functions to initialise a vertex model.
"""

from cells.bind import VertexModel
from cells.read import Read
from cells.plot import plot
from cells.keratin import plot_keratin
from cells.init import A0
from cells import __path__

import pickle

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,\
    BooleanOptionalAction

import numpy as np
import matplotlib.pyplot as plt

import os
import sys

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-N", type=int, default=400)
parser.add_argument("-p0", type=float, default=3.81,help="dimensionless target perimeter of cell")
 # active forces
parser.add_argument("-taup", type=float, default=1e0,help="persistence time of active force")
# active Brownian force
parser.add_argument("-v0", type=float, default=0.1,help="vertex self-propulsion velocity")
parser.add_argument("-Fpull", type=float, default=1.0,help="outer pulling force")


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
for k in range(10):
    print(k)
    vm.nintegrate(100, 0.01, 0.1, 0.11)
    plot(vm, fig=fig, ax=ax)
plt.close(fig)

vm.removeVertexForce("abp")
vm.removeVertexForce("perimeter")
vm.removeVertexForce("area")

"""
 |      addKeratinModel(self: cells.bind.VertexModel, name: str, K: float, taur: float, Gamma: float, p0: float, l0: float, alpha: float, kth: float, tau: float, sigma: float, ron: float, k0: float, pr0: float) -> None
 |      
 |      Add keratin model.
 |      
 |      Parameters
 |      ----------
 |      name : str
 |          Unique name for the force.
 |      K : float
 |          Area elasticity.
 |      taur : float
 |          Target area relaxation time scale.
 |      Gamma : float
 |          Perimeter elasticity.
 |      p0 : float
 |          Target shape index.
 |      l0 : float
 |          Bond rest length.
 |      alpha : float
 |          Bond elasticity per keratin concentration above threshold.
 |      kth : float
 |          Bond keratin concentration threshold.
 |      tau : float
 |          Keratin concentration evolution time scale.
 |      sigma : float
 |          Keratin concentration noise standard deviation.
 |      ron : float
 |          Keratin concentration on-rate evolution time rate
 |          (=1/tauon).
 |      k0 : float
 |          Keratin concentration off-rate inverse pressure constant.
 |      pr0 : float
 |          Keratin concentration off-rate inflection pressure.
 """

vm.addKeratinModel("keratin",1,100,1,4.2,1,0.2,0,10,0.05,0,5,0.3)
vm.addEdgePullForce("pull",args.Fpull)
fig, ax = plot_keratin(vm)
plt.ion()
plt.show()

#for k in range(100):
while True:
    vm.nintegrate(100, 0.01, 0.1, 0.11)
    plot_keratin(vm, fig=fig, ax=ax)

