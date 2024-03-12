"""
Define functions to initialise a vertex model.
"""

from cells.bind import VertexModel
from cells.read import Read
from cells.plot import plot, plot_forces
from cells.keratin import plot_keratin
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
parser.add_argument("-Fpull", type=float, default=1.0,help="outer pulling force")


args = parser.parse_args()
vm = VertexModel(0)
vm.initOpenRegularHexagonalLattice(nCells=args.N)
vm.addPerimeterForce("perimeter",1, args.p0*np.sqrt(A0))
vm.addAreaForce("area",1, A0)
vm.addActiveBrownianForce("abp",args.v0, args.taup)
vm.addEdgePullForce("pull",-0.1)

#fig, ax = plot(vm)
#plt.ion()
#plt.show()
for k in range(100):
    print(k)
    vm.nintegrate(100, 0.01, 0.1, 0.11)
    #plot(vm, fig=fig, ax=ax)

vm.removeVertexForce("abp")
vm.removeVertexForce("perimeter")
vm.removeVertexForce("area")
#plt.close(fig)

"""addKeratinModel(...)
 |      addKeratinModel(self: cells.bind.VertexModel, name: str, K: float, A0: float, Gamma: float, P0: float, l0: float, alpha: float, kth: float, tau: float, sigma: float, tauon: float, k0: float, p0: float) -> None
 |      
 |      Add keratin model.
 |      
 |      Parameters
 |      ----------
 |      name : str
 |          Unique name for the force.
 |      K : float
 |          Area elasticity.
 |      A0 : float
 |          Target area.
 |      Gamma : float
 |          Perimeter elasticity.
 |      P0 : float
 |          Target perimeter.
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
 |      tauon : float
 |          Keratin concentration on-rate evolution time scale.
 |      k0 : float
 |          Keratin concentration off-rate inverse pressure constant.
 |      p0 : float
 |          Keratin concentration off-rate inflection pressure.
 """

vm.addKeratinModel("keratin",1,A0,1,4.2*np.sqrt(A0),1,0.2,0,10,0.05,1000000,5,0.3)
vm.addEdgePullForce("pull",args.Fpull)
fig, ax = plot_keratin(vm)
plt.ion()
plt.show()

for k in range(100):
    vm.nintegrate(100, 0.01, 0.1, 0.11)
    plot_keratin(vm, fig=fig, ax=ax)