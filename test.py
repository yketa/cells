from cells.bind import VertexModel
from cells.vm import plot

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,\
    BooleanOptionalAction

import numpy as np

import matplotlib.pyplot as plt

if __name__ == "__main__":

    # PARAMETERS

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    # physics
    parser.add_argument("-N", type=int, default=9,
        help="[open] total number of cells (! square number)")
    parser.add_argument("-n", type=int, default=6,
        help="[close] number of vertices in each direction (! multiple of 6)")
    parser.add_argument("-v0", type=float, default=1e-1,
        help="vertex self-propulsion velocity")
    parser.add_argument("-Dr", type=float, default=1e-1,
        help="vertex propulsion rotational diffusion constant")
    parser.add_argument("-p0", type=float, default=3.81,
        help="dimensionless target perimeter of cell")
    parser.add_argument("-open", "-o",
        action=BooleanOptionalAction,
        help="turn on specific checks for boundary vertices")
    # algorithm
    parser.add_argument("-seed", type=int, default=0,
        help="random number generator seed")
    # integration
    parser.add_argument("-dt", type=float, default=1e-3,
        help="intergation time step")
    parser.add_argument("-delta", type=float, default=0.1,
        help="length below which to perform T1")
    parser.add_argument("-epsilon", type=float, default=0.1,
        help="create junction with length epsilon above threshold after T1")
    parser.add_argument("-iterations", type=int, default=1000,
        help="number of iterations between plots")

    args = parser.parse_args()

    # INITIALISATION

    vm = VertexModel(args.seed)

    if args.open:
        vm.initOpenRegularHexagonalLattice(nCells=args.N)
#         vm.initOpenRegularTriangularLattice(size=args.n)
    else:
        vm.initRegularTriangularLattice(size=args.n)
    fig, ax = plot(vm)

    # FORCES

    A0 = (3./2.)/np.tan(np.pi/6.)   # area of a regular hexagon with edge length 1
    vm.addPerimeterForce("perimeter", 1, args.p0*np.sqrt(A0))
    vm.addAreaForce("area", 1, A0)
    vm.addActiveBrownianForce("abp", args.v0, 1./args.Dr)

    # RUN

    plt.ion()
    plt.show()
    while True:
        # integrate
        vm.nintegrate(args.iterations,
            dt=args.dt, delta=args.delta, epsilon=args.epsilon)
        #m.checkMesh()
        # plot
        plot(vm, fig=fig, ax=ax)

