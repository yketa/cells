from cells.bind import VertexModel
from cells.vm import plot

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import matplotlib.pyplot as plt

if __name__ == "__main__":

    # PARAMETERS

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    # physics
    parser.add_argument("-n", type=int, default=6,
        help="number of vertices in each direction (must be multiple of 6)")
    parser.add_argument("-v0", type=float, default=1e-1,
        help="vertex self-propulsion velocity")
    parser.add_argument("-Dr", type=float, default=1e-1,
        help="vertex propulsion rotational diffusion constant")
    parser.add_argument("-p0", type=float, default=3.81,
        help="dimensionless target perimeter of cell")
    # algorithm
    parser.add_argument("-seed", type=int, default=0,
        help="random number generator seed")
    # integration
    parser.add_argument("-dt", type=float, default=1e-3,
        help="intergation time step")
    parser.add_argument("-delta", type=float, default=0.5,
        help="length below which to perform T1")
    parser.add_argument("-epsilon", type=float, default=0.1,
        help="create junction with length epsilon above threshold after T1")
    parser.add_argument("-iterations", type=int, default=1000,
        help="number of iterations between plots")

    args = parser.parse_args()

    # INITIALISATION

    vm = VertexModel(args.seed, args.v0, args.Dr, args.p0)
    vm.initRegularTriangularLattice(size=args.n)
    fig, ax = plot(vm)

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

