"""
Initialise vertex model.
"""

from cells.bind import VertexModel

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,\
    BooleanOptionalAction

import numpy as np

import os
import sys

A0 = (3./2.)/np.tan(np.pi/6.)           # area of a regular hexagon with edge length 1
script = os.path.basename(sys.argv[0])  # name of invoking script

def init_vm():
    """
    Initilise VertexModel object with command line arguments.

    Returns
    -------
    args : argparse.Namespace
        Arguments.
    vm : VertexModel
        New instance from command line arguments.
    """

    args = parse_args() # command line arguments

    # INITIALISATION

    vm = VertexModel(args.seed)

    if args.periodic:
        vm.initRegularTriangularLattice(size=args.N)
    else:
        vm.initOpenRegularHexagonalLattice(nCells=args.N)
#         vm.initOpenRegularTriangularLattice(size=args.N)

    # FORCES

    if args.perimeter:
        vm.addPerimeterForce("perimeter",
            1, args.p0*np.sqrt(A0))
    if args.area:
        vm.addAreaForce("area",
            1, A0)
    if args.abp:
        vm.addActiveBrownianForce("abp",
            args.v0, args.taup)
    if args.out:
        vm.addOrnsteinUhlenbeckTension("out",
            args.t0, args.st0, args.taup)
    if args.model0:
        vm.addModel0("model0",
            args.Gamma, args.p0*np.sqrt(A0), args.sigma, args.taup)

    return args, vm

def parse_args():
    """
    Parse command line arguments for vertex model.

    Returns
    -------
    args : argparse.Namespace
        Arguments.
    """

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    # BOX
    parser.add_argument("-N", type=int, default=36,
        help=
            "[--no-periodic (default)] number of cells (! square number) "
            "[-periodic] number of vertices in each direction"
            "(! multiple of 6)")
    parser.add_argument("--periodic", "-periodic", "-p",
        action=BooleanOptionalAction,
        help="periodic boundary conditions")

    # FORCES
    # area force
    parser.add_argument("-area",
        action=BooleanOptionalAction,
        help="add area force")
    # perimeter force
    parser.add_argument("-perimeter",
        action=BooleanOptionalAction,
        help="add perimeter force")
    parser.add_argument("-p0", type=float, default=3.81,
        help="dimensionless target perimeter of cell")
    # active forces
    parser.add_argument("-taup", type=float, default=1e0,
        help="persistence time of active force")
    # active Brownian force
    parser.add_argument("-abp",
        action=BooleanOptionalAction,
        help="add active Brownian force")
    parser.add_argument("-v0", type=float, default=5e-1,
        help="vertex self-propulsion velocity")
    # Ornstein-Uhlenbeck tension
    parser.add_argument("-out",
        action=BooleanOptionalAction,
        help="add perimeter force")
    parser.add_argument("-t0", type=float, default=1,
        help="active tension mean")
    parser.add_argument("-st0", type=float, default=5e-1,
        help="active tension standard deviation")
    # MODELS 0-4
    parser.add_argument("-Gamma", type=float, default=1,
        help="perimeter elasticity constant")
    parser.add_argument("-sigma", type=float, default=1e-1,
        help="noise amplitude")
    # model 0
    parser.add_argument("-model0",
        action=BooleanOptionalAction,
        help="add model 0")

    # ALGORITHM
    parser.add_argument("-seed", type=int, default=0,
        help="random number generator seed")

    # INTEGRATION
    parser.add_argument("-dt", type=float, default=1e-2,
        help="intergation time step")
    parser.add_argument("-delta", type=float, default=0.1,
        help="length below which to perform T1")
    parser.add_argument("-epsilon", type=float, default=0.1,
        help="create junction with length epsilon above threshold after T1")
    if script == "run.py":
        parser.add_argument("-iterations", type=int, default=100,
            help="number of iterations between plots")

    # SAVING
    if script == "vm.py":
        parser.add_argument("-init", type=int, default=0,
            help="number of initial iterations")
        parser.add_argument("-niter", type=int, default=1000,
            help="number of production iterations")
        parser.add_argument("-dtmin", type=int, default=1,
            help=
                "[--no-linear-frames (default)] lag time between each frame "
                "[-linear-frames] min. lag time")
        parser.add_argument("-dtmax", type=int, default=500,
            help="[--no-linear-frames (default)] max. lag time")
        parser.add_argument("-nmax", type=int, default=50,
            help="[--no-linear-frames (default)] max. number of lag times")
        parser.add_argument("-intmax", type=int, default=20,
            help="[--no-linear-frames (default)] max. number of initial times")
        parser.add_argument("--linear-frames", "-linear-frames", "-linear",
            action=BooleanOptionalAction,
            help="compute logarithmically spaced frames")
        parser.add_argument("-id", type=int, default=0,
            help="numerical identifier for simulation file")

    # PARSE
    args = parser.parse_args()
    assert(not(args.abp and args.out))  # competing active forces

    return args

