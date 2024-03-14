"""
Define functions to initialise a vertex model.
"""

from cells.bind import VertexModel
from cells.read import Read
from cells import __path__

import pickle

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,\
    BooleanOptionalAction

import numpy as np

import os
import sys

A0 = (3./2.)/np.tan(np.pi/6.)           # area of a regular hexagon with edge length 1
script = os.path.basename(sys.argv[0])  # name of invoking script

def init_vm(parser=None):
    """
    Initilise VertexModel object with command line arguments.

    Parameters
    ----------
    parser : argparse.ArgumentParser or None
        Argument parser. (default: None)

    Returns
    -------
    args : argparse.Namespace
        Arguments.
    vm : VertexModel
        New instance from command line arguments.
    """

    args = parse_args(parser=parser)    # command line arguments

    # INITIALISATION

    if args.input is None:

        # regular grid

        vm = VertexModel(args.seed)

        if args.periodic:
            vm.initRegularTriangularLattice(size=args.N)
        else:
            vm.initOpenRegularHexagonalLattice(nCells=args.N)
#             vm.initOpenRegularTriangularLattice(size=args.N)

    else:

        # input file

        try:                    # simulation file

            r = Read(args.input)
            vm = r[r.frames[-1] if args.frame < 0 else args.frame]
            del r

        except AssertionError:  # pickle of VertexModel

            with open(args.input, "rb") as dump:
                vm = pickle.load(dump)
            assert(type(vm) == VertexModel)

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
    if args.model1:
        vm.addModel1("model1",
            args.Gamma, args.p0*np.sqrt(A0), args.sigma, args.taup)
    if args.model2:
        vm.addModel2("model2",
            args.Gamma, args.taur, args.sigma, args.taup)
    if args.model3:
        vm.addModel3("model3",
            args.Gamma, args.sigma, args.taup)
    if args.model4:
        vm.addModel4("model4",
            args.Gamma, args.taur, args.sigma, args.taup)

    # initialise velocities and forces
    # this should in principle be called each time a force is added or deleted
    vm.nintegrate(1, 0)

    return args, vm

def parse_args(parser=None):
    """
    Parse command line arguments for vertex model.

    Parameters
    ----------
    parser : argparse.ArgumentParser or None
        Argument parser. (default: None)

    Returns
    -------
    args : argparse.Namespace
        Arguments.
    """

    if parser is None:
        parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    # INITIALISATION
    # regular grid
    parser.add_argument("-N", type=int, default=36,
        help=
            "[--no-periodic (default)] number of cells (! square number) "
            "[-periodic] number of vertices in each direction"
            "(! multiple of 6)")
    parser.add_argument("--periodic", "-periodic", "-p",
        action=BooleanOptionalAction,
        help="periodic boundary conditions")
    # input file
    parser.add_argument("-input", "-input-name", "-i", type=str, default=None,
        help="input file name (! discards grid args and seed if != None)")
    parser.add_argument("-frame", "-input-frame", "-if", type=int, default=0,
        help="input frame (< 0 gets last frame) (! used when -input != None)")

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
        help="junction/perimeter elasticity constant")
    parser.add_argument("-sigma", type=float, default=1e-1,
        help="noise amplitude")
    parser.add_argument("-taur", type=float, default=1e0,
        help="relaxation time")
    # model 0
    parser.add_argument("-model0",
        action=BooleanOptionalAction,
        help="add model 0")
    # model 1
    parser.add_argument("-model1",
        action=BooleanOptionalAction,
        help="add model 1")
    # model 2
    parser.add_argument("-model2",
        action=BooleanOptionalAction,
        help="add model 2")
    # model 3
    parser.add_argument("-model3",
        action=BooleanOptionalAction,
        help="add model 3")
    # model 2
    parser.add_argument("-model4",
        action=BooleanOptionalAction,
        help="add model 4")

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
    if not(script == "vm.py"):
        parser.add_argument("-iterations", type=int, default=100,
            help="number of iterations between plots")

    # DISPLAY
    if not(script == "vm.py"):
        parser.add_argument("-velocities", action=BooleanOptionalAction,
            help="display velocities on vertices")

    # SAVING
    if not(script == "vm.py"):
        parser.add_argument("-movie", "-m",
            action=BooleanOptionalAction,
            help="make movie from displayed frames (see %s)"
                % os.path.join(__path__[0], "movie.sh"))
#         parser.add_argument("-ffmpeg", type=str, default="",
#             help="FFmpeg executable")
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
            help="compute linearly spaced frames")
        parser.add_argument("-filename", "-f", type=str, default=None,
            help="prefix to the simulation file name")
        parser.add_argument("-id", type=int, default=0,
            help="numerical identifier for simulation file")

    # PARSE
    args = parser.parse_args()
    assert(not(args.abp and args.out))  # competing active forces

    return args

