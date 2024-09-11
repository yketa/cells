"""
Define functions to initialise a vertex model.
"""

from cells import __path__
from cells import read
from cells.bind import VertexModel

import pickle

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter,\
    BooleanOptionalAction

import numpy as np

import os
import sys

K = 1                                                   # area elasticity
A0 = 1                                                  # area of a regular hexagon
script = os.path.basename(sys.argv[0])                  # name of invoking script
out_fname = "out.p"                                     # default saving file name
movie_sh_fname = os.path.join(__path__[-1], "movie.sh") # movie making shell script file name

def init_vm(user_args=None, parser=None, **kwargs):
    """
    Initilise VertexModel object with command line arguments.

    Parameters
    ----------
    user_args : list of str or None
        Override command line arguments with these arguments when not None.
        (default: None) (see argparse.ArgumentParser().parse_args)
    parser : argparse.ArgumentParser or None
        Argument parser. (default: None)

    Additional keyword arguments are passed to the grid initialisation
    function.

    Returns
    -------
    args : argparse.Namespace
        Arguments.
    vm : VertexModel
        New instance from command line arguments.
    """

    args = parse_args(user_args, parser=parser) # command line arguments

    # INITIALISATION

    if args.input is None:

        # regular grid

        vm = VertexModel(args.seed)

        if args.periodic:
            vm.initRegularTriangularLattice(size=args.N, hexagonArea=A0,
                **kwargs)
        else:
            vm.initOpenRegularHexagonalLattice(nCells=args.N, hexagonArea=A0,
                **kwargs)
#             vm.initOpenRegularTriangularLattice(size=args.N, hexagonArea=A0)

    else:

        # input file

        try:                            # simulation file

            r = read.Read(args.input)
            vm = r[r.frames[-1] if args.frame < 0 else args.frame]
            del r

        except AssertionError:          # pickle of VertexModel

            with open(args.input, "rb") as dump:
                vm = pickle.load(dump)
            assert(type(vm) == VertexModel)

        vm = VertexModel(vm, args.seed) # change random seed

    # FORCES

    if args.perimeter:
        vm.addPerimeterForce("perimeter",
            args.Gamma, args.p0*np.sqrt(A0))
    if args.area:
        vm.addAreaForce("area",
            K, A0)
    if args.volume:
        vm.addVolumeForce("volume",
            K, args.h0*np.sqrt(A0), A0)
    if args.linear_volume:
        vm.addLinearVolumeForce("linear_volume",
            K, A0, args.Gamma, args.p0*np.sqrt(A0),
            args.taur, args.h0*np.sqrt(A0), args.taua)
    if args.gamma != 0:
        vm.addBoundaryTension("boundary_tension",
            args.gamma)
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

    # INTEGRATORS

    if args.eta != 0:
        vm.setPairFrictionIntegrator(args.eta)

    # initialise velocities and forces
    # this should in principle be called each time a force is added or deleted
    if args.input is None: vm.nintegrate(1, 0)

    return args, vm

def parse_args(user_args=None, parser=None):
    """
    Parse command line arguments for vertex model.

    Parameters
    ----------
    user_args : list of str or None
        Override command line arguments with these arguments when not None.
        (default: None) (see argparse.ArgumentParser().parse_args)
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
        help="input file name (! discards grid args)")
    parser.add_argument("-frame", "-input-frame", "-if", type=int, default=-1,
        help="input frame (< 0 gets last frame) (! used when -input != None)")

    # FORCES
    # perimeter force
    parser.add_argument("-perimeter",
        action=BooleanOptionalAction,
        help="add perimeter force")
    parser.add_argument("-Gamma", type=float, default=K*A0,
        help="junction/perimeter elasticity constant")
    parser.add_argument("-p0", type=float, default=3.81,
        help="dimensionless target perimeter of cell")
    # area force
    parser.add_argument("-area",
        action=BooleanOptionalAction,
        help="add area force")
    # volume force
    parser.add_argument("-volume",
        action=BooleanOptionalAction,
        help="add volume force")
    parser.add_argument("-h0", type=float, default=1,
        help="dimensionless target area-to-volume ratio of cell")
    # linear volume force
    parser.add_argument("-linear_volume",
        action=BooleanOptionalAction,
        help="add linear volume force")
    parser.add_argument("-taua", type=float, default=1,
        help="alignment time")
    # boundary line tension
    parser.add_argument("-gamma", type=float, default=0,
        help="line tension on open boundary")
    # active forces
    parser.add_argument("-taup", type=float, default=1e0,
        help="persistence time of active force")
    # active Brownian force
    parser.add_argument("-abp",
        action=BooleanOptionalAction,
        help="add active Brownian force")
    parser.add_argument("-v0", type=float, default=2e-1,
        help="vertex self-propulsion velocity")
    # Ornstein-Uhlenbeck tension
    parser.add_argument("-out",
        action=BooleanOptionalAction,
        help="add Ornstein-Uhlenbeck tension")
    parser.add_argument("-t0", type=float, default=4e-1,
        help="active tension mean")
    parser.add_argument("-st0", type=float, default=2e-1,
        help="active tension standard deviation")
    # MODELS 0-4
    parser.add_argument("-sigma", type=float, default=4e-2,
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

    # INTEGRATORS
    parser.add_argument("-eta", type=float, default=0,
        help="ratio of vertex-vertex friction over vertex-substrate friction")

    # ALGORITHM
    parser.add_argument("-seed", type=int, default=0,
        help="random number generator seed")

    # INTEGRATION
    parser.add_argument("-dt", type=float, default=1e-2,
        help="intergation time step")
    parser.add_argument("-delta", type=float, default=0.02,
        help="length below which to perform T1")
    parser.add_argument("-epsilon", type=float, default=0.002,
        help="create junction with length epsilon above threshold after T1")
    if not(script == "vm.py"):
        parser.add_argument("-iterations", type=int, default=100,
            help="number of iterations between plots")

    # DISPLAY
    if not(script == "vm.py"):
        parser.add_argument("-rainbow", action=BooleanOptionalAction,
            help="display rainbow plot")
        parser.add_argument("-velocities", action=BooleanOptionalAction,
            help="display velocities on vertices")
        parser.add_argument("-neighbours", action=BooleanOptionalAction,
            help="display number of neighbours on cells")

    # SAVING
    if not(script == "vm.py"):
        parser.add_argument("-movie", "-m",
            action=BooleanOptionalAction,
            help="make movie from displayed frames (see %s)" % movie_sh_fname)
#         parser.add_argument("-ffmpeg", type=str, default="",
#             help="FFmpeg executable")
        parser.add_argument("-save", "-s",
            action=BooleanOptionalAction,
            help="pickle last system state to %s" % out_fname)
    if script[-5:] == "vm.py":
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
    args = parser.parse_args(args=user_args)
    assert(not(args.abp and args.out))  # competing active forces

    return args

def get_perimeters(vm):
    """
    Return all cell perimeters.

    Parameters
    ----------
    vm : VertexModel
        Instance of vertex model.

    Returns
    -------
    perimeters : (*,) float Numpy array
        Perimeters of all cells in instance of vertex model.
    """

    return np.array(list(map(
        lambda i: vm.getVertexToNeighboursPerimeter(i),
        vm.getVertexIndicesByType("centre"))))

def get_areas(vm):
    """
    Return all cell areas.

    Parameters
    ----------
    vm : VertexModel
        Instance of vertex model.

    Returns
    -------
    areas : (*,) float Numpy array
        Areas of all cells in instance of vertex model.
    """

    return np.array(list(map(
        lambda i: vm.getVertexToNeighboursArea(i),
        vm.getVertexIndicesByType("centre"))))

