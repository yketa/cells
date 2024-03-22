"""
Routine to run and plot in real time a simulation of the vertex model. This
does not save data.
"""

from cells.init import init_vm, out_fname, movie_sh_fname
from cells.plot import plot, plot_velocities

import matplotlib.pyplot as plt

import subprocess
import os
from shutil import rmtree
import sys
import signal
from tempfile import mkdtemp
import atexit
import traceback

import pickle

def movie_on_exit():
    """
    Cleanup function to make movie from frames in `frames_dir` at script exit.
    """

    def exit_handler(*_args, **_kwargs):
        # make movie on exit
        if "frames_dir" in globals():
            try:
                subprocess.call([movie_sh_fname,
                    "-d", frames_dir, "-p", sys.executable, # "-F", args.ffmpeg,
                    "-y"])
            except:
                print(traceback.format_exc(), file=sys.stderr)  # print traceback
            rmtree(frames_dir)
        # exit
        os._exit(0)

    signal.signal(signal.SIGINT, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    atexit.register(exit_handler)

def run(args, vm, plot_function=None):
    """
    Infinite loop to run and plot vertex model simulation.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments. (see cells.init)
    vm : cells.bind.VertexModel
        Vertex model object.
    plot_function : function
        Plotting function. (see cells.plot, default: None)
        NOTE: if plot_function is None then
              plot_function = cells.plot.plot_velocities if args.velocities,
                              cells.plot.plot            otherwise.
    """

    if plot_function is None:
        plot_function = plot_velocities if args.velocities else plot

    # frames directory
    if args.movie:
        global frames_dir
        frames_dir = mkdtemp()

    # initialise plot
    fig, ax = plot_function(vm)

    # infinite loop
    plt.ion()
    plt.show()
    while True:
        # save frame
        if args.movie:
            try: count += 1
            except NameError: count = 0
            fig.savefig(os.path.join(frames_dir, "%05d.png" % count))
        # save system state
        if args.save:
            with open(out_fname, "wb") as dump:
                pickle.dump(vm, dump)
        # integrate
        try:
            vm.nintegrate(args.iterations,
                dt=args.dt, delta=args.delta, epsilon=args.epsilon)
#         vm.checkMesh(["junction"])
        except:
            print(traceback.format_exc(), file=sys.stderr)  # print traceback
            exit(0)                                         # exit with handler
        # plot
        plot_function(vm, fig=fig, ax=ax)

if __name__ == "__main__":

    movie_on_exit()

    # INITIALISATION

    args, vm = init_vm()

    # RUN

    run(args, vm, plot_function=None)

