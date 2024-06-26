"""
Routine to run and plot in real time a simulation of the vertex model. This
does not save data.
"""

from cells.init import init_vm, out_fname, movie_sh_fname
from cells.plot import plot, plot_velocities, plot_neighbours,\
    WindowClosedException

import matplotlib.pyplot as plt

from copy import deepcopy
import subprocess
import os
from shutil import rmtree
import sys
import signal
from tempfile import mkdtemp
import atexit
import traceback
import __main__

import pickle

def _exit_handler(*_args, **_kwargs):
    # make movie on exit
    if "_frames_dir" in globals():
        # make movie from frames with ffmpeg
        try:
            subprocess.call([movie_sh_fname,
                "-d", _frames_dir, "-p", sys.executable, # "-F", args.ffmpeg,
                "-y"])
        except:
            print(traceback.format_exc(), file=sys.stderr)  # print traceback
        # delete frames
        rmtree(_frames_dir)
    # close matplotlib figures
    plt.close("all")
    # exit
    if hasattr(__main__, "__file__"): os._exit(0)           # not a python console: exit

def run(args, vm, plot_function=None, **kwargs):
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
                              cells.plot.plot_neighbours if args.neighbours,
                              cells.plot.plot            otherwise.

    Additional keywords arguments are passed to plot_function.
    """

    vm0 = deepcopy(vm)  # initial vertex model object

    # plotting function
    if plot_function is None:
        assert(not(args.velocities and args.neighbours))    # not compatible
        if args.velocities:
            plot_function = plot_velocities
        elif args.neighbours:
            plot_function = plot_neighbours
        else:
            plot_function = plot
    def _plot(*_args, **_kwargs):
        return plot_function(*_args, **{**kwargs, **_kwargs},
            rainbow=vm0 if args.rainbow else None)

    # frames directory
    global _frames_dir
    if args.movie:
        _frames_dir = mkdtemp()
        print("Saving frames to temporary directory \"%s\"." % _frames_dir,
            file=sys.stderr)
    else:
        try: del _frames_dir
        except NameError: pass

    # initialise plot
    fig, ax = _plot(vm)

    # infinite loop
    plt.ion()
    plt.show()
    while True:
        # save frame
        if args.movie:
            try: count += 1
            except NameError: count = 0
            fig.savefig(os.path.join(_frames_dir, "%05d.png" % count))
        # save system state
        if args.save:
            with open(out_fname, "wb") as dump:
                pickle.dump(vm, dump)
        # integrate and plot
        try:
            vm.nintegrate(args.iterations,
                dt=args.dt, delta=args.delta, epsilon=args.epsilon)
#         vm.checkMesh(["junction"])
            _plot(vm, fig=fig, ax=ax)
        except Exception as e:
            if not(type(e) is WindowClosedException):   # print traceback
                print(traceback.format_exc(), file=sys.stderr)
            _exit_handler()                             # exit handler
            break                                       # break loop

if __name__ == "__main__":

    signal.signal(signal.SIGINT, _exit_handler)
    signal.signal(signal.SIGTERM, _exit_handler)
    atexit.register(_exit_handler)

    # INITIALISATION

    args, vm = init_vm()

    # RUN

    run(args, vm, plot_function=None)

