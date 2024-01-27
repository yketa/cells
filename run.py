"""
Routine to run and plot in real time a simulation of the vertex model. This
does not save data.
"""

from cells.init import init_vm
from cells.read import plot
from cells import __path__

import matplotlib.pyplot as plt

import subprocess
import os
import sys
import signal
from tempfile import TemporaryDirectory
import atexit
import traceback

if __name__ == "__main__":

    def exit_handler(*_args, **_kwargs):
        # make movie on exit
        if "args" in globals() and args.movie:
            subprocess.call([os.path.join(__path__[0], "movie.sh"),
                "-d", tmpdir.name, "-p", sys.executable, # "-f", args.ffmpeg,
                "-y"])
            tmpdir.cleanup()
        # exit
        os._exit(0)
    signal.signal(signal.SIGINT, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    atexit.register(exit_handler)

    # INITIALISATION

    args, vm = init_vm()
    fig, ax = plot(vm)

    if args.movie: tmpdir = TemporaryDirectory()

    # RUN

    plt.ion()
    plt.show()
    while True:
        # save frame
        if args.movie:
            try: count += 1
            except NameError: count = 0
            fig.savefig(os.path.join(tmpdir.name, "%05d.png" % count))
        # integrate
        try:
            vm.nintegrate(args.iterations,
                dt=args.dt, delta=args.delta, epsilon=args.epsilon)
#         vm.checkMesh(["junction"])
        except:
            print(traceback.format_exc(), file=sys.stderr)  # print traceback
            exit_handler()                                  # exit with handler
        # plot
        plot(vm, fig=fig, ax=ax)

