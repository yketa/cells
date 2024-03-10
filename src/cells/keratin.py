"""
Routine to run and plot in real time a simulation of the keratin vertex model.
This does not save data.
"""

import sys
sys.argv[0] = "run.py"  # appear as "run.py"

from cells.init import init_vm, A0
from cells.plot import plot
from cells import __path__

import matplotlib.pyplot as plt

import subprocess
import os
import signal
from tempfile import TemporaryDirectory
import atexit
import traceback

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np

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

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    # K is 1 by default
    # A0 is (3./2.)/np.tan(np.pi/6.) by default
    # Gamma is 1 by default
    # P0 is defined by sqrt(A0)*{p0}
    parser.add_argument("-l0", type=float, default=1,
        help="bond rest length")
    parser.add_argument("-alpha", type=float, default=1,
        help="bond elasticity per keratin concentration above threshold")
    parser.add_argument("-kth", type=float, default=1,
        help="bond keratin concentration threshold")
    # tau is defined by {taur}
    # sigma is defined by {sigma}
    parser.add_argument("-tauon", type=float, default=1,
        help="keratin concentration on-rate evolution time scale")
    parser.add_argument("-k0", type=float, default=1,
        help="keratin concentration off-rate inverse pressure constant")
    parser.add_argument("-pr0", type=float, default=1,
        help="keratin concentration off-rate inflection pressure")

    args, vm = init_vm(parser)
    fig, ax = plot(vm)

    if args.movie: tmpdir = TemporaryDirectory()

    # KERATIN

    assert len(vm.vertexForces) == 0
    assert len(vm.halfEdgeForces) == 0
    vm.addKeratinModel("keratin",
        1, A0, 1, args.p0*np.sqrt(A0), args.l0, args.alpha, args.kth,
        args.taur, args.sigma, args.tauon, args.k0, args.pr0)

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

