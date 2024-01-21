"""
Routine to run and save a simulation of the vertex model.
"""

from cells.init import init_vm
from cells.exponents import float_to_letters

from datetime import datetime
import atexit, signal
import os

import numpy as np
from collections import OrderedDict

import pickle

def filename(N, identifier, prefix=None):
    """
    Standard filename for simulation file.

    Parameters
    ----------
    N : int
        Number of vertices.
    identifier : int
        Unique integer identifier for file.
    prefix : str or None
        Prefix to the simulation file name. (default: None)

    Returns
    -------
    name : str
        File name.
    """

    return (("" if prefix is None else "%s_" % prefix) + "N%s_%%%s.p"
        % tuple(map(float_to_letters, (N, identifier))))

if __name__ == "__main__":

    # COMPUTATION TIME

    start_t = datetime.now()
    print("Started on %s." % start_t,
        flush=True)

    def exit_handler(*_args, **_kwargs):
        end_t = datetime.now()
        print("Stopped on %s (elapsed: %s)." % (end_t, end_t - start_t),
            flush=True)
        os._exit(0)
    # print elapsed time on exit
    signal.signal(signal.SIGINT, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    atexit.register(exit_handler)

    # INITIALISATION

    args, vm = init_vm()

    # METADATA

    metadata = {        # metadata for simulation file
        "filename": filename(len(vm.vertices), args.id, prefix=args.filename),
        "dt": args.dt,  # used to check computed times
        "args": args,   # save all arguments
    }

    # CHOOSE FRAMES

    if args.linear_frames:  # LINEARLY SPACED FRAMES
        if args.niter%args.dtmin:
            raise ValueError(
                "No integer multiple of %i frames in a total of %i frames."
                    % (args.dtmin, args.niter))
        else:
            metadata["t0"] = np.array(                                      # initial times
                [args.init                                                  # initialisation frames
                    + i*args.dtmin for i in range(args.niter//args.dtmin)], # frames spaced with dtmin
                dtype=int)
            metadata["t"] = np.array([0, args.dtmin])                       # lag times
    else:                       # LOGARITHMICALLY SPACED FRAMES
        dtmax = np.min([args.dtmax, args.niter - (args.intmax - 1)])        # maximum lag time
        metadata["t0"] = np.array(                                          # initial times
            [args.init                                                      # initialisation frames
                + i*(args.niter - dtmax)/(args.intmax - 1)
                    for i in range(args.intmax)],
            dtype=int)
        metadata["t"] = np.array(                                           # lag times
            [args.dtmin]                                                    # minimum lag time
                + [args.dtmin
                    + np.exp(i*np.log(dtmax - args.dtmin)/(args.nmax - 1))  # frames logarithmically spaced
                        for i in range(1, args.nmax - 1)]
                + [dtmax],                                                  # maximum lag time
            dtype=int)

    metadata["frames"] = np.array(list(OrderedDict().fromkeys(sorted(       # ensemble of unique frames
        [0, *[t0 + t
            for t0 in metadata["t0"]
            for t in [0, *metadata["t"]]]]))))

    # DUMP METADATA

    with open(metadata["filename"], "wb") as dump:
        pickle.dump(metadata, dump)
    print("Writing to \"%s\"." % metadata["filename"])

    # SIMULATION

    for t in np.diff(metadata["frames"], prepend=0):
        vm.nintegrate(t, args.dt, args.delta, args.epsilon)
        vm.checkMesh(["junction"])
        with open(metadata["filename"], "ab") as dump:
            pickle.dump(vm, dump)

