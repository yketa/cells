"""
Routine to run and save a simulation of the vertex model.
"""

from cells.init import init_vm
from cells.read import ReadWYC
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
        # print elapsed time on exit
        end_t = datetime.now()
        print("Stopped on %s (elapsed: %s)." % (end_t, end_t - start_t),
            flush=True)
        # exit
        os._exit(0)
    signal.signal(signal.SIGINT, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    atexit.register(exit_handler)

    # INITIALISATION

    args, vm = init_vm()

    # --- resume from last frame of incomplete file
    if args.resume:

        # METADATA

        r = ReadWYC(args.input)             # read
        metadata = r.metadata               # copy metadata
        metadata["filename"] = r.filename   # change filename

        # VERTEX MODEL OBJECT

        assert r.skip.size > 0                          # check there is at least one frame
        assert r.skip.size < r.frames.size              # check there remains frame to compute
        last_saved_frame = r.frames[r.skip.size - 1]    # last frame in the file
        vm = r[last_saved_frame]                        # vertex model object corresponding to frame

    # --- start new simulation
    else:

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
        else:                   # LOGARITHMICALLY SPACED FRAMES
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
        print("Writing to \"%s\"." % metadata["filename"],
            flush=True)

    # SIMULATION

    frames = metadata["frames"]
    if args.resume: # --- resume from last frame of incomplete file
        time_increments = np.diff(frames[frames >= last_saved_frame])
    else:           # --- start new simulation
        time_increments = np.diff(frames, prepend=0)
    for t in time_increments:
        vm.nintegrate(t, args.dt, args.delta, args.epsilon)
#         vm.checkMesh(["junction"])
        with open(metadata["filename"], "ab") as dump:
            pickle.dump(vm, dump)

