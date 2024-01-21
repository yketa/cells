"""
Routine to run and plot in real time a simulation of the vertex model. This
does not save data.
"""

from cells.init import init_vm
from cells.read import plot

import matplotlib.pyplot as plt

import signal
import os

if __name__ == "__main__":

    signal.signal(signal.SIGINT, lambda _, __: os._exit(0)) # avoid traceback when interrupting

    # INITIALISATION

    args, vm = init_vm()
    fig, ax = plot(vm)

    # RUN

    plt.ion()
    plt.show()
    while True:
        # integrate
        vm.nintegrate(args.iterations,
            dt=args.dt, delta=args.delta, epsilon=args.epsilon)
#         vm.checkMesh(["junction"])
        # plot
        plot(vm, fig=fig, ax=ax)

