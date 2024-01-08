"""
Run and plot in real time a simulation of the vertex model. This does not save
data.
"""

from cells.init import init_vm
from cells.vm import plot

import matplotlib.pyplot as plt

if __name__ == "__main__":

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

