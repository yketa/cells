"""
Define objects and functions to read data file and plot vertex model object.
When executed this checks the consistency of a simulation file.
"""

from cells.bind import VertexModel
from cells.plot import WindowClosedException, plot

import os
import sys
from tempfile import gettempdir
import traceback

from math import ceil
import numpy as np

import matplotlib.pyplot as plt

import pickle

# TOOLS

# progress bar (https://stackoverflow.com/a/3160819/7385044)
_toolbar_width = 40
def _progressbar(p):
    out = "[%s] (%2i%%) " % (
        "="*ceil(p*_toolbar_width)
            + " "*(_toolbar_width - ceil(p*_toolbar_width)),
        ceil(100*p))
    sys.stderr.write(out)
    sys.stderr.flush()
    sys.stderr.write("\b"*len(out))
    if p == 1: sys.stderr.write("\n")   # end the progress bar

# READ SIMULATION

class Read:
    """
    Simulation file reading object.
    """

    def __init__(self, fname, check=None):
        """
        Extract metadata and check consistency of file.

        Parameters
        ----------
        fname : str
            Path to input file.
        check : None or bool
            Check file consistency. This check must be done at least once in
            order to access data, otherwise only metadata is read and
            accessible. (default: None)
            NOTE: After a succesful check, a corresponding data file is saved
                  in a temporary directory for later faster loads. If check is
                  None then this file is used if it exists, otherwise or if
                  check is True then it is computed and saved. If check is
                  False then it is not checked at all.
            NOTE: Use cells.read.ReadWYC to access as much data as possible in
                  the case of corrupted files.
        """

        self.filename = fname
        self.fig, self.ax = None, None  # used for plotting

        with open(self.filename, "rb") as dump:
            self.metadata = pickle.load(dump)
            assert type(self.metadata) == dict
            current_pointer = dump.tell()   # current position of the read pointer

            # default metadata
            self.t0 = self.metadata["t0"]           # initial times
            self.t = self.metadata["t"]             # lag times
            self.frames = self.metadata["frames"]   # array of computed frames
            self.dt = self.metadata["dt"]           # integration time step
            if check is False: return

            # saved skipped dictionary
            data_fname = os.path.join(
                gettempdir(), "%s.data" % os.path.basename(self.filename))
            if check is None:
                try:
                    with open(data_fname, "rb") as data_dump:
                        self.skip = pickle.load(data_dump)
                    print("Data loaded from \"%s\"." % data_fname,
                        file=sys.stderr)
                    return
                except: pass

            # check file consistency and build skip directory
            self.skip = np.array([], dtype=int)     # position of each object in file
            _max_diff_t = 0                         # check consistency with metadata in computed times
            for i, time in enumerate(self.frames*self.dt):
                _progressbar(i/self.frames.size)    # display progress bar
                vm = pickle.load(dump)              # load vertex model object
                assert(type(vm) == VertexModel)     # check object has correct type
                if time == 0:                       # initial time (might be != 0 if started from input)
                    time0 = vm.time
                else:                               # compute relative difference in time
                    _max_diff_t = (
                        max(_max_diff_t, np.abs(time0 + time - vm.time)/time))
                self.skip = np.append(self.skip, current_pointer)
                current_pointer = dump.tell()       # current position of the read pointer
            _progressbar(1)
            try:
                pickle.load(dump)   # this should raise an EOFError if the file was read completely
                raise ValueError("File size is not consistent with metadata.")
            except EOFError:
                pass

        assert _max_diff_t < 1e-8   # relative difference between python and C++ times

        # save skipped dictionary
        with open(data_fname, "wb") as dump:
            pickle.dump(self.skip, dump)
        print("Data saved to \"%s\"." % data_fname,
            file=sys.stderr)

    def plot(self, frame, rainbow=None, override=None, **kwargs):
        """
        Plot vertex model state corresponding to frame.

        Parameters
        ----------
        frame : int
            Index of frame.
        rainbow : int or None
            Index of previous frame of the system with respect to which colour
            cells. (default: None)
            NOTE: if rainbow != None then this overrides all other cell
                  colouring.
        override : function or None
            Override plotting function. (default: None)
            NOTE: This is expected to return a figure and an axes subplot, and
                  to support keyword argument `rainbow'. (see cells.plot.plot)

        Additional keywords arguments are passed to the plotting function.
        """

        self.fig, self.ax = (
            # plotting function
            override if not(override is None) else plot)(
            # plotting argument
            self[frame], fig=self.fig, ax=self.ax,
            rainbow=rainbow if rainbow is None else self[rainbow],
            **kwargs)

    def play(self, frames=None, loop=True, **kwargs):
        """
        Plot frames.

        Parameters
        ----------
        frames : int array-like or None
            Frames to plot.
            NOTE: if frames == None then all frames in self.frames are plotted.
        loop : bool
            Loop frames. (default: True)

        Additional keywords arguments are passed to self.plot.
        """

        if frames is None:
            frames = self.frames
        frames = list(frames)

        _progressbar(0)
        self.plot(frames[0], **kwargs)
        if not(plt.isinteractive()):
            plt.ion()
            plt.show()
        for frame in frames[1:]:
            _progressbar(frames.index(frame)/len(frames))
            try:
                self.plot(frame, **kwargs)
            except WindowClosedException:       # when window is closed
                plt.close(self.fig)             # close figure
                self.fig, self.ax = None, None  # and reset figure and axis
                _progressbar(1)
                return
        if loop:
            self.play(frames=frames, loop=loop, **kwargs)   # play again
        else:
            self.fig, self.ax = None, None      # reset figure and axis
            _progressbar(1)
            return

    def __getitem__(self, frame):
        """
        Returns saved vertex model state corresponding to frame.

        Parameters
        ----------
        frame : int
            Index of frame.

        Returns
        -------
        vm : cells.bind.VertexModel
            Vertex model state.
        """

        assert frame in self
        with open(self.filename, "rb") as dump:
            dump.seek(self.skip[self.frames.tolist().index(frame)]) # skip until object (https://stackoverflow.com/questions/76252112)
            vm = pickle.load(dump)
            assert type(vm) == VertexModel
            if False:
                # initialise velocities and forces
                vm.nintegrate(1, 0)

        return vm

    def __contains__(self, frame):  # is frame in saved frames?
        return self.frames.__contains__(frame)

    def __iter__(self):             # iterator over frames
        return self.frames.__iter__()

    def __next__(self):             # iterator over frames
        return self.frames.__next__()

    def __len__(self):              # number of frames
        return self.frames.__len__()

    def __del__(self):
        if not(self.fig is None):
            plt.close(self.fig)
            del self.fig

# READ CORRUPTED SIMULATION (experimental)

class ReadWYC(Read):
    """
    CAUTION: This object is expected to be prone to errors.

    Derived class of simulation file reading object cells.read.Read which
    bypasses exceptions due to file corruptions in object initialiser.
    https://stackoverflow.com/questions/51254859/
    """

    def __init__(self, *args, **kwargs):
        """
        CAUTION: This object is expected to be prone to errors.

        All arguments are passed to cells.read.Read.__init__.

        ==== [cells.read.Read.__init__.__doc__] ====

        {0}
        """

        try: super().__init__(*args, **kwargs)
        except: print(traceback.format_exc(), file=sys.stderr)  # print traceback

ReadWYC.__init__.__doc__ = (    # append docstring https://stackoverflow.com/questions/58256881/
    ReadWYC.__init__.__doc__.format(Read.__init__.__doc__))

# CHECK

if __name__ == "__main__":

    if len(sys.argv) == 2:
        try:
            Read(sys.argv[1], check=None)
            print("true")   # print lowercase to be treated as bash boolean
        except:
            print("false")  # print lowercase to be treated as bash boolean

