# Cells

Python3 library written in C++20 to integrate vertex models.

![polygonal tiling](docs/cells.svg)

This readme assumes that `python` is the Python3 command.

## Installing

**PLEASE READ THIS SECTION ENTIRELY BEFORE EXECUTING ANY COMMAND.**

### Using `pip`

This requires a C++20 compiler and [`eigen3`](https://eigen.tuxfamily.org/index.php).

Configuration file `pyproject.toml` defines instructions to build the libraries and install the `cells` package using `pip`. This is achieved by executing *from the directory containing this readme file* the following command
```sh
python -m pip install . --verbose --break-system-packages
```
then -- after successful completion -- package information should be given by `python -m pip show cells`.

Package `cells` can be uninstalled by running `python -m pip uninstall cells --break-system-packages`.

Package can be installed in [editable mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html)
```sh
python -m pip install --editable . --verbose --break-system-packages
```
which necessitates however that the libraries be manually compiled (see section "Compiling manually [...]").

### Compiling manually and adding to `$PYTHONPATH`

This requires a C++20 compiler and [`eigen3`](https://eigen.tuxfamily.org/index.php), as well as `pybind11` (`python -m pip install -r requirements.txt --break-system-packages`).

Makefile `src/cells/Makefile` defines instructions to build the libraries. Compilation is achieved by executing *from the source code directory (`src/cells`)* the command `make`.

Addition of the `cells` package to the `$PYTHONPATH` (i.e. making it discoverable to `python`) is achieved by executing *from the directory containing this readme file* the following commands
```sh
echo "export PYTHONPATH=\$PYTHONPATH:${PWD}/src" >> ~/.${0}rc   # append to shell configuration file
source ~/.${0}rc                                                # reload shell configuration file
```
then -- after successful completion -- the module path should be given by `python -c "import cells; print(cells.__path__)"`.

Addition to the `$PYTHONPATH` in the shell configuration file can be reverted by executing *from the directory containing this readme file* the following command
```sh
sed -i '/'"export PYTHONPATH=\$PYTHONPATH:${PWD//'/'/'\/'}"'\/src/d' ~/.${0}rc
```
even though it is preferable to perform this operation manually (e.g. `vim + ~/.${0}rc`).

### Building a `singularity` container

This requires [`singularity`](https://docs.sylabs.io/guides/latest/user-guide/) and admin (`sudo`) privilege.

Configuration file `container.def` and makefile `container.makefile` define instructions to build the container. This is achieved by executing *from the directory containing this readme file* the following command
```sh
make -f container.makefile
```
then -- after successful completion -- package information should be given by `./container.sif python -m pip show cells`.

This operation is reverted by `rm -iv container.sif`.

Note that package `cells` is then available through the `python` interpreter of the container with `./container.sif python -m cells` or `singularity exec container.sif python -m cells`.

### Installing `eigen3`

Library [`eigen3`](https://eigen.tuxfamily.org/index.php) should preferably be installed through your distribution package manager (e.g. `sudo apt install libeigen3-dev` on Debian Linux), or through environment modules (e.g. check it is available with `module avail eigen` then load it with `module load ...`).

Makefile `src/cells/Makefile` calls `pkg-config` to append include paths which are not already part of [`$CPATH`](https://gcc.gnu.org/onlinedocs/gcc/Environment-Variables.html). As a backup for when this fails (e.g. because `pkg-config` does not exist and `Eigen` is not found in any of the directories in `$CPATH`) the makefile also adds `src/cells/extern` as an include path, therefore directory [`Eigen`](https://gitlab.com/libeigen/eigen/-/tree/master/Eigen) should be placed there (e.g. with a symbolic link).

It is possible to clone the complete library by executing *from the directory containing this readme file* the following commands
```sh
git clone https://gitlab.com/libeigen/eigen.git src/cells/extern/eigen.git  # clone repository to extern
ln -s eigen.git/Eigen src/cells/extern                                      # symlink to library
```
then -- after successful completion -- the library header files should be given by `ls src/cells/Eigen`.

## Running

After installation, use `python -c "import cells; print(cells.__path__)"` to find the location of the package files.

Python routines `read.py`, `run.py`, and `vm.py` require `numpy` and `matplotlib` (`python -m pip install -r requirements.txt --break-system-packages`), as well as the built libraries compiled in `bind.so` (see section "Installing").

### Routines and modules

There are two default routines to simulate vertex models: `run.py` runs and plots in real time a simulation of a vertex model, and `vm.py` runs and saves a simulation of a vertex model as a pickle file. These routines rely on modules `init.py`, `plot.py` and `read.py`.

Module `init.py` defines functions to parse command line arguments and initialise vertex models as function of these arguments. A list of these arguments can be displayed with `python -m cells.run -h` and `python -m cells.vm -h`.

Module `plot.py` defines functions to plot vertex models.

Module `read.py` defines objects and functions to access and plot vertex model data. Executed as a routine (`python -m cells.read`) , with a simulation file name as a command line argument, this prints `true` (respectively `false`) if the file is consistent (respectively not consistent).

Module `run.py` defines functions to run and plot simultaneously vertex model simulations.

### Additional scripts

Script `movie.sh` is a quick tool to make movies and requires `plot.py`, `read.py`, and [`ffmpeg`](https://ffmpeg.org/download.html). Calling module `run.py` with command line argument `-m` (or `-movie`) will save displayed frames and make a movie from these frames when exited.

### Examples

Vertex models real-time simulations can be run either from the shell command line
```sh
python -m cells.run -abp -area -perimeter
python -m cells.run -abp -area -perimeter -velocities
python -m cells.run -abp -area -perimeter -velocities -m
python -m cells.run -out -area -periodic -N 12
```
or from a `python` interpreter
```python
from cells.init import init_vm
from cells.run import run
args, vm = init_vm(user_args=["-abp", "-area", "-perimeter"]); run(args, vm)
args, vm = init_vm(user_args=["-abp", "-area", "-perimeter", "-velocities"]); run(args, vm)
args, vm = init_vm(user_args=["-abp", "-area", "-perimeter", "-velocities", "-m"]); run(args, vm)
args, vm = init_vm(user_args=["-out", "-area", "-periodic", "-N", "12"]); run(args, vm)
```
and in both cases the script halts when the plotting window is closed. It is possible to zoom in on parts of the figure using the `matplotlib` magnifier tool, and to reset zoom by double clicking on the figure.

## C++ source files

These are located in `src/cells`.

### Vertex model

This vertex model implementation is separated in two parts. First `mesh.*pp` contains the implementation of the geometrical features of the model, which relies on vertices which are linked together by directed half-edges enabling to move across the mesh (see `docs/mesh.pdf`). Second `system.*pp` provides a general `VertexModel` class to integrate the dynamics.

`VertexModel` methods to initialise configurations are defined in `initialisation.cpp`.

### Python binding

`VertexModel` is exposed to Python through `pybind11` in `bind.cpp`. Importantly, this file also contains additional methods to provide control of the forces and faster access to information within `VertexModel` in Python.

### Forces

Forces are defined separately and are attached to the `VertexModel` class following a class factory method (see `class_factory.hpp`) which enables to add and remove different forces. There are two general classes of forces enabled (see `base_forces.hpp`): forces which are computed for all vertices of a given type (`VertexForce` class), and forces which are computed for all half-edges of a given type (`HalfEdgeForce` class).

Definitions of forces belong in `forces.hpp`. Script `forces.cpp` provides definitions of `VertexModel` methods to add them to the simulation object. Header `forces_pickle.hpp` provides definitions of own and `VertexModel` methods to pickle and unpickle these forces and their state. Finally, forces are exposed to Python in `bind.cpp`. On the Python side, forces are initialised in `init.py` and forces-dependent plotting are defined in `plot.py`.

Descriptions of some forces can be found in `docs/forces.pdf` and `docs/active_junction.pdf`.

### Integrators

Integrators compute, from the velocities and forces at a given time, the velocities at the next time step. Header `base_integrators.hpp` provides their base class. They are defined separately and are attached to the `VertexModel` through a smart pointer which enables to switch between integrators.

Definitions of integrators belong in `integrators.hpp`. Script `integrators.cpp` provides definitions of `VertexModel` methods to set them in the simulation object. Header `integrators_pickle.hpp` provides definitions of own and `VertexModel` methods to pickle and unpickle these integrators and their state. Finally, integrators are exposed to Python in `bind.cpp`. On the Python side, integrators are initialised in `init.py`.

### Pickling

`VertexModel` instances (including forces and integrators) can be pickled via Python (see `*_pickle.hpp`).

## Authors

- [Yann-Edwin Keta](mailto:keta@lorentz.leidenuniv.nl), Instituut-Lorentz, Universiteit Leiden
- [Silke Henkes](mailto:shenkes@lorentz.leidenuniv.nl), Instituut-Lorentz, Universiteit Leiden

