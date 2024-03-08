# Cells

Python library written in C++ to integrate vertex models.

![polygonal tiling](docs/cells.svg)

## Compiling and running

### Compiling

#### Directly on the machine

**PLEASE READ THIS PART ENTIRELY BEFORE EXECUTING ANY COMMAND.**

Compilation of the shared library with `make bind.so` requires `pybind11` (`python -m pip install -r requirements.txt`) and a C++20 compiler.

Python scripts are written for `python3` and import the `cells` package which necessitates the directory containing this repository to be added to the `$PYTHONPATH`. This can be achieved e.g. by executing *from the directory containing this README file* the following commands
```sh
( [[ "${PWD##*/}" == "cells" ]] && echo "export PYTHONPATH=\$PYTHONPATH:${PWD}/.." >> ~/.${0}rc && source ~/.${0}rc && echo "Success." ) || echo "Error: this directory is not 'cells'."
```
which output `Success.` if succesful.

Python routines `read.py`, `run.py` and `vm.py` require `numpy` and `matplotlib` (`python -m pip install -r requirements.txt`).

#### In a container

It is possible to compile this library in a `singularity` container (with `sudo` privilege) with `make container.sif`. Package `cells` is then available through the `python` interpreter of the container with `./container.sif python -m cells` or `singularity exec container.sif python -m cells`.

### Routines and modules

There are two default routines to simulate vertex models: `run.py` runs and plots in real time a simulation of a vertex model, and `vm.py` runs and saves a simulation of a vertex model. These routines rely on modules `init.py`, `plot.py` and `read.py`.

Module `init.py` defines functions to parse command line arguments and initialise vertex models as function of these arguments. A list of these arguments can be displayed with `python run.py -h` (respectively `python vm.py -h`) or `python -m cells.run -h` (respectively `python -m cells.vm -h`).

Module `plot.py` defines functions to plot vertex models.

Module `read.py` defines objects and functions to access and plot vertex model data. Executed as a routine, with a simulation file name as a command line argument, this prints `true` (respectively `false`) if the file is consistent (respectively not consistent).

### Additional scripts

Script `movie.sh` is a quick tool to make movies and requires `plot.py`, `read.py`, and [`ffmpeg`](https://ffmpeg.org/download.html). Calling module `run.py` with command line argument `-m` (or `-movie`) will save displayed frames and make a movie from these frames when exited.

### Examples

```bash
python -m cells.run -abp -area -perimeter
python -m cells.run -abp -area -perimeter -forces
python -m cells.run -abp -area -perimeter -forces -m
python -m cells.run -out -area -periodic -N 12
```

## C++ source files

### Vertex model

This vertex model implementation is separated in two parts. First `mesh.*pp` contains the implementation of the geometrical features of the model, which relies on vertices which are linked together by directed half-edges enabling to move accross the mesh (see `docs/mesh.pdf`). Second `system.*pp` provides a general `VertexModel` class to integrate the dynamics.

`VertexModel` methods to initialise configurations are defined in `initialisation.cpp`.

### Python binding

`VertexModel` is exposed to Python through `pybind11` in `bind.cpp`. Importantly, this file also contains additional methods to provide control of the forces and faster access to information within `VertexModel` in Python.

### Forces

Forces are defined separately and are attached to the `VertexModel` class following a class factory method (see `class_factory.hpp`) which enables to add and remove different forces. There are two general classes of forces enabled (see `base_forces.hpp`): forces which are computed for all vertices of a given type (`VertexForce` class), and forces which are computed for all half-edges of a given type (`HalfEdgeForce` class).

Definitions of forces belong in `forces.hpp`. Script `forces.cpp` provides definitions of `VertexModel` methods to add them to the simulation object. Header `forces_pickle.hpp` provides definitions of own and `VertexModel` methods to pickle and unpickle these forces and their state. Finally, forces are exposed to Python in `bind.cpp`. On the Python side, forces are initialised in `init.py` and forces-dependent plotting are defined in `plot.py`.

Descriptions of some forces can be found in `docs/forces.pdf` and `docs/active_junction.pdf`.

### Pickling

`VertexModel` instances including the associated forces can be pickled via Python (see `forces_pickle.hpp` and `system_pickle.hpp`).

## Authors

- [Yann-Edwin Keta](keta@lorentz.leidenuniv.nl), Instituut-Lorentz, Universiteit Leiden
- [Silke Henkes](shenkes@lorentz.leidenuniv.nl), Instituut-Lorentz, Universiteit Leiden

