# Cells

Python library written in C++ to integrate vertex models.

![polygonal tiling](docs/cells.svg)

## Compiling and running

### Directly on the machine

Compilation of the shared library with `make bind.so` requires `pybind11` (`python -m pip install -r requirements.txt`) and a C++20 compiler.

Python scripts are written for `python3` and import the `cells` package which necessitates the directory containing this repository to be added to the `$PYTHONPATH`, e.g. by executing
```
echo "export PYTHONPATH=\$PYTHONPATH:${PWD}/.." >> ~/.bashrc
```
from this directory.

Python routines `run.py` and `vm.py` require `numpy` and `matplotlib` (`python -m pip install -r requirements.txt`).

### In a container

It is possible to compile this library in a `singularity` container (with `sudo` privilege) with `make container.sif`. Package `cells` is then available through the `python` interpreter of the container with `./container.sif python -m cells` or `singularity exec container.sif python -m cells`.

### Routines

There are two default routines to simulate vertex models (VMs): `run.py` plots in real time a simulation of a vertex model, and `vm.py` saves simulations of vertex models. It is noteworthy that `vm.py` also defines objects and functions to access and plot vertex model data (and on which `run.py` relies).

Both routines rely on `init.py` which parses command line arguments and initialise vertex models as function of these arguments. A list of these arguments can be displayed with `python run.py -h` (respectively `python vm.py -h`) or `python -m cells.run -h` (respectively `python -m cells.vm -h`).

### Examples

```
python -m cells.run -abp -area -perimeter
python -m cells.run -out -area -periodic -N 12
```

## C++ scripts

### Vertex model

This vertex model implementation is separated in two parts. First `mesh.*pp` contains the implementation of the geometrical features of the model, which relies on vertices which are linked together by directed half-edges enabling to move accross the mesh (see `docs/mesh.pdf`). Second `system.*pp` provides a general `VertexModel` class to integrate the dynamics.

`VertexModel` methods to initialise configurations are defined in `initialisation.cpp`.

### Python binding

`VertexModel` is exposed to Python through `pybind11` in `bind.cpp`. Importantly, this file also contains additional methods to provide control of the forces in Python.

### Forces

Forces are defined separately and are attached to the `VertexModel` class following a class factory method (see `class_factory.hpp`) which enables to add and remove different forces. There are two general classes of forces enabled (see `base_forces.hpp`): forces which are computed for all vertices of a given type (`VertexForce` class), and forces which are computed for all half-edges of a given type (`HalfEdgeForce` class).

Proper definitions of forces belong in `forces.hpp`, with `forces.cpp` providing definitions of `VertexModel` methods to add them. Finally, these are exposed to Python in `bind.cpp`.

Descriptions of some forces can be found in `docs/forces.pdf`.

### Pickling

`VertexModel` instances can be pickled via Python (see `system_pickle.hpp`). However this pickling does not save the state nor the definitions of the forces, which then have to be redefined after unpickling.

## Authors

- Yann-Edwin Keta, Instituut-Lorentz, Universiteit Leiden (keta@lorentz.leidenuniv.nl)
- Silke Henkes, Instituut-Lorentz, Universiteit Leiden (shenkes@lorentz.leidenuniv.nl)

