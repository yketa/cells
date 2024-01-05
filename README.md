# Cells
*Yann-Edwin Keta (Instituut-Lorentz, Universiteit Leiden — 2023)*

Python/C++ integration of vertex models. Developed under the supervision of Silke Henkes.

## Compilation and running

Compilation of the shared library with `make bind.so` requires [`pybind11`](https://pybind11.readthedocs.io) (`python -m pip install pybind11`) and `g++` or an equivalent C++ compiler.

Python scripts are written for `python3.*` and import the `cells` package which necessitates the directory containing this repository to be added to the `$PYTHONPATH`, e.g. by executing
```
echo "export PYTHONPATH=\$PYTHONPATH:${PWD}/.." >> ~/.bashrc
```
from this directory.

Test with `python test.py` requires `numpy` and `matplotlib`. These can be installed with `python -m pip install -r requirements.txt`.

## C++ scripts structure

### Vertex model

This vertex model implementation is separated in two parts. First `mesh.*pp` contains the implementation of the geometrical features of the model, which relies on vertices which are linked together by half-edges enabling to move accross the mesh (see `docs/avm.pdf`). Second `system.*pp` provides a general `VertexModel` class to integrate the dynamics.

`VertexModel` methods to initialise configurations are defined in `initialisation.cpp`.

### Forces

Forces are defined separately and are attached to the `VertexModel` class following a class factory method (see `class_factory.hpp`) which enables to add and remove different forces. There are two general classes of forces enabled (see `base_forces.hpp`): forces which are computed for all vertices (`VertexForce` class), and forces which are computed for all half-edges (`HalfEdgeForce` class).

Proper definitions of forces belong in `forces.hpp`, with `forces.cpp` providing definitions of `VertexModel` methods to add them.

### Python binding

`VertexModel` is exposed to Python through pybind11 in `bind.cpp`. Importantly, this file also contains additional methods to provide control of the forces in Python.

### Pickling

Work in progress...

