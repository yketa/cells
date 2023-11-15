# Cells
### Yann-Edwin Keta (Instituut-Lorentz, Universiteit Leiden — 2023)

Compilation of the shared library with `make bind.so` requires [`pybind11`](https://pybind11.readthedocs.io) (`python -m pip install pybind11`) and `g++` or an equivalent compiler.

Python scripts are written for `python3.*` and import the `cells` package which necessitates the directory containing this repository to be added to the `$PYTHONPATH`, e.g. by executing
```
echo "export PYTHONPATH=\$PYTHONPATH:${PWD}/.." >> ~/.bashrc
```
from this directory.

Test with `python test.py` requires `numpy` and `matplotlib`.
