SHELL:=/bin/bash

# COMPILATION PARAMETERS

# -- python executable --
PYTHON=python
# -- compiler --
# default C++ compiler (g++, clang++, etc.) should already be available in Makefile through variable $(CXX)
# CXX=g++
# -- compiler flags --
# -g flag should not affect performance (https://stackoverflow.com/questions/10988318)
CFLAGS=-std=c++20 -O3 -Wall -g
# -- linker flags --
LDFLAGS=
# -- openMP flags --
MPFLAGS=-fopenmp

# find actual python executable by resolving symlinks
# (https://stackoverflow.com/a/42918/7385044)
PYTHON:=$(shell perl -MCwd -le 'print Cwd::abs_path(shift)' "`which $(PYTHON)`")
# python libraries
PYTHONFLAGS=$(shell $(PYTHON) -m pybind11 --includes) $(shell $(PYTHON)-config --includes)

.PHONY: clean

all: bind.so

clean:
	rm -rf *.o *.so *.sif

# OBJECT FILES
# dependencies to be checked with `$(CXX) -MM $<`

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

bind.o: forces.hpp base_forces.hpp mesh.hpp random.hpp tools.hpp forces_pickle.hpp base_pickle.hpp class_factory.hpp system.hpp plot.hpp system_pickle.hpp

forces.o: forces.hpp base_forces.hpp mesh.hpp random.hpp tools.hpp system.hpp class_factory.hpp

initialisation.o: system.hpp class_factory.hpp forces.hpp base_forces.hpp mesh.hpp random.hpp tools.hpp

mesh.o: mesh.hpp tools.hpp

system.o: system.hpp class_factory.hpp forces.hpp base_forces.hpp mesh.hpp random.hpp tools.hpp

tools.o: tools.hpp

bind.o forces.o initialisation.o system.o: CFLAGS+=$(PYTHONFLAGS)

# LIBRARIES

bind.so: CFLAGS+=-fPIC
bind.so: bind.o forces.o initialisation.o mesh.o system.o tools.o
	$(CXX) -o $@ -shared $^ $(LDFLAGS)

# CONTAINER

container.sif: container.def requirements.txt Makefile $(wildcard *.*pp) $(wildcard *.py)
	sudo singularity build $@ $<

