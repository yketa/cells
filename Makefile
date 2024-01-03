SHELL:=/bin/bash

# COMPILATION PARAMETERS

# default C++ compiler should already be available in Makefile through variable $(CXX)
PYTHON=python   # python executable
CFLAGS=-std=c++20 -O3 -Wall
LDFLAGS=
MPIFLAGS=-fopenmp

# find actual python executable by resolving symlinks
# (https://stackoverflow.com/a/42918/7385044)
PYTHON:=$(shell perl -MCwd -le 'print Cwd::abs_path(shift)' "`which $(PYTHON)`")
# python libraries
PYTHONFLAGS=$(shell $(PYTHON) -m pybind11 --includes) $(shell $(PYTHON)-config --includes)

all: bind.so

clean:
	rm -f *.o *.so

# OBJECT FILES

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

bind.o: mesh.hpp system.hpp tools.hpp pickle.hpp plot.hpp

forces.o: forces.hpp system.hpp

initialisation.o: system.hpp tools.hpp

mesh.o: mesh.hpp tools.hpp

system.o: system.hpp base_forces.hpp class_factory.hpp forces.hpp random.hpp tools.hpp

tools.o: tools.hpp

# LIBRARIES

bind.o: CFLAGS+=$(PYTHONFLAGS)
bind.so: CFLAGS+=-fPIC
bind.so: bind.o forces.o initialisation.o mesh.o system.o tools.o
	$(CXX) -o $@ -shared $^ $(LDFLAGS)

