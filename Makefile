SHELL:=/bin/bash

# COMPILATION PARAMETERS

# default C++ compiler should already be available in Makefile through variable $(CXX)
PYTHON=python
CFLAGS=-std=c++20 -O3 -Wall
LDFLAGS=
MPIFLAGS=-fopenmp

# find actual python executable by resolving symlinks
# (https://stackoverflow.com/a/42918/7385044)
PYTHON:=$(shell perl -MCwd -le 'print Cwd::abs_path(shift)' "`which $(PYTHON)`")

all: bind.so

clean:
	rm -f *.o *.so

# OBJECT FILES

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

bind.o: mesh.hpp system.hpp tools.hpp pickle.hpp plot.hpp

initialisation.o: system.hpp tools.hpp

mesh.o: mesh.hpp tools.hpp

system.o: system.hpp tools.hpp

tools.o: tools.hpp

# LIBRARIES

bind.so: CFLAGS+=-fPIC `$(PYTHON) -m pybind11 --includes` `$(PYTHON)-config --includes`
bind.so: bind.o initialisation.o system.o mesh.o tools.o
	$(CXX) -o $@ -shared $^ $(LDFLAGS)

