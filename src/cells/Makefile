SHELL:=/bin/bash

# COMPILATION PARAMETERS

# -- python executable --
PYTHON=python
# -- compiler --
# default C++ compiler (g++, clang++, etc.) should already be available in Makefile as $(CXX)
# CXX=g++

# -- compiler flags --
# -g flag should not affect performance (https://stackoverflow.com/questions/10988318/g-is-using-the-g-flag-for-production-builds-a-good-idea)
CFLAGS=-std=c++20 -O3 -Wall -g
# -- linker flags --
# linker flags may be OS-dependent (https://stackoverflow.com/questions/714100/os-detecting-makefile)
# macOS needs specific linker flag (https://pybind11.readthedocs.io/en/stable/compiling.html#building-manually)
LDFLAGS=$(shell [[ `uname -s` == Darwin ]] && echo -undefined dynamic_lookup)
# -- openMP flags --
MPFLAGS=-fopenmp
# -- python flags --
# find actual python executable by resolving symlinks
# (https://stackoverflow.com/questions/7665/how-to-resolve-symbolic-links-in-a-shell-script/42918#42918)
PYTHON:=$(shell perl -MCwd -le 'print Cwd::abs_path(shift)' "`which $(PYTHON)`")
PYTHONFLAGS=
# python libraries
# PYTHONFLAGS+=$(shell $(PYTHON)-config --includes)
PYTHONFLAGS+=$(shell $(PYTHON) -m pybind11 --includes)

.PHONY: clean

all: bind.so

clean:
	rm -rf *.o *.so

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

