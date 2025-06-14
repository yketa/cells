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
# https://stackoverflow.com/a/59074274
PYTHONFLAGS+=$(shell $(PYTHON)-config --includes 2> /dev/null)
PYTHONFLAGS+=$(shell $(PYTHON) -m pybind11 --includes 2> /dev/null)
# -- external libraries flags --
EXTERNFLAGS=
# use pkg-config to find libraries
EXTERNFLAGS+=$(shell pkg-config --cflags eigen3 2> /dev/null)
# as a backup search in extern
EXTERNFLAGS+=-Iextern

.PHONY: clean

all: bind.so

clean:
	rm -rf *.o *.so

# OBJECT FILES
# dependencies to be checked with `$(CXX) -MM $<`

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

assert.o: assert.cpp assert.hpp

bind.o: bind.cpp analysis.hpp system.hpp class_factory.hpp forces.hpp \
 base_forces.hpp base.hpp mesh.hpp random.hpp tools.hpp integrators.hpp \
 assert.hpp base_integrators.hpp forces_pickle.hpp base_pickle.hpp \
 integrators_pickle.hpp plot.hpp system_pickle.hpp

forces.o: forces.cpp forces.hpp base_forces.hpp base.hpp mesh.hpp \
 random.hpp tools.hpp system.hpp class_factory.hpp integrators.hpp \
 assert.hpp base_integrators.hpp

initialisation.o: initialisation.cpp assert.hpp system.hpp \
 class_factory.hpp forces.hpp base_forces.hpp base.hpp mesh.hpp \
 random.hpp tools.hpp integrators.hpp base_integrators.hpp

integrators.o: integrators.cpp integrators.hpp assert.hpp \
 base_integrators.hpp base.hpp mesh.hpp system.hpp class_factory.hpp \
 forces.hpp base_forces.hpp random.hpp tools.hpp

mesh.o: mesh.cpp assert.hpp mesh.hpp tools.hpp

system.o: system.cpp assert.hpp system.hpp class_factory.hpp forces.hpp \
 base_forces.hpp base.hpp mesh.hpp random.hpp tools.hpp integrators.hpp \
 base_integrators.hpp

tools.o: tools.cpp tools.hpp

bind.o forces.o initialisation.o integrators.o system.o: CFLAGS+=$(PYTHONFLAGS) $(EXTERNFLAGS)

# LIBRARIES

bind.so: CFLAGS+=-fPIC
bind.so: assert.o bind.o forces.o initialisation.o mesh.o integrators.o system.o tools.o
	$(CXX) -o $@ -shared $^ $(LDFLAGS)

