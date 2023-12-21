SHELL:=/bin/bash

CC=g++
PYTHON=python
CFLAGS=-std=gnu++20 -O3 -Wall
LDFLAGS=
MPIFLAGS=

all: bind.so

clean:
	rm -f *.o *.so

# object files

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

bind.o: mesh.hpp system.hpp tools.hpp pickle.hpp plot.hpp

initialisation.o: system.hpp tools.hpp

mesh.o: mesh.hpp tools.hpp

system.o: system.hpp tools.hpp

tools.o: tools.hpp

# libraries

bind.so: CFLAGS+=-fPIC `python -m pybind11 --includes`
bind.so: bind.o initialisation.o system.o mesh.o tools.o
	$(CC) -o $@ -shared $^ $(LDFLAGS)

