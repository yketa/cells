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

mesh.o: mesh.hpp tools.hpp

system.o: system.hpp tools.hpp

tools.o: tools.hpp

# libraries

bind.so: CFLAGS+=-fPIC `$(PYTHON) -m pybind11 --includes`
bind.so: bind.o system.o mesh.o tools.o
	$(CC) -o $@ -shared $^ $(LDFLAGS)

