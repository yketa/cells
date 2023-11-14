SHELL:=/bin/bash

CC=g++
CFLAGS=-std=gnu++20 -O3 -Wall
LDFLAGS=
MPIFLAGS=

all: test

clean:
	rm *.o

# object files

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

system.o: system.hpp tools.hpp

mesh.o: mesh.hpp tools.hpp

test.o: system.hpp

tools.o: tools.hpp

bind.o: mesh.hpp system.hpp tools.hpp

# executables

test: system.o mesh.o test.o tools.o
	$(CC) -o $@ $^ $(LDFLAGS)

# libraries

bind.so: CFLAGS+=-fPIC `python -m pybind11 --includes`
bind.so: bind.o system.o mesh.o tools.o
	$(CC) -o $@ -shared $^ $(LDFLAGS)

