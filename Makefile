CFLAGS=--std=c++20

all: test

test: system.*pp mesh.*pp test.*pp
	g++ $(CFLAGS) -c system.cpp
	g++ $(CFLAGS) -c mesh.cpp
	g++ $(CFLAGS) -c test.cpp
	g++ $(CFLAGS) -c tools.cpp
	g++ $(CFLAGS) -o test system.o mesh.o test.o tools.o
