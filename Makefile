# Simple Makefile for testing projFuncs.cpp

all: projFuncs

clean:
	rm -f projFuncs
	rm -f *.out

projFuncs:
	g++ projFuncs.cpp -I./ -std=c++17 -o projFuncs -O3