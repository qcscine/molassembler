.PHONY: all

all:
	g++ -Wall -Wpedantic -std=c++14 -Istatic_math/include test.cpp
	./a.out
