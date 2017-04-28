.PHONY: all

all:
	g++ -Wall -Wpedantic -std=c++14 test.cpp
	./a.out
