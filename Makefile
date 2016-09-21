.PHONY:	test 

test:
	g++ -std=c++14 -otest -L/usr/lib/x86_64-linux-gnu/ test_AdjacencyList.cpp -lboost_unit_test_framework
	./test
clangparse:
	clang++ -I/usr/include/boost -I/home/jan-grimo/PhD/code/delib2-0/include -fsyntax-only -Wall -Wextra -Wpedantic -std=c++14 test_AdjacencyList.cpp

