.PHONY: all, test, clang, clean
all:
	g++ -O2 -std=c++14 -Wall -Wpedantic -Wextra -obuild/SymmetryInformation.o -c SymmetryInformation.cpp
	ar rvs build/libAssignment.a build/SymmetryInformation.o
	rm build/SymmetryInformation.o

test:
	g++ -g -O2 -std=c++14 -Wall -Wpedantic -Wextra -obuild/bin_test -L/usr/lib/x86_64-linux-gnu/ -Lbuild/ tests.cpp -lboost_unit_test_framework -lAssignment
	./build/bin_test --report_level=detailed

clang:
	clang++ -std=c++14 -Wall -Wpedantic -Wextra -fsyntax-only tests.cpp 
