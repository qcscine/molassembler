.PHONY: all, clang
all:
	g++ -g -O2 -std=c++14 -Wall -Wpedantic -Wextra -obin_test -L/usr/lib/x86_64-linux-gnu/ tests.cpp -lboost_unit_test_framework
	./bin_test --report_level=detailed

clang:
	clang++ -std=c++14 -Wall -Wpedantic -Wextra -fsyntax-only tests.cpp 

