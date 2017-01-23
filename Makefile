.PHONY: test, clang 

test:
	g++ -g -std=c++14 -Wall -Wpedantic -Wextra -obin_test.out -L/usr/lib/x86_64-linux-gnu/ tests.cpp -lboost_unit_test_framework 
	./bin_test.out --report_level=detailed

clang:
	clang++ -std=c++14 -Wall -Wpedantic -Wextra -fsyntax-only tests.cpp 
