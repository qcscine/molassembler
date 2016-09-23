.PHONY:	test perf

test:
	g++ -std=c++14 -otest -L/usr/lib/x86_64-linux-gnu/ test_ConnectivityManager.cpp -lboost_unit_test_framework
	./test --report_level=detailed

perf:
	g++ -std=c++14 -operf -L/usr/lib/x86_64-linux-gnu/ perf_ConnectivityManager.cpp -lboost_unit_test_framework
	./perf

clangparse:
	clang++ -I/usr/include/boost -I/home/jan-grimo/PhD/code/delib2-0/include -fsyntax-only -Wall -Wextra -Wpedantic -std=c++14 ConnectivityManager.hpp

clangparsetest:
	clang++ -I/usr/include/boost -I/home/jan-grimo/PhD/code/delib2-0/include -fsyntax-only -Wall -Wextra -Wpedantic -std=c++14 test_ConnectivityManager.cpp

clangparseperf:
	clang++ -I/usr/include/boost -I/home/jan-grimo/PhD/code/delib2-0/include -fsyntax-only -Wall -Wextra -Wpedantic -std=c++14 perf_ConnectivityManager.cpp

