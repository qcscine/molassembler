.PHONY: all 
all:
	g++ -O2 -std=c++14 -otest -L/usr/lib/x86_64-linux-gnu/ tests.cpp -lboost_unit_test_framework
	./test --report_level=detailed
