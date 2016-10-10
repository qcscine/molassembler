#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "numeric.h"

template<typename T>
std::ostream& operator << (
    std::ostream& os,
    const std::vector<T>& vector
) {
    for(const auto& element: vector) {
        os << element << std::endl;
    }
    return os;
}


BOOST_AUTO_TEST_CASE ( average ) {
    std::vector<
        std::vector<unsigned>
    > test_cases = {
        {0},
        {0, 2},
        {1},
        {1,2,3},
        {1,2,3,4,5}
    };

    std::vector<unsigned> expected_result = {
        0,
        1,
        1,
        2,
        3
    };

    for(unsigned i = 0; i < test_cases.size(); i++) {
        BOOST_CHECK( numeric::average(test_cases[i]) == expected_result[i] );
    }
}

BOOST_AUTO_TEST_CASE ( stdev ) {
    std::vector<
        std::vector<unsigned>
    > test_cases = {
        {1,5,7,3,2,9},
        {0},
        {10,20,30,50,70,90},
        {0,0,1},
        {4,4,4}
    };

    std::vector<unsigned> expected_result = {
        3, 
        0,
        28,
        0, 
        0
    };

    for(unsigned i = 0; i < test_cases.size(); i++) {
        BOOST_CHECK( 
            numeric::population_stddev(test_cases[i]) == expected_result[i] 
        );
    }
}
