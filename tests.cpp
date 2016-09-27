#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

#include "algorithm.hpp"

template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& vector) {
    for(unsigned i = 0; i < vector.size(); i++) {
        os << vector[i];
        if(i != vector.size() - 1) {
            os << ", ";
        }
    }

    return os;
}

void run_tests(
    const std::vector<Assignment>& test_cases,
    const std::vector<unsigned>& expected
) {
    for(unsigned i = 0; i < test_cases.size(); i++) {
        auto unique = unique_assignments(test_cases[i]);
        BOOST_CHECK ( unique.size() == expected[i] );
        if(unique.size() != expected[i]) {
            std::cout << "Mismatch: Expected " << expected[i]
                << " assignments for: " << test_cases[i] << ", got " 
                << unique.size() << std::endl;
            std::cout << "Unique assignments: " << std::endl;
            for(const auto& assignment: unique) {
                std::cout << assignment << std::endl;
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( M_A ) {
    run_tests(
        {
            {'A', 'A', 'A', 'A', 'A', 'A'}
        },{
            1
        }
    );
}

BOOST_AUTO_TEST_CASE( M_AB ) {
    run_tests(
        {
            {'A', 'A', 'A', 'A', 'A', 'B'},
            {'A', 'A', 'A', 'A', 'B', 'B'},
            {'A', 'A', 'A', 'B', 'B', 'B'}
        }, {
            1,
            2, // cis, trans
            2 // fac, mer
        }
    );
}

BOOST_AUTO_TEST_CASE( M_ABC ) {
    /* A4 B C
     * A3 B2 C
     * not A3 B C2 == A3 B2 C
     * A2 B2 C2
     * not A2 B3 C == A3 B2 C
     */
    run_tests(
        {
            {'A', 'A', 'A', 'A', 'B', 'C'}, // A4 B C
            {'A', 'A', 'A', 'B', 'B', 'C'}, // A3 B2 C
            {'A', 'A', 'B', 'B', 'C', 'C'} // A2 B2 C2
        }, {
            2, // cis, trans
            3, // fac, mer-cis, mer-trans
            6 // ccc (2), cct, ctc, tcc, ttt
        }
    );
}

BOOST_AUTO_TEST_CASE( M_ABCD ) {
    /* A3 B C D
     * A2 B2 C D
     */
    run_tests(
        {
            {'A', 'A', 'A', 'B', 'C', 'D'}, // A3 B C D
            {'A', 'A', 'B', 'B', 'C', 'D'} // A2 B2 C D
        }, {
            5,
            8
        }
    );
}

BOOST_AUTO_TEST_CASE( M_ABCDE ) {
    /* A2 B C D E
     */
    run_tests(
        {
            {'A', 'A', 'B', 'C', 'D', 'E'} // A2 B C D
        }, {
            15
        }
    );
}

BOOST_AUTO_TEST_CASE( M_ABCDEF ) {
    run_tests(
        {
            {'A', 'B', 'C', 'D', 'E', 'F'} 
        }, {
            30
        }
    );
}
