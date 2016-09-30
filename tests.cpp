#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

#include "algorithm.hpp"

/* TODO
 * - investigate failures!
 * - Issues hindering generalizability
 *   Custom rotation algorithm with the symmetry operations
 */

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
    const std::vector<
        std::vector<char>
    >& test_cases,
    const std::vector<unsigned>& expected
) {
    assert(test_cases.size() == expected.size());
    for(unsigned i = 0; i < test_cases.size(); i++) {
        auto unique = unique_assignments(
            Assignment(test_cases[i])
        );
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

void run_test(
    const Assignment& instance,
    const unsigned& expected
) {
    auto unique = unique_assignments(
        instance
    );
    BOOST_CHECK ( unique.size() == expected );
    if(unique.size() != expected) {
        std::cout << "Mismatch: Expected " << expected
            << " assignments for: " << std::endl;
        std::cout << instance << std::endl;
        std::cout << ", got " << unique.size() << std::endl;
        std::cout << "Unique assignments: " << std::endl;
        for(const auto& assignment: unique) {
            std::cout << assignment << std::endl;
        }
    }
}


BOOST_AUTO_TEST_CASE( assignment_column ) {
    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {true, false, true}
        ) == AssignmentColumn(
            'A',
            {true, false, true}
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {true, false, true}
        ) != AssignmentColumn(
            'A',
            {false, false, true}
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {true, false, true}
        ) != AssignmentColumn(
            'B',
            {true, false, true}
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {true, false, true}
        ) < AssignmentColumn(
            'B',
            {true, false, true}
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {true, false, false}
        ) < AssignmentColumn(
            'A',
            {true, false, true}
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {false, false, true}
        ) < AssignmentColumn(
            'A',
            {true, false, false}
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            std::vector<bool>()
        ) == AssignmentColumn(
            'A',
            std::vector<bool>()
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            std::vector<bool>()
        ) < AssignmentColumn(
            'B',
            std::vector<bool>()
        )
    );

    BOOST_CHECK(
        AssignmentColumn(
            'A',
            {true}
        ) < AssignmentColumn(
            'B',
            {false}
        )
    );
    BOOST_CHECK(
        (
            AssignmentColumn(
                'B',
                {false}
            ) < AssignmentColumn(
                'A',
                {true}
            )
        ) == false
    );
}

BOOST_AUTO_TEST_CASE( assignment_basics ) {
    // Constructors
    Assignment instance(
        {'A', 'A', 'A', 'A', 'A', 'A'}
    );
    Assignment instance_with_bonded_ligands(
        {'A', 'A', 'A', 'A', 'A', 'A'},
        {
            {true, true, false, false, false, false},
            {false, false, true, true, false, false},
            {false, false, false, false, true, true}
        }
    );
    Assignment instance_with_ordered_ligands(
        {'A', 'A', 'A', 'A', 'A', 'A'},
        {
            {false, false, false, false, true, true},
            {false, false, true, true, false, false},
            {true, true, false, false, false, false}
        }
    );

    // all constructs sort the ligands
    BOOST_CHECK(
        instance.ligand_connections_are_ordered() == true
    );
    BOOST_CHECK(
        instance_with_bonded_ligands.ligand_connections_are_ordered() == true
    );
    BOOST_CHECK(
        instance_with_ordered_ligands.ligand_connections_are_ordered() == true
    );

    BOOST_CHECK(
        instance_with_bonded_ligands == instance_with_ordered_ligands
    );

    BOOST_CHECK(
        instance_with_bonded_ligands.is_rotationally_superimposable(
            instance_with_ordered_ligands
        )
    );
}

// all cases

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

// and (some) multidentate cases

BOOST_AUTO_TEST_CASE( M_AA3 ) {
    run_test(
        Assignment(
            {'A', 'A', 'A', 'A', 'A', 'A'},
            {
                {true, true, false, false, false, false},
                {false, false, true, true, false, false},
                {false, false, false, false, true, true}
            }
        ),
        2
        /* besides the two cis-cis-cis enantiomers, there are
         * cis-cis-trans and trans-trans-trans isomers, these are
         * removed by default!
         */
    );
}

BOOST_AUTO_TEST_CASE( M_AB3 ) {
    run_test(
        Assignment(
            {'A', 'B', 'A', 'B', 'A', 'B'},
            {
                {true, true, false, false, false, false},
                {false, false, true, true, false, false},
                {false, false, false, false, true, true}
            }
        ),
        4
    );
}

BOOST_AUTO_TEST_CASE( M_AB2_C_D ) {
    run_test(
        Assignment(
            {'A', 'B', 'A', 'B', 'C', 'D'},
            {
                {true, true, false, false, false, false},
                {false, false, true, true, false, false}
            }
        ),
        11
    );
}

BOOST_AUTO_TEST_CASE( M_AA_BC_D_E ) {
    run_test(
        Assignment(
            {'A', 'A', 'B', 'C', 'D', 'E'},
            {
                {true, true, false, false, false, false},
                {false, false, true, true, false, false}
            }
        ),
        10
    );
}

BOOST_AUTO_TEST_CASE( M_AB_CD_E_F ) {
    run_test(
        Assignment(
            {'A', 'B', 'C', 'D', 'E', 'F'},
            {
                {true, true, false, false, false, false},
                {false, false, true, true, false, false}
            }
        ),
        20
    );
}

BOOST_AUTO_TEST_CASE( M_ABA_C_D_E ) {
    run_test(
        Assignment(
            {'A', 'B', 'A', 'C', 'D', 'E'},
            {
                {true, true, false, false, false, false},
                {false, true, true, false, false, false}
            }
        ),
        9
    );
}

BOOST_AUTO_TEST_CASE( M_ABC2 ) {
    run_test(
        Assignment(
            {'A', 'B', 'C', 'A', 'B', 'C'},
            {
                {true, true, false, false, false, false},
                {false, true, true, false, false, false},
                {false, false, false, true, true, false},
                {false, false, false, false, true, true}
            }
        ),
        11 // fails, get 9
    );
}

BOOST_AUTO_TEST_CASE( M_ABBA_C_D ) {
    run_test(
        Assignment(
            {'A', 'B', 'B', 'A', 'C', 'D'},
            {
                {true, true, false, false, false, false},
                {false, true, true, false, false, false},
                {false, false, true, true, false, false}
            }
        ),
        7
    );
}

BOOST_AUTO_TEST_CASE( M_ABCBA_D ) {
    run_test(
        Assignment(
            {'A', 'B', 'C', 'B', 'A', 'D'},
            {
                {true, true, false, false, false, false},
                {false, true, true, false, false, false},
                {false, false, true, true, false, false},
                {false, false, false, true, true, false}
            }
        ),
        7 // fails, get 5
    );
}
