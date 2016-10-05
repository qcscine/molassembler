#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

#include "algorithm.h"
#include "SymmetryInformation.h"

/* TODO
 * - investigate failures!
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

BOOST_AUTO_TEST_CASE( assignment_instantiation ) {
  // can make instances of all symmetries
  Assignment<PermSymmetry::Tetrahedral> tetr(
    std::vector<char>(4, 'A')
  );
  Assignment<PermSymmetry::SquarePlanar> sqpl(
    std::vector<char>(4, 'A')
  );
  Assignment<PermSymmetry::SquarePyramidal> sqpy(
    std::vector<char>(5, 'A')
  );
  Assignment<PermSymmetry::TrigonalBiPyramidal> trigbipy(
    std::vector<char>(5, 'A')
  );
  Assignment<PermSymmetry::Octahedral> octa(
    std::vector<char>(6, 'A')
  );

  // and use them polymorphically
  std::vector<
    std::shared_ptr<AbstractAssignment>
  > container = {
    std::make_shared<
      Assignment<PermSymmetry::Tetrahedral>
    >(std::vector<char>(4, 'A')),
    std::make_shared<
      Assignment<PermSymmetry::SquarePlanar>
    >(std::vector<char>(4, 'A')),
    std::make_shared<
      Assignment<PermSymmetry::SquarePyramidal>
    >(std::vector<char>(5, 'A')),
    std::make_shared<
      Assignment<PermSymmetry::TrigonalBiPyramidal>
    >(std::vector<char>(5, 'A')),
    std::make_shared<
      Assignment<PermSymmetry::Octahedral>
    >(std::vector<char>(6, 'A')),
  };
}

BOOST_AUTO_TEST_CASE( assignment_basics ) {
  // Constructors
  Assignment<PermSymmetry::Octahedral> instance(
    {'A', 'A', 'A', 'A', 'A', 'A'}
  );
  Assignment<PermSymmetry::Octahedral> instanceWithBondedLigands(
    {'A', 'A', 'A', 'A', 'A', 'A'},
    {
      std::make_pair(0,1),
      std::make_pair(2,3),
      std::make_pair(4,5)
    }
  );
  // all constructors sort the ligands
  BOOST_CHECK(
    instance.ligandConnectionsAreOrdered() == true
  );
  BOOST_CHECK(
    instanceWithBondedLigands.ligandConnectionsAreOrdered() == true
  );
}

template<
  template<typename T>
  class Symmetry
>
void run_tests(
  const std::vector<
    std::tuple<
      std::vector<char>, // characters
      std::vector<
        std::pair<unsigned, unsigned>
      >, // pairs
      unsigned // expectedUnique
    >
  >& test_cases
) {
  for(const auto& tuple: test_cases) {
    std::vector<char> characters;
    std::vector<
      std::pair<unsigned, unsigned>
    > pairs;
    unsigned expectedUnique;

    std::tie(characters, pairs, expectedUnique) = tuple;

    // instantiate
    Assignment<Symmetry> assignment = (pairs.size() == 0)
      ? Assignment<Symmetry>(characters)
      : Assignment<Symmetry>(characters, pairs);

    auto unique = uniqueAssignments(assignment);
    BOOST_CHECK(unique.size() == expectedUnique );
    if(unique.size() != expectedUnique) {
      std::cout << "Mismatch: Expected " << expectedUnique
        << " assignments for: " << assignment << ", got " 
        << unique.size() << std::endl;
      std::cout << "Unique assignments: " << std::endl;
      for(const auto& uniqueAssignment: unique) {
        std::cout << uniqueAssignment << std::endl;
      }
    }
  }
}

/* Octahedral tests */
BOOST_AUTO_TEST_CASE( octahedral_monodentate ) {
  run_tests<PermSymmetry::Octahedral>(
    {
      // M_A
      std::make_tuple(
        std::vector<char>(6, 'A'),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        1
      ),
      // M_AB
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'A', 'B'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        1
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'B', 'B'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        2
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'B', 'B', 'B'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        2
      ),
      // M_ABC
      /* A4 B C
       * A3 B2 C
       * not A3 B C2 == A3 B2 C
       * A2 B2 C2
       * not A2 B3 C == A3 B2 C
       */
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'B', 'C'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        2
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'B', 'B', 'C'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        3
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'B', 'C', 'C'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        6
      ),
      // M_ABCD
      /* A3 B C D
       * A2 B2 C D
       */
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'B', 'C', 'D'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        5
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'B', 'C', 'D'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        8
      ),
      // M_ABCDE
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'C', 'D', 'E'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        15
      ),
      // M_ABCDEF
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'D', 'E', 'F'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >(),
        30
      )
    }
  );
}

BOOST_AUTO_TEST_CASE( octahedral_multidentate ) {
  run_tests<PermSymmetry::Octahedral>(
    {
      // M(A-A)_3
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'A', 'A'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3),
          std::make_pair(4, 5)
        }),
        2
        /* NOTE: besides the two cis-cis-cis enantiomers, there are
         * cis-cis-trans and trans-trans-trans isomers, these are removed by
         * default!
         */
      ),
      // M(A-B)_3
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'B', 'A', 'B'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3),
          std::make_pair(4, 5)
        }),
        4
      ),
      // M(A-B)_2 CD
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'B', 'C', 'D'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        11
      ),
      // M(A-A)(B-C)DE
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'C', 'D', 'E'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        10
      ),
      // M(A-B)(C-D)EF
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'D', 'E', 'F'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        20
      ),
      // M(A-B-A)CDE
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'C', 'D', 'E'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2)
        }),
        9
      ),
      // M(A-B-C)_2
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'A', 'B', 'C'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2),
          std::make_pair(3, 4),
          std::make_pair(4, 5)
        }),
        11 
      ),
      // M(A-B-B-A)CD
      std::make_tuple(
        std::vector<char>({'A', 'B', 'B', 'A', 'C', 'D'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2),
          std::make_pair(2, 3)
        }),
        7
      ),
      // M(A-B-C-B-A)D
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'B', 'A', 'D'}),
        std::vector<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2),
          std::make_pair(2, 3),
          std::make_pair(3, 4)
        }),
        7
      )
    }
  );
}
