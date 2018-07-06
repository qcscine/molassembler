#define BOOST_TEST_MODULE StdlibTypeAlgorithmTests
#include <boost/test/unit_test.hpp>

#include "detail/StdlibTypeAlgorithms.h"
#include "temple/Containers.h"

/* Algorithms to test
 *
 * t # Name
 *   1 mergeOverlappingSetsInplace
 *   2 mergeOverlappingSets
 *   3 makeIndividualSets
 *   4 copyMerge
 *   5 minMaxAdaptor
 *   6 makeFunction
 */

bool firstSmallerSecond(
  const unsigned& a,
  const unsigned& b
) {
  return a < b;
}

BOOST_AUTO_TEST_CASE( stdlibTypeAlgorithms ) {
  using namespace StdlibTypeAlgorithms;

  /* 1, 2 (2 just calls 1) */
  std::vector<
    std::set<unsigned>
  > testSetList {
    {5, 2, 3, 9, 11, 4},
    {2, 1, 0, 12},
    {13, 6}
  };

  BOOST_CHECK(
    StdlibTypeAlgorithms::vectorOfSetsEqual(
      mergeOverlappingSets(
        testSetList
      ),
      {
        {0, 1, 2, 3, 4, 5, 9, 11, 12},
        {6, 13}
      }
    )
  );

  /* 3 */
  std::set<
    std::pair<unsigned, unsigned>
  > pairsSet {
    {1, 2},
    {2, 3},
    {4, 5},
    {5, 7},
    {7, 6} // can't remember if order is important or not, try it out
  };

  BOOST_CHECK(
    StdlibTypeAlgorithms::vectorOfSetsEqual(
      makeIndividualSets(
        pairsSet
      ),
      {
        {1, 2, 3},
        {4, 5, 6, 7}
      }
    )
  );

  /* 4 */
  std::vector<unsigned>
    a {1, 4, 7},
    b {2, 9, 3},
    expectedMerge {1, 4, 7, 2, 9, 3};
  auto merged = copyMerge(a, b);
  BOOST_CHECK(
    std::equal(
      merged.begin(),
      merged.end(),
      expectedMerge.begin(),
      expectedMerge.end()
    )
  );


  /* 5, 6 */
  BOOST_CHECK(
    minMaxAdaptor(
      makeFunction(firstSmallerSecond),
      7u,
      5u
    )
  );
}

BOOST_AUTO_TEST_CASE( combinationPermutation ) {
  const std::vector<unsigned> testLimits {4, 1, 3, 6, 9, 2};
  std::vector<unsigned> combination (testLimits.size(), 0);

  bool alwaysSmallerOrEqual = true;
  while(StdlibTypeAlgorithms::nextCombinationPermutation(combination, testLimits)) {
    if(
      !temple::all_of(
        temple::zipMap(
          combination,
          testLimits,
          [](const unsigned& index, const unsigned& limit) -> bool {
            return index <= limit;
          }
        )
      )
    ) {
      alwaysSmallerOrEqual = false;
      std::cout << "Falsified for combination {"
        << temple::condenseIterable(combination)
        << "}" << std::endl;
      break;
    }
  }

  BOOST_CHECK(alwaysSmallerOrEqual);
}
