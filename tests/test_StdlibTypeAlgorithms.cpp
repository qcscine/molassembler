#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include "StdlibTypeAlgorithms.h"

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

using namespace StdlibTypeAlgorithms;

bool firstSmallerSecond(
  const unsigned& a,
  const unsigned& b
) {
  return a < b;
}

template<typename T>
bool vectorOfSetsEqual(
  const std::vector<
    std::set<T>
  >& a,
  const std::vector<
    std::set<T>
  >& b
) {
  return std::all_of(
    a.begin(),
    a.end(),
    [&b](const auto& setI) {
      return std::accumulate(
        b.begin(),
        b.end(),
        false,
        [&setI](const bool& carry, const auto& setJ) {
          return (
            carry 
            || std::equal(
              setI.begin(),
              setI.end(),
              setJ.begin(),
              setJ.end()
            )
          );
        }
      );
    }
  );
}

BOOST_AUTO_TEST_CASE( stdlibTypeAlgorithms ) {
  /* 1, 2 (2 just calls 1) */
  std::vector<
    std::set<unsigned>
  > testSetList {
    {5, 2, 3, 9, 11, 4},
    {2, 1, 0, 12},
    {13, 6}
  };

  BOOST_CHECK(
    vectorOfSetsEqual(
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
    vectorOfSetsEqual(
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
