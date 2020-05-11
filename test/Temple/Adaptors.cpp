/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Adaptors/All.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include <iostream>

using namespace Scine::Molassembler;

template<class Container>
std::size_t iteratorDistance(const Container& container) {
  return std::distance(std::begin(container), std::end(container));
}

BOOST_AUTO_TEST_CASE(pairAdaptorTests) {
  const std::vector<unsigned> i {5, 3, 9, 11};
  const std::vector<unsigned> j {3, 4};

  auto adjacents = Temple::Adaptors::sequentialPairs(i);

  BOOST_CHECK(adjacents.size() == 3);
  BOOST_CHECK(iteratorDistance(adjacents) == adjacents.size());

  BOOST_CHECK(
    Temple::sum(
      Temple::Adaptors::transform(adjacents, std::plus<>())
    ) == 8 + 12 + 20
  );

  auto singlePairs = Temple::Adaptors::allPairs(i);

  BOOST_CHECK(singlePairs.size() == 6);
  BOOST_CHECK(
    iteratorDistance(singlePairs) == singlePairs.size()
  );

  BOOST_CHECK(
    Temple::sum(
      Temple::Adaptors::transform(singlePairs, std::plus<>())
    ) == 8 + 14 + 16 + 12 + 14 + 20
  );

  auto twoPairs = Temple::Adaptors::allPairs(i, j);

  BOOST_CHECK(twoPairs.size() == 8);
  BOOST_CHECK(
    iteratorDistance(twoPairs) == twoPairs.size()
  );

  BOOST_CHECK(
    Temple::sum(
      Temple::Adaptors::transform(twoPairs, std::plus<>())
    ) == 8 + 9 + 6 + 7 + 12 + 13 + 14 + 15
  );
}

BOOST_AUTO_TEST_CASE(iotaAdaptorTests) {
  auto a = Temple::Adaptors::range(5U);

  BOOST_CHECK(a.size() == 5);
  BOOST_CHECK(iteratorDistance(a) == a.size());
  BOOST_CHECK(Temple::sum(a) == 10U);

  auto b = Temple::Adaptors::range(4U, 7U);

  BOOST_CHECK(b.size() == 3);
  BOOST_CHECK(iteratorDistance(b) == b.size());
  BOOST_CHECK(Temple::sum(b) == 15U);
}

BOOST_AUTO_TEST_CASE(zipAdaptorTests) {
  const std::vector<unsigned> i {5, 3, 9, 11}, j {3, 4};

  auto zipRange = Temple::Adaptors::zip(i, j);

  BOOST_CHECK(zipRange.size() == 2);
  BOOST_CHECK(iteratorDistance(zipRange) == zipRange.size());
  BOOST_CHECK(
    Temple::sum(
      Temple::Adaptors::transform(zipRange, std::plus<>())
    ) == 15U
  );
}

BOOST_AUTO_TEST_CASE(transformAdaptorTests) {
  const std::vector<unsigned> i {5, 3, 9, 11};

  auto transformRange = Temple::Adaptors::transform(
    i,
    [](unsigned x) -> int {return static_cast<int>(x) - 10;}
  );

  BOOST_CHECK(transformRange.size() == 4);
  BOOST_CHECK(iteratorDistance(transformRange) == transformRange.size());
  BOOST_CHECK(
    Temple::sum(transformRange) == static_cast<int>(Temple::sum(i)) - 4 * 10
  );
}

BOOST_AUTO_TEST_CASE( enumerateTests) {
  std::vector<unsigned> testVec {5, 2, 3, 4};

  bool pass = true;
  for(const auto& enumPair : Temple::Adaptors::enumerate(testVec)) {
    if(testVec.at(enumPair.index) != enumPair.value) {
      pass = false;
      break;
    }
  }

  BOOST_CHECK(pass);

  auto weirdSum = Temple::sum(
    Temple::map(
      Temple::Adaptors::enumerate(testVec),
      [](const auto& enumPair) -> unsigned {
        return enumPair.index + enumPair.value;
      }
    )
  );

  BOOST_CHECK(weirdSum == 5 + 3 + 5 + 7);
}


BOOST_AUTO_TEST_CASE(compoundAdaptorOwnership) {
  auto pairsOfRange = Temple::Adaptors::allPairs(
    Temple::Adaptors::range(4U)
  );

  auto selfOwningRange = Temple::Adaptors::range(4U);
  auto referenceOwningPairs = Temple::Adaptors::allPairs(selfOwningRange);

  auto checkPairs = [](const auto& rangeObject) -> void {
    BOOST_CHECK(rangeObject.size() == 6);
    BOOST_CHECK(iteratorDistance(rangeObject) == rangeObject.size());
    BOOST_CHECK(
      Temple::invoke(std::plus<>(), *std::begin(rangeObject)) == 1U
    );

    BOOST_CHECK(
      Temple::sum(
        Temple::Adaptors::transform(
          Temple::Adaptors::allPairs(
            Temple::Adaptors::range(4U)
          ),
          [](const unsigned i, const unsigned j) -> unsigned {
            return i * j;
          }
        )
      ) == 11U
    );
  };

  checkPairs(pairsOfRange);
  checkPairs(referenceOwningPairs);

  std::vector<unsigned> i {1, 4, 9}, j {5, 2};

  auto pairFromTwoReferences = Temple::Adaptors::allPairs(i, j);
  auto pairFromTwoRValues = Temple::Adaptors::allPairs(
    std::vector<unsigned> {1, 4, 9},
    std::vector<unsigned> {5, 2}
  );
  auto pairFromMixed = Temple::Adaptors::allPairs(
    std::vector<unsigned> {1, 4, 9},
    j
  );

  auto checkTwoPairs = [](const auto& rangeObject) -> void {
    BOOST_CHECK(rangeObject.size() == 6);
    BOOST_CHECK(iteratorDistance(rangeObject) == rangeObject.size());
    BOOST_CHECK(
      Temple::invoke(std::plus<>(), *std::begin(rangeObject)) == 6
    );

    BOOST_CHECK(
      Temple::sum(
        Temple::Adaptors::transform(
          rangeObject,
          std::plus<>()
        )
      ) == 6 + 3 + 9 + 6 + 14 + 11
    );
  };

  checkTwoPairs(pairFromTwoReferences);
  checkTwoPairs(pairFromTwoRValues);
  checkTwoPairs(pairFromMixed);

}

BOOST_AUTO_TEST_CASE(adaptorShortRanges) {
  auto checkRangeLength = [](
    const auto& rangeObject,
    const unsigned expectedSize,
    const std::string& description
  ) {
    BOOST_CHECK_MESSAGE(
      rangeObject.size() == expectedSize,
      description << " size is " << rangeObject.size() << ", expected "
        << expectedSize
    );
    BOOST_CHECK_MESSAGE(
      iteratorDistance(rangeObject) == rangeObject.size(),
      description << " iterator distance is " << iteratorDistance(rangeObject)
        << ", expected equal to size (" << rangeObject.size() << ")"
    );
  };

  checkRangeLength(
    Temple::Adaptors::allPairs(std::vector<unsigned> {4}),
    0,
    "single-element all-pairs"
  );

  checkRangeLength(
    Temple::Adaptors::allPairs(std::vector<unsigned> {}),
    0,
    "no-element all-pairs"
  );

  checkRangeLength(
    Temple::Adaptors::allPairs(
      std::vector<unsigned> {4},
      std::vector<unsigned> {6}
    ),
    1,
    "one-one all-pairs"
  );

  checkRangeLength(
    Temple::Adaptors::allPairs(
      std::vector<unsigned> {},
      std::vector<unsigned> {6}
    ),
    0,
    "none-one all-pairs"
  );

  checkRangeLength(
    Temple::Adaptors::allPairs(
      std::vector<unsigned> {},
      std::vector<unsigned> {}
    ),
    0,
    "none-none all-pairs"
  );

  checkRangeLength(
    Temple::Adaptors::sequentialPairs(std::vector<unsigned> {4}),
    0,
    "one-element sequential pairs"
  );

  checkRangeLength(
    Temple::Adaptors::sequentialPairs(std::vector<unsigned> {}),
    0,
    "no-element sequential pairs"
  );
}

template<typename Range>
void checkRangeLengthTempl(
  const Range& rangeObject,
  const unsigned expectedSize,
  const std::string& description
) {
  BOOST_CHECK_MESSAGE(
    rangeObject.size() == expectedSize,
    description << " size is " << rangeObject.size() << ", expected "
      << expectedSize
  );
  BOOST_CHECK_MESSAGE(
    iteratorDistance(rangeObject) == rangeObject.size(),
    description << " iterator distance is " << iteratorDistance(rangeObject)
      << ", expected equal to size (" << rangeObject.size() << ")"
  );
}

BOOST_AUTO_TEST_CASE(frameAdaptorTest) {
  checkRangeLengthTempl(
    Temple::Adaptors::cyclicFrame<1>(std::vector<unsigned> {}),
    0,
    "no-element cyclic frame of size 1"
  );

  checkRangeLengthTempl(
    Temple::Adaptors::cyclicFrame<1>(std::vector<unsigned> {1}),
    1,
    "single-element cyclic frame of size 1"
  );

  checkRangeLengthTempl(
    Temple::Adaptors::cyclicFrame<1>(std::vector<unsigned> {1, 2}),
    2,
    "two-element cyclic frame of size 1"
  );

  checkRangeLengthTempl(
    Temple::Adaptors::cyclicFrame<2>(std::vector<unsigned> {1, 2}),
    2,
    "two-element cyclic frame of size 2"
  );

  checkRangeLengthTempl(
    Temple::Adaptors::cyclicFrame<2>(std::vector<unsigned> {1, 2, 3}),
    3,
    "three-element cyclic frame of size 2"
  );

  checkRangeLengthTempl(
    Temple::Adaptors::cyclicFrame<4>(std::vector<unsigned> {1, 2, 3}),
    0,
    "three-element cyclic frame of size 4"
  );

  BOOST_CHECK(
    Temple::sum(
      Temple::map(
        Temple::Adaptors::cyclicFrame<2>(std::vector<unsigned> {1, 2, 3}),
        [](unsigned i, unsigned j) -> unsigned {
          return i * j;
        }
      )
    ) == 2U + 6U + 3U
  );
}

BOOST_AUTO_TEST_CASE(filterAdaptorTests) {
  const auto filterDistance = iteratorDistance(
    Temple::Adaptors::filter(
      std::vector<unsigned> {1, 2, 3},
      [](const unsigned x) -> bool {return x % 2 == 0;}
    )
  );

  BOOST_CHECK_MESSAGE(
    filterDistance == 1,
    "Filter of 1,2,3 applying is_even isn't length one, but " << filterDistance
  );

  // BOOST_CHECK_MESSAGE(
  //   Temple::sum(
  //     Temple::map(
  //       Temple::Adaptors::allPairs(
  //         Temple::Adaptors::filter(std::vector<unsigned> {1, 2, 3, 4, 5}, [](const unsigned x) -> bool {return x < 4;})
  //       ),
  //       std::plus<>{}
  //     )
  //   ) == 3U + 4U + 5U,
  //   "AllPairs and filter do not work together!"
  // );
}
