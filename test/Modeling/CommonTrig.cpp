/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "boost/test/unit_test.hpp"

#include "Molassembler/Modeling/CommonTrig.h"
#include "Molassembler/Options.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"

template<typename T>
bool isApproxAbs(
  const T& a,
  const T& b
) {
  return std::fabs(a - b) < 1e-10;
}


BOOST_AUTO_TEST_CASE(CommonTrigExampleTests, *boost::unit_test::label("Molassembler")) {
  using namespace Scine;
  using namespace Molassembler;
  using namespace CommonTrig;

  BOOST_CHECK(
    isApproxAbs(
      lawOfCosines(1.0, 1.0, M_PI / 2),
      std::sqrt(2)
    )
  );

  BOOST_CHECK_MESSAGE(
    isApproxAbs(
      dihedralLength(1, 1, 1, M_PI / 2, M_PI / 2, 0),
      1.0
    ),
    "It's not 1, it's "
      << dihedralLength(1, 1, 1, M_PI / 2, M_PI / 2, 0)
  );

  BOOST_CHECK_MESSAGE(
    isApproxAbs(
      dihedralLength(1, 1, 1, M_PI / 2, M_PI / 2, M_PI),
      std::sqrt(5)
    ),
    "It's not " << std::sqrt(5) <<", it's "
      << dihedralLength(1, 1, 1, M_PI / 2, M_PI / 2, M_PI)
  );
}

BOOST_AUTO_TEST_CASE(DihedralZeroAlwaysSmallerDihedralPi, *boost::unit_test::label("Molassembler")) {
  using namespace Scine;
  using namespace Molassembler;
  using namespace CommonTrig;

  for(unsigned testNum = 0; testNum < 100; ++testNum) {
    const auto sideLengths = Temple::Random::getN<double>(1.4, 5.6, 3, randomnessEngine());
    const auto angles = Temple::Random::getN<double>(0.0, M_PI, 2, randomnessEngine());

    BOOST_CHECK_MESSAGE(
      dihedralLength(
        sideLengths.at(0),
        sideLengths.at(1),
        sideLengths.at(2),
        angles.at(0),
        angles.at(1),
        0
      ) < dihedralLength(
        sideLengths.at(0),
        sideLengths.at(1),
        sideLengths.at(2),
        angles.at(0),
        angles.at(1),
        M_PI
      ),
      "Disproved by side lengths {"
        << Temple::condense(sideLengths)
        << "} and angles {"
        << Temple::condense(angles)
        << "}: 0 -> "
        << dihedralLength(
          sideLengths.at(0),
          sideLengths.at(1),
          sideLengths.at(2),
          angles.at(0),
          angles.at(1),
          0
        ) << ", π -> "
        << dihedralLength(
          sideLengths.at(0),
          sideLengths.at(1),
          sideLengths.at(2),
          angles.at(0),
          angles.at(1),
          M_PI
        )
    );
  }
}
