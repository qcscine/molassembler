#define BOOST_TEST_MODULE CommonTrigTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "CommonTrig.h"
#include "template_magic/Random.h"
#include "template_magic/TemplateMagic.h"

template<typename T>
bool isApproxAbs(
  const T& a,
  const T& b
) {
  return std::fabs(a - b) < 1e-10;
}


BOOST_AUTO_TEST_CASE(randomExamples) {
  using namespace MoleculeManip::CommonTrig;

  BOOST_CHECK(
    isApproxAbs(
      lawOfCosines(
        1.0,
        1.0,
        M_PI / 2
      ),
      std::sqrt(2)
    )
  );

  BOOST_CHECK_MESSAGE(
    isApproxAbs(
      dihedralLength(
        1,
        1,
        1,
        M_PI / 2,
        M_PI / 2,
        0
      ),
      1.0
    ),
    "It's not 1, it's "
      << dihedralLength(
        1,
        1,
        1,
        M_PI / 2,
        M_PI / 2,
        0
      )
  );

  BOOST_CHECK_MESSAGE(
    isApproxAbs(
      dihedralLength(
        1,
        1,
        1,
        M_PI / 2,
        M_PI / 2,
        M_PI
      ),
      std::sqrt(5)
    ),
    "It's not " << std::sqrt(5) <<", it's " 
      << dihedralLength(
        1,
        1,
        1,
        M_PI / 2,
        M_PI / 2,
        M_PI
      )
  );
}

BOOST_AUTO_TEST_CASE(dihedralZeroAlwaysSmallerDihedralPi) {
  using namespace MoleculeManip::CommonTrig;

  for(unsigned testNum = 0; testNum < 100; ++testNum) {
    const auto sideLengths = TemplateMagic::random.getN<double>(1.4, 5.6, 3);
    const auto angles = TemplateMagic::random.getN<double>(0.0, M_PI, 2);

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
        << TemplateMagic::condenseIterable(sideLengths)
        << "} and angles {" 
        << TemplateMagic::condenseIterable(angles)
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
