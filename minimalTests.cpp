#define BOOST_TEST_MODULE CyclicPolygonsTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Minimal.h"
#include <iostream>
#include "temple/VectorView.h"
#include "temple/Random.h"

BOOST_AUTO_TEST_CASE(centralAngleRootFinding) {
  const double upperLimit = 5.6; // Fr-Fr single
  const double lowerLimit = 0.7; // H-H single

  const unsigned nTests = 100;

  for(unsigned nSides = 5; nSides < 10; nSides++) {
    for(unsigned testNumber = 0; testNumber < nTests; testNumber++) {
      std::vector<double> edgeLengths = temple::random.getN<double>(lowerLimit, upperLimit, nSides);
      while(!CyclicPolygons::exists(edgeLengths)) {
        edgeLengths = temple::random.getN<double>(lowerLimit, upperLimit, nSides);
      }

      double circumradius;
      auto assignCircumradius = [&]() -> void {
        circumradius = CyclicPolygons::detail::convexCircumradius(
          edgeLengths
        );
      };

      BOOST_CHECK_NO_THROW(assignCircumradius());
      BOOST_CHECK(
        std::fabs(
          CyclicPolygons::detail::centralAnglesDeviation(circumradius, edgeLengths)
        ) < 1e-6
      );
      BOOST_CHECK(
        std::fabs(
          temple::sum(
            CyclicPolygons::detail::generalizedInternalAngles(edgeLengths, circumradius)
          ) - 3 * M_PI
        ) < 1e-6
      );
    }
  }
}
