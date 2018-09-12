#include <boost/test/unit_test.hpp>

#include "temple/constexpr/Numeric.h"

#include <iostream>

BOOST_AUTO_TEST_CASE(numericAverageStdDev) {
  const std::vector<double> values {29, 30, 31, 32, 33};

  BOOST_CHECK(temple::average(values) == 31);
  BOOST_CHECK(
    std::fabs(temple::stddev(values) - std::sqrt(2))
    < 1e-10
  );
}
