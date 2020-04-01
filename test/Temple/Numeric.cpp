/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Temple/constexpr/Numeric.h"

#include <iostream>

using namespace Scine;

BOOST_AUTO_TEST_CASE(numericAverageStdDev) {
  const std::vector<double> values {29, 30, 31, 32, 33};

  BOOST_CHECK(temple::average(values) == 31);
  BOOST_CHECK(
    std::fabs(temple::stddev(values) - std::sqrt(2))
    < 1e-10
  );
}
