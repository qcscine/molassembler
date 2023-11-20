/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/constexpr/Numeric.h"

#include <iostream>

using namespace Scine::Molassembler;

BOOST_AUTO_TEST_CASE(NumericAverageStdDev, *boost::unit_test::label("Temple")) {
  const std::vector<double> values {29, 30, 31, 32, 33};

  BOOST_CHECK(Temple::average(values) == 31);
  BOOST_CHECK(
    std::fabs(Temple::stddev(values) - std::sqrt(2))
    < 1e-10
  );
}
