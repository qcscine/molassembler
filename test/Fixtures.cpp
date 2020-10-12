/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Fixtures.h"

#include "boost/test/unit_test.hpp"

namespace Scine {
namespace Molassembler {

LowTemperatureFixture::LowTemperatureFixture() {
  BOOST_TEST_MESSAGE("Setting low temperature regime");
  swap();
}
LowTemperatureFixture::~LowTemperatureFixture() {
  BOOST_TEST_MESSAGE("Restoring previous temperature regime");
  swap();
}

} // namespace Molassembler
} // namespace Scine
