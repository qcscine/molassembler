/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "Fixtures.h"

#include "boost/test/unit_test.hpp"

namespace Scine {
namespace molassembler {

LowTemperatureFixture::LowTemperatureFixture() {
  BOOST_TEST_MESSAGE("Setting low temperature regime");
  swap();
}
LowTemperatureFixture::~LowTemperatureFixture() {
  BOOST_TEST_MESSAGE("Restoring previous temperature regime");
  swap();
}

} // namespace molassembler
} // namespace Scine
