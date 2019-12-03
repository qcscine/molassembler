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
