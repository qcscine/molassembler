/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEST_FIXTURES_H
#define INCLUDE_MOLASSEMBLER_TEST_FIXTURES_H

#include "molassembler/Options.h"

namespace Scine {
namespace Molassembler {

/* Fixture enforcing low temperature approximation (no thermalization of atom
 * stereopermutators)
 */
struct LowTemperatureFixture {
  void swap() { std::swap(Options::temperatureRegime, priorRegime); }

  LowTemperatureFixture();
  ~LowTemperatureFixture();

  TemperatureRegime priorRegime = TemperatureRegime::Low;
};


} // namespace Molassembler
} // namespace Scine

#endif
