/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEST_FIXTURES_H
#define INCLUDE_MOLASSEMBLER_TEST_FIXTURES_H

#include "Molassembler/Options.h"

namespace Scine {
namespace Molassembler {

/* Fixture enforcing low temperature approximation (no thermalization of atom
 * stereopermutators)
 */
struct LowTemperatureFixture {
  void swap() {
    std::swap(Options::Thermalization::pyramidalInversion, pyramidalInversion);
    std::swap(Options::Thermalization::berryPseudorotation, berryPseudorotation);
    std::swap(Options::Thermalization::bartellMechanism, bartellMechanism);
  }

  LowTemperatureFixture();
  ~LowTemperatureFixture();


  bool pyramidalInversion = false;
  bool berryPseudorotation = false;
  bool bartellMechanism = false;
};


} // namespace Molassembler
} // namespace Scine

#endif
