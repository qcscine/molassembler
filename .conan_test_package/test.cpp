/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Types.h"

using namespace Scine::Molassembler;

int main() {
  BondIndex f {4, 3};

  if(f.first == 3) {
    return 0;
  }

  return 1;
}
