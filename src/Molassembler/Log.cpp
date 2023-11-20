/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Log.h"

namespace Scine {
namespace Molassembler {
namespace Log {
namespace Detail {

int NullBuffer::overflow(int c) { return c; }

NullBuffer nullBuffer {};
std::ostream nullStream(&nullBuffer);

} // namespace Detail

std::ostream& log(const Level& decisionLevel) {
  if(decisionLevel >= level) {
    return std::cout;
  }

  return Detail::nullStream;
}

std::ostream& log(const Particulars& particular) {
  if(particulars.count(particular) == 1) {
    return std::cout;
  }

  return Detail::nullStream;
}

bool isSet(const Particulars particular) {
  return particulars.count(particular) > 0;
}

Level level = Level::Trace;
std::unordered_set<Particulars> particulars {};

} // namespace Log
} // namespace Molassembler
} // namespace Scine
