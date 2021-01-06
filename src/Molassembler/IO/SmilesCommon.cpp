/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/IO/SmilesCommon.h"
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace IO {
namespace {

template<typename ... Ts>
constexpr unsigned valenceFillHandleDifferences(Ts ... is) {
  constexpr unsigned N = sizeof...(is);
  const std::array<int, N> values {is...};

  for(unsigned i = 0; i < N; ++i) {
    int x = values[i];

    if(x == 0) {
      return 0;
    }

    if(x > 0) {
      return x;
    }
  }

  return 0;
}

static_assert(valenceFillHandleDifferences(-1, 1) == 1, "Nope");
static_assert(valenceFillHandleDifferences(-1, 0) == 0, "Nope");
static_assert(valenceFillHandleDifferences(-3, -1) == 0, "Nope");
static_assert(valenceFillHandleDifferences(-3, -1, 1) == 1, "Nope");

} // namespace

bool isValenceFillElement(const Utils::ElementType e) {
  const unsigned Z = Utils::ElementInfo::Z(e);
  if(5 <= Z && Z <= 9) {
    // B, C, N, O, F
    return true;
  }

  if(15 <= Z && Z <= 17) {
    // P, S, Cl
    return true;
  }

  if(Z == 35 || Z == 53) {
    // Br, I
    return true;
  }

  return false;
}

unsigned valenceFillElementImplicitHydrogenCount(
  int valence,
  Utils::ElementType e,
  const bool aromatic
) {
  assert(valence >= 0);
  assert(isValenceFillElement(e));

  if(aromatic) {
    valence += 1;
  }

  /* Quoting from the spec:
   *
   * The implicit hydrogen count is determined by summing the bond orders of
   * the bonds connected to the atom. If that sum is equal to a known valence
   * for the element or is greater than any known valence then the implicit
   * hydrogen count is 0. Otherwise the implicit hydrogen count is the
   * difference between that sum and the next highest known valence.
   */

  switch(Utils::ElementInfo::Z(e)) {
    case 5: return std::max(0, 3 - valence); // B
    case 6: return std::max(0, 4 - valence); // C
    case 7: { // N
      return valenceFillHandleDifferences(3 - valence, 5 - valence);
    }
    case 8: return std::max(0, 2 - valence); // O
    case 15: { // P
      return valenceFillHandleDifferences(3 - valence, 5 - valence);
    }
    case 16: { // S
      return valenceFillHandleDifferences(2 - valence, 4 - valence, 6 - valence);
    }
    default: return std::max(0, 1 - valence); // F, Cl, Br, I are the remaining cases
  }
}

} // namespace IO
} // namespace Molassembler
} // namespace Scine
