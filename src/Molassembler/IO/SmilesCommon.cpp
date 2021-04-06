/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/IO/SmilesCommon.h"
#include "Molassembler/Modeling/BondDistance.h"
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
  switch(e) {
    case Utils::ElementType::B:
    case Utils::ElementType::C:
    case Utils::ElementType::N:
    case Utils::ElementType::O:
    case Utils::ElementType::F:
    case Utils::ElementType::P:
    case Utils::ElementType::S:
    case Utils::ElementType::Cl:
    case Utils::ElementType::Br:
    case Utils::ElementType::I:
      return true;
    default:
      return false;
  }
}

int vertexValence(const PrivateGraph::Vertex i, const PrivateGraph& g) {
  int valence = 0;
  for(const PrivateGraph::Edge edge : g.edges(i)) {
    valence += Bond::bondOrderMap.at(
      static_cast<unsigned>(g.bondType(edge))
    );
  }
  return valence;
}

unsigned valenceFillElementImplicitHydrogenCount(
  int valence,
  Utils::ElementType e
) {
  assert(valence >= 0);
  assert(isValenceFillElement(e));

  /* Quoting from the spec:
   *
   * The implicit hydrogen count is determined by summing the bond orders of
   * the bonds connected to the atom. If that sum is equal to a known valence
   * for the element or is greater than any known valence then the implicit
   * hydrogen count is 0. Otherwise the implicit hydrogen count is the
   * difference between that sum and the next highest known valence.
   */

  switch(Utils::ElementInfo::base(e)) {
    case Utils::ElementType::B: return std::max(0, 3 - valence);
    case Utils::ElementType::C: return std::max(0, 4 - valence);
    case Utils::ElementType::N: {
      return valenceFillHandleDifferences(3 - valence, 5 - valence);
    }
    case Utils::ElementType::O: return std::max(0, 2 - valence);
    case Utils::ElementType::P: {
      return valenceFillHandleDifferences(3 - valence, 5 - valence);
    }
    case Utils::ElementType::S: {
      return valenceFillHandleDifferences(2 - valence, 4 - valence, 6 - valence);
    }
    // F, Cl, Br, I are the remaining cases, handled identically
    default: return std::max(0, 1 - valence);
  }
}

} // namespace IO
} // namespace Molassembler
} // namespace Scine
