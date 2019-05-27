/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Types.h"

#include "boost/functional/hash.hpp"

#include <tuple>

namespace Scine {

namespace molassembler {

BondIndex::BondIndex() = default;

BondIndex::BondIndex(AtomIndex a, AtomIndex b) noexcept : first(a), second(b) {
  if(second < first) {
    std::swap(first, second);
  }
}

bool BondIndex::contains(const AtomIndex a) const {
  return a == first || a == second;
}

bool BondIndex::operator < (const BondIndex& other) const {
  return std::tie(first, second) < std::tie(other.first, other.second);
}

bool BondIndex::operator == (const BondIndex& other) const {
  return std::tie(first, second) == std::tie(other.first, other.second);
}

std::size_t hash_value(const BondIndex& bond) {
  std::size_t seed = 0;
  boost::hash_combine(seed, bond.first);
  boost::hash_combine(seed, bond.second);
  return seed;
}

} // namespace molassembler

} // namespace Scine
