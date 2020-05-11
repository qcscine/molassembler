/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Types.h"

#include "boost/functional/hash.hpp"

#include <tuple>
#include <iterator>

namespace Scine {
namespace Molassembler {

static_assert(
  static_cast<std::underlying_type<BondType>::type>(BondType::Eta) == nBondTypes - 1,
  "Did you add a bond type and not alter nBondTypes?"
);

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

BondIndex::const_iterator BondIndex::begin() const {
  return &first;
}

//! Returns the address past second
BondIndex::const_iterator BondIndex::end() const {
  return std::next(&second);
}

std::size_t hash_value(const BondIndex& bond) {
  std::size_t seed = 0;
  boost::hash_combine(seed, bond.first);
  boost::hash_combine(seed, bond.second);
  return seed;
}

} // namespace Molassembler
} // namespace Scine
