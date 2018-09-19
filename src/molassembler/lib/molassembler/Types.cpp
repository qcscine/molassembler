// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Types.h"

#include <tuple>

namespace molassembler {

BondIndex::BondIndex() = default;
BondIndex::BondIndex(AtomIndex a, AtomIndex b) noexcept : first(a), second(b) {
  if(b < a) {
    std::swap(a, b);
  }
}

bool BondIndex::operator < (const BondIndex& other) const {
  return std::tie(first, second) < std::tie(other.first, other.second);
}

bool BondIndex::operator == (const BondIndex& other) const {
  return std::tie(first, second) == std::tie(other.first, other.second);
}

} // namespace molassembler
