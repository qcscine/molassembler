#include "Stereocenter.h"

#include "AtomStereocenter.h"
#include "BondStereocenter.h"

namespace molassembler {

namespace Stereocenters {

/* Static constants */
constexpr AtomIndexType Stereocenter::removalPlaceholder;

bool compareStereocenterEqual(
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& b
) {
  if(a -> type() == b -> type()) {
    if(a -> type() == Type::AtomStereocenter) {
      auto aDerived = std::dynamic_pointer_cast<AtomStereocenter>(a);
      auto bDerived = std::dynamic_pointer_cast<AtomStereocenter>(b);

      return *aDerived == *bDerived;
    }

    // Remaining case is BondStereocenter
    auto aDerived = std::dynamic_pointer_cast<BondStereocenter>(a);
    auto bDerived = std::dynamic_pointer_cast<BondStereocenter>(b);

    return *aDerived == *bDerived;
  }

  // Differing types are unequal
  return false;
}

bool compareStereocenterLessThan(
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& b
) {
  if(a -> type() < b -> type()) {
    return true;
  }

  if(a -> type() > b -> type()) {
    return false;
  }

  // Types are equal now
  if(a -> type() == Type::AtomStereocenter) {
    auto aDerived = std::dynamic_pointer_cast<AtomStereocenter>(a);
    auto bDerived = std::dynamic_pointer_cast<AtomStereocenter>(b);

    return *aDerived < *bDerived;
  }

  // Remaining case is BondStereocenter
  auto aDerived = std::dynamic_pointer_cast<BondStereocenter>(a);
  auto bDerived = std::dynamic_pointer_cast<BondStereocenter>(b);

  return *aDerived < *bDerived;
}

} // namespace Stereocenters

} // namespace molassembler
