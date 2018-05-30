#include "Stereocenter.h"

#include "CNStereocenter.h"
#include "EZStereocenter.h"

namespace molassembler {

namespace Stereocenters {

/* Static constants */
constexpr AtomIndexType Stereocenter::removalPlaceholder;

bool compareStereocenterEqual(
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& b
) {
  if(a -> type() == b -> type()) {
    if(a -> type() == Type::CNStereocenter) {
      auto aDerived = std::dynamic_pointer_cast<CNStereocenter>(a);
      auto bDerived = std::dynamic_pointer_cast<CNStereocenter>(b);

      return *aDerived == *bDerived;
    }

    // Remaining case is EZStereocenter
    auto aDerived = std::dynamic_pointer_cast<EZStereocenter>(a);
    auto bDerived = std::dynamic_pointer_cast<EZStereocenter>(b);

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
  if(a -> type() == Type::CNStereocenter) {
    auto aDerived = std::dynamic_pointer_cast<CNStereocenter>(a);
    auto bDerived = std::dynamic_pointer_cast<CNStereocenter>(b);

    return *aDerived < *bDerived;
  }

  // Remaining case is EZStereocenter
  auto aDerived = std::dynamic_pointer_cast<EZStereocenter>(a);
  auto bDerived = std::dynamic_pointer_cast<EZStereocenter>(b);

  return *aDerived < *bDerived;
}

} // namespace Stereocenters

} // namespace molassembler
