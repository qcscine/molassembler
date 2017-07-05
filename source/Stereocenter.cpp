#include "Stereocenter.h"

#include "CNStereocenter.h"
#include "EZStereocenter.h"

namespace MoleculeManip {

std::basic_ostream<char>& operator << (
  std::basic_ostream<char>& os,
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
) {
  auto assignedOption = stereocenterPtr -> assigned();
  os << stereocenterPtr -> info();
  return os;
}

namespace Stereocenters {

bool strictComparePtr(
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& b
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

} // namespace Stereocenters

} // namespace MoleculeManip
