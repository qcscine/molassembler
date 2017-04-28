#include "Stereocenter.h"

#include "template_magic/TemplateMagic.h"
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
  using namespace MoleculeManip::Stereocenters;

  if(a -> type() == b -> type()) {
    if(a -> type() == Type::CNStereocenter) {
      auto aDerived = std::dynamic_pointer_cast<CNStereocenter>(a);
      auto bDerived = std::dynamic_pointer_cast<CNStereocenter>(b);

      return *aDerived == *bDerived;
    } else { // EZStereocenter
      auto aDerived = std::dynamic_pointer_cast<EZStereocenter>(a);
      auto bDerived = std::dynamic_pointer_cast<EZStereocenter>(b);

      return *aDerived == *bDerived;
    }
  } else {
    return false;
  }
}

} // namespace Stereocenters

} // namespace MoleculeManip
