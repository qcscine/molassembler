#include "Stereocenter.h"

namespace MoleculeManip {

std::basic_ostream<char>& operator << (
  std::basic_ostream<char>& os,
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
) {
  auto assignedOption = stereocenterPtr -> assigned();
  os << stereocenterPtr -> type()
    << " at ";
  // call global namespace set output operator
  ::operator<<(os, stereocenterPtr -> involvedAtoms());
  os << ", assigned ";
  if((bool) assignedOption) os << assignedOption.value();
  else os << "u";
  os << "/" << stereocenterPtr -> assignments();
  return os;
}

} // eo namespace MoleculeManip
