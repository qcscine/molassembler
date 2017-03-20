#include "Stereocenter.h"

namespace MoleculeManip {

std::basic_ostream<char>& operator << (
  std::basic_ostream<char>& os,
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
) {
  auto assignedOption = stereocenterPtr -> assigned();
  os << stereocenterPtr -> info();
  return os;
}

} // eo namespace MoleculeManip
