#include "molassembler/Patterns.h"

#include "molassembler/Molecule.h"

namespace Scine {
namespace molassembler {
namespace patterns {

PlugType methyl() {
  Molecule mol {
    Utils::ElementType::C,
    Utils::ElementType::H
  };

  for(unsigned i = 0; i < 2; ++i) {
    mol.addAtom(Utils::ElementType::H, 0);
  }

  return {mol, {0}};
}

} // namespace patterns
} // namespace molassembler
} // namespace Scine
