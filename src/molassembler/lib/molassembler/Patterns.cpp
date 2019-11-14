/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Patterns.h"

#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"

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

Molecule alkane(const unsigned N) {
  if(N <= 1) {
    // Methane
    Molecule methane {Utils::ElementType::C, Utils::ElementType::H};
    methane.addAtom(Utils::ElementType::H, 0);
    methane.addAtom(Utils::ElementType::H, 0);
    methane.addAtom(Utils::ElementType::H, 0);
    methane.setShapeAtAtom(0, Shapes::Shape::Tetrahedron);
    return methane;
  }

  Molecule mol {Utils::ElementType::C, Utils::ElementType::C};
  constexpr unsigned firstTerminus = 0;
  unsigned secondTerminus = 1;
  while(mol.graph().N() != N) {
    secondTerminus = mol.addAtom(Utils::ElementType::C, secondTerminus);
  }

  for(unsigned i = 0; i <= secondTerminus; ++i) {
    // Add two hydrogens to each carbon
    mol.addAtom(Utils::ElementType::H, i);
    mol.addAtom(Utils::ElementType::H, i);
  }

  // Add one more hydrogen to the terminii
  mol.addAtom(Utils::ElementType::H, firstTerminus);
  mol.addAtom(Utils::ElementType::H, secondTerminus);

  // Set the shapes at all carbons
  for(unsigned i = 0; i <= secondTerminus; ++i) {
    mol.setShapeAtAtom(0, Shapes::Shape::Tetrahedron);
  }

  return mol;
}

} // namespace patterns
} // namespace molassembler
} // namespace Scine
