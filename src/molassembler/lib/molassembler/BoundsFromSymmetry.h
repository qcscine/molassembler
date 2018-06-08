#ifndef INCLUDE_TESTING_DISTANCE_BOUNDS_FROM_SYMMETRY_H
#define INCLUDE_TESTING_DISTANCE_BOUNDS_FROM_SYMMETRY_H

#include "Molecule.h"

/*! @file
 *
 * Some testing helper functions to generate prototypical molecules
 */

namespace DGDBM { // Distance Geometry Distance Bounds Matrix

enum class DistancesOption {
  Uniform, // all 1.5
  Incremented, // all in 1/Nth steps from 1 to 2 ( N is symmetry size )
  Random // randomly chosen between 1 and 2
};

/*!
 * Constructs a prototypical symmetric (i.e. all ligands identical) molecule
 * with a Ruthenium center and Hydrogen substituents otherwise.
 */
[[deprecated]]
molassembler::Molecule symmetricMolecule(
  const Symmetry::Name& symmetry
);

molassembler::Molecule symmetricMolecule(
  const Symmetry::Name& symmetry
) {
  using namespace molassembler;

  Molecule molecule(
    Delib::ElementType::Ru,
    Delib::ElementType::H,
    BondType::Single
  );

  while(molecule.numAtoms() - 1 < Symmetry::size(symmetry)) {
    molecule.addAtom(
      Delib::ElementType::H,
      0,
      BondType::Single
    );
  }

  return molecule;
}

/*!
 * Constructs a prototypical asymmetric (i.e. all ligands different) molecule
 * of a specific size which is not VSEPR compliant past size 5.
 *
 * Does not set the central stereocenter symmetry.
 */
[[deprecated]]
molassembler::Molecule asymmetricMolecule(
  const Symmetry::Name& symmetry
);
molassembler::Molecule asymmetricMolecule(
  const Symmetry::Name& symmetry
) {
  const std::vector<Delib::ElementType> elements {
    Delib::ElementType::F,
    Delib::ElementType::Cl,
    Delib::ElementType::Br,
    Delib::ElementType::I,
    Delib::ElementType::N,
    Delib::ElementType::C,
    Delib::ElementType::O,
    Delib::ElementType::S,
    Delib::ElementType::P
  };

  using namespace molassembler;
  Molecule molecule(
    Delib::ElementType::Ru,
    Delib::ElementType::H,
    BondType::Single
  );

  for(
    unsigned elementIndex = 0;
    molecule.numAtoms() - 1 < Symmetry::size(symmetry);
    elementIndex++
  ) {
    molecule.addAtom(
      elements.at(elementIndex),
      0,
      BondType::Single
    );
  }

  return molecule;
}

} // namespace DGDBM

#endif
