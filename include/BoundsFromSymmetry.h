#ifndef INCLUDE_TESTING_DISTANCE_BOUNDS_FROM_SYMMETRY_H
#define INCLUDE_TESTING_DISTANCE_BOUNDS_FROM_SYMMETRY_H

#include "symmetry_information/Symmetries.h"
#include "Molecule.h"

namespace DGDBM { // Distance Geometry Distance Bounds Matrix

enum class DistancesOption {
  Uniform, // all 1.5
  Incremented, // all in 1/Nth steps from 1 to 2 ( N is symmetry size )
  Random // randomly chosen between 1 and 2
};

MoleculeManip::Molecule symmetricMolecule(
  const Symmetry::Name& symmetry
) {
  using namespace MoleculeManip;

  Molecule molecule(
    Delib::ElementType::Ru,
    Delib::ElementType::H,
    BondType::Single
  );

  while(molecule.getNumAtoms() - 1 < Symmetry::size(symmetry)) {
    molecule.addAtom(
      Delib::ElementType::H,
      0,
      BondType::Single
    );
  }

  auto rankResultPair = molecule.getAdjacencyList().rankPriority(0);

  auto stereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
    symmetry,
    0,
    rankResultPair.first,
    rankResultPair.second
  );

  stereocenterPtr -> assign(0);

  molecule.stereocenters.add(stereocenterPtr);

  return molecule;
}

MoleculeManip::Molecule asymmetricMolecule(
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

  using namespace MoleculeManip;
  Molecule molecule(
    Delib::ElementType::Ru,
    Delib::ElementType::H,
    BondType::Single
  );

  for(
    unsigned elementIndex = 0;
    molecule.getNumAtoms() - 1 < Symmetry::size(symmetry);
    elementIndex++
  ) {
    molecule.addAtom(
      elements.at(elementIndex),
      0,
      BondType::Single
    );
  }

  auto rankResultPair = molecule.getAdjacencyList().rankPriority(0);

  auto stereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
    symmetry,
    0,
    rankResultPair.first,
    rankResultPair.second
  );

  stereocenterPtr -> assign(0);

  molecule.stereocenters.add(stereocenterPtr);

  return molecule;
}

}

#endif
