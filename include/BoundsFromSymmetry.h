#ifndef INCLUDE_TESTING_DISTANCE_BOUNDS_FROM_SYMMETRY
#define INCLUDE_TESTING_DISTANCE_BOUNDS_FROM_SYMMETRY

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "symmetry_information/Symmetries.h"
#include "RNG/RNG.h"
#include "CommonTrig.h"
#include "Molecule.h"
#include "DistanceGeometry/BFSConstraintCollector.h"

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

MoleculeManip::DistanceGeometry::DistanceBoundsMatrix distanceBoundsFromSymmetry(
  const Symmetry::Name& symmetry,
  const DistancesOption& distancesOption __attribute__ ((unused))
) __attribute__ ((deprecated));

MoleculeManip::DistanceGeometry::DistanceBoundsMatrix distanceBoundsFromSymmetry(
  const Symmetry::Name& symmetry,
  const DistancesOption& distancesOption __attribute__ ((unused))
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

  return molecule.getDistanceBoundsMatrix();
}

}

#endif
