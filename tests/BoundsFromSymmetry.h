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

// Essentially mimics BFSConstraintVisitor in its functionality
/* MoleculeManip::DistanceGeometry::DistanceBoundsMatrix distanceBoundsFromSymmetry(
  const Symmetry::Name& symmetry,
  const DistancesOption& distancesOption
) {
  const unsigned N = Symmetry::size(symmetry) + 1; // Central atom needed
  MoleculeManip::DistanceGeometry::DistanceBoundsMatrix bounds(N);
  double variance = 0.025;

  switch(distancesOption) {
    case (DistancesOption::Uniform): {
      // Fill in 1-2 bounds to center
      for(unsigned i = 0; i < Symmetry::size(symmetry); i++) {
        bounds.upperBound(0, 1 + i) = 1.5 + variance;
        bounds.lowerBound(0, 1 + i) = 1.5 - variance;
      }
    } break;
    case (DistancesOption::Incremented): {
      std::vector<double> distanceTo(N - 1);

      std::iota(
        distanceTo.begin(),
        distanceTo.end(),
        1.0
      );

      double trueDistance;

      // Fill in 1-2 bounds to center
      for(unsigned i = 0; i < N - 1; i++) {
        trueDistance = 1 + distanceTo[i] / (N - 1);
        bounds.upperBound(0, 1 + i) = trueDistance + variance;
        bounds.lowerBound(0, 1 + i) = trueDistance - variance;
      }
    } break;
    case (DistancesOption::Random): {
      // Fill in 1-2 bounds to center
      double randomDistance;
      for(unsigned i = 0; i < N - 1; i++) {
        randomDistance = RNG::rng.getSingle<double>(1, 2);
        bounds.upperBound(0, 1 + i) = randomDistance + variance;
        bounds.lowerBound(0, 1 + i) = randomDistance - variance;
      }
    } break;
  }

  // Generate 1-3 distance bounds
  for(unsigned i = 0; i < N - 1; i++) {
    for(unsigned j = i + 1; j < N - 1; j++) {
      bounds.upperBound(1 + i, 1 + j) = MoleculeManip::CommonTrig::lawOfCosines(
        bounds.upperBound(0, 1 + i),
        bounds.upperBound(0, 1 + j),
        Symmetry::angleFunction(symmetry)(i, j) * M_PI / 180.0
      ) + variance;

      bounds.lowerBound(1 + i, 1 + j) = MoleculeManip::CommonTrig::lawOfCosines(
        bounds.lowerBound(0, 1 + i),
        bounds.lowerBound(0, 1 + j),
        Symmetry::angleFunction(symmetry)(i, j) * M_PI / 180.0
      ) - variance;
    }
  }

  return bounds;
} */

}

#endif
