#include "SymmetryFit.h"

// Implementation
#include "template_magic/TemplateMagic.h"
#include "CommonTrig.h"
#include <iomanip>

#include "DistanceGeometry/generateConformation.h"
#include <Eigen/Geometry>

/* TODO
 * - Some sanity checks ought to be beneficial, i.e. fitting to seesaw on a 
 *   tetravalent carbon center is abject nonsense and wasted effort
 * - Optimize -> No need to carry all fits in memory! Keeping track of the 
 *   (perhaps several) lowest is enough.
 */

namespace MoleculeManip {

/* Helper class Fit implementation */
double SymmetryFit::Fit::_calculateAngleDeviation(
  const std::vector<AtomIndexType>& adjacentAtoms,
  const Delib::PositionCollection& positions,
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
) {
  return TemplateMagic::numeric::sum(
    TemplateMagic::allPairsMap(
      adjacentAtoms,
      [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
        return std::fabs(
          _getAngle( // The angle from the positions
            positions,
            i,
            CNStereocenterPtr -> centerAtom,
            k
          ) - CNStereocenterPtr -> angle(
            i,
            CNStereocenterPtr -> centerAtom,
            k
          )
        );
      }
    )
  );
}

double SymmetryFit::Fit::_calculateOneThreeDeviation(
  const std::vector<AtomIndexType>& adjacentAtoms,
  const Delib::PositionCollection& positions,
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
) {
  return TemplateMagic::numeric::sum(
    TemplateMagic::allPairsMap(
      adjacentAtoms,
      [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
        return std::fabs(
          _getDistance( // i-k 1-3 distance from positions
            positions,
            i,
            k
          ) - CommonTrig::lawOfCosines( // idealized 1-3 distance from
            _getDistance( // i-j 1-2 distance from positions
              positions,
              i,
              CNStereocenterPtr -> centerAtom
            ),
            _getDistance( // j-k 1-2 distance from positions
              positions,
              CNStereocenterPtr -> centerAtom,
              k
            ),
            CNStereocenterPtr -> angle( // idealized Stereocenter angle
              i,
              CNStereocenterPtr -> centerAtom,
              k
            )
          )
        );
      }
    )
  );
}

double SymmetryFit::Fit::_calculateChiralityDeviation(
  const Delib::PositionCollection& positions,
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
) {
  const auto prototypes = CNStereocenterPtr -> chiralityConstraints();

  if(prototypes.empty()) { // early return
    return 0;
  } 

  return TemplateMagic::numeric::sum(
    TemplateMagic::map(
      prototypes,
      [&positions](const auto& constraintPrototype) -> double {
        using TargetEnumType = Stereocenters::ChiralityConstraintTarget;

        double volume = _getVolume(
          positions,
          constraintPrototype.indices[0],
          constraintPrototype.indices[1],
          constraintPrototype.indices[2],
          constraintPrototype.indices[3]
        );

        // If the target is flat, then the "error" is continuous:
        if(constraintPrototype.target == TargetEnumType::Flat) {
          return std::fabs(volume);
        }

        /* Otherwise, no bounds, error is some arbitrary penalty if the sign is
         * wrong
         */
        if(
          (
            constraintPrototype.target == TargetEnumType::Positive
            && volume < 0
          ) || (
            constraintPrototype.target == TargetEnumType::Negative
            && volume > 0
          )
        ) {
          return 1; // Arbitrary penalty
        }

        return 0;
      }
    )
  );
}

SymmetryFit::Fit::Fit(
  const Symmetry::Name& symmetryName,
  const unsigned& assignment,
  const std::vector<AtomIndexType>& adjacentAtoms,
  const Delib::PositionCollection& positions,
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
) : symmetryName(symmetryName),
    assignment(assignment) 
{
  // Calculate the deviation from the positions
  /* Three contributions to deviations
   * - i-j-k angles
   * - 1-3 distances (via 1-2 bond distances from positions and 
   *   symmetry-ideal angles)
   * - chirality constraints (if applicable)
   */
  angleDeviation = _calculateAngleDeviation(
    adjacentAtoms,
    positions,
    CNStereocenterPtr
  );

  oneThreeDeviation = _calculateOneThreeDeviation(
    adjacentAtoms,
    positions,
    CNStereocenterPtr
  );

  chiralityDeviation = _calculateChiralityDeviation(
    positions,
    CNStereocenterPtr
  );
}

double SymmetryFit::Fit::_getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k
) {
  auto a = positions[i].asEigenVector() - positions[j].asEigenVector(),
       b = positions[k].asEigenVector() - positions[j].asEigenVector();

  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

double SymmetryFit::Fit::_getDistance(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j
) {
  return (
    positions[i].asEigenVector()
    - positions[j].asEigenVector()
  ).norm();
}

double SymmetryFit::Fit::_getVolume(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
) {
  return (
    positions[i].asEigenVector()
    - positions[l].asEigenVector()
  ).dot(
    (
      positions[j].asEigenVector()
      - positions[l].asEigenVector()
    ).cross(
      positions[k].asEigenVector()
      - positions[l].asEigenVector()
    )
  );
}

double SymmetryFit::Fit::totalDeviation() const {
  return angleDeviation + oneThreeDeviation + chiralityDeviation + symmetryPenalty;
}

bool SymmetryFit::Fit::operator < (const Fit& other) const {
  return TemplateMagic::componentSmaller(
    totalDeviation(),
    other.totalDeviation()
  ).value_or(
    assignment < other.assignment
  );
}

/* SymmetryFit implementation */
SymmetryFit::SymmetryFit(
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr,
  const std::vector<AtomIndexType>& adjacentAtoms,
  const Delib::PositionCollection& positions,
  const boost::optional<Symmetry::Name>& expectedSymmetry
) {
  auto initialAssignment = CNStereocenterPtr -> assigned();
  auto initialSymmetry = CNStereocenterPtr -> symmetry;

  /* STEPS
   * - Reduce the local geometry to some combination of 
   *   internal coordinates, 1-3 distances and signed tetrahedron values
   *   (chirality constraints) so that it can be compared easily to gathered
   *   distance geometry constraints
   * - Cycle through all Symmetries of appropriate size, checking the 
   *   gathered 3D data against supposed DG constraints of the candidate
   *   Symmetry and assignment (if present)
   * - Cycle through the Stereocenter's assignments and calculate total
   *   deviation for all angles
   * - Pick the assignment that has lowest deviation, if it has multiplicity 1
   */

  for(const auto& symmetryName : Symmetry::allNames) {
    if( // Skip any Symmetries of different size
      Symmetry::size(symmetryName) != Symmetry::size(
        CNStereocenterPtr -> symmetry
      )
    ) {
      continue;
    }

    // Change the symmetry of the CNStereocenter
    CNStereocenterPtr -> changeSymmetry(symmetryName);

    std::vector<double> deviations;

    for(
      unsigned assignment = 0;
      assignment < (CNStereocenterPtr -> numAssignments());
      assignment++
    ) {
      // Assign the stereocenter
      CNStereocenterPtr -> assign(assignment);

      // Do the fit
      Fit currentFit {
        symmetryName,
        assignment,
        adjacentAtoms,
        positions,
        CNStereocenterPtr
      };

      if(
        expectedSymmetry 
        && expectedSymmetry.value() != symmetryName
      ) {
        currentFit.symmetryPenalty = 0.5;
      }

      // Add to fits vector
      _fits.insert(currentFit);

    }
  }

  // Return the CNStereocenter to its initial state
  CNStereocenterPtr -> changeSymmetry(initialSymmetry);
  if(initialAssignment) {
    CNStereocenterPtr -> assign(initialAssignment.value());
  }

  lowestDeviation = _fits.begin() -> totalDeviation();
  bestSymmetry = _fits.begin() -> symmetryName;

  for(const auto& fit : _fits) {
    if(fit.totalDeviation() == lowestDeviation) {
      /* We assume that equally good best fits must be from the same symmetry,
       * fail if this is ever untrue
       */
      assert(fit.symmetryName == bestSymmetry);
      
      assignmentsWithLowestDeviation.push_back(fit.assignment);
    }
  }
} 

} // namespace MoleculeManip

// Global namespace ostream operator
std::ostream& operator << (std::ostream& os, const MoleculeManip::SymmetryFit& symmetryFit) {
  /* Make a functor satisfying std::set::Compare that sorts by symmetry Index 
   * and then by assignment number
   */
  struct Sorter {
    bool operator() (
      const MoleculeManip::SymmetryFit::Fit& a,
      const MoleculeManip::SymmetryFit::Fit& b
    ) {
      return TemplateMagic::componentSmaller<unsigned>(
        Symmetry::nameIndex(a.symmetryName),
        Symmetry::nameIndex(b.symmetryName)
      ).value_or(
        a.assignment < b.assignment
      );
    }
  };

  // Resort the fit set using the previous functor
  std::set<MoleculeManip::SymmetryFit::Fit, Sorter> resortedSet (
    symmetryFit._fits.begin(),
    symmetryFit._fits.end()
  );

  // Output
  for(const auto& fit : resortedSet) {
    os << Symmetry::nameIndex(fit.symmetryName)
      << ", " << fit.assignment
      << ", " << std::setprecision(4) << std::fixed
      << fit.angleDeviation << ", "
      << fit.oneThreeDeviation << ", "
      << fit.chiralityDeviation 
      << std::endl;
  }

  return os;
}
