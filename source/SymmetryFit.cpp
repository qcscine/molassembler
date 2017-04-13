#include "SymmetryFit.h"

// Implementation
#include "template_magic/templateMagic.h"
#include "CommonTrig.h"
#include "VectorView.h"
#include <iomanip>

#include "DistanceGeometry/generateConformation.h"

namespace MoleculeManip {

/* Helper class Fit implementation */
double SymmetryFit::Fit::_calculateAngleDeviation(
  const std::vector<AtomIndexType> adjacentAtoms,
  const Delib::PositionCollection& positions,
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
) {
  return TemplateMagic::sum(
    TemplateMagic::allPairsMap(
      adjacentAtoms,
      [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
        return std::fabs(
          _getAngle( // The angle from the positions
            positions,
            i,
            CNStereocenterPtr -> centerAtom,
            k
          ) - _toRadians( // The ideal angle from the Stereocenter
            CNStereocenterPtr -> angle(
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

double SymmetryFit::Fit::_calculateOneThreeDeviation(
  const std::vector<AtomIndexType> adjacentAtoms,
  const Delib::PositionCollection& positions,
  std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
) {
  return TemplateMagic::sum(
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
            _toRadians(
              CNStereocenterPtr -> angle( // idealized Stereocenter angle
                i,
                CNStereocenterPtr -> centerAtom,
                k
              )
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

  if(prototypes.size() == 0) {
    return 0;
  } else {
    /* Propagate all prototypes defined by the stereocenter into constraints
     * using distances from the positions
     */
    const auto constraints = TemplateMagic::map(
      prototypes,
      DistanceGeometry::detail::makePropagator(
        [&positions](const unsigned& i, const unsigned& j) {
          return (
            positions[i] - positions[j]
          ).norm();
        }
      )
    );

    /* Calculate deviations for each constraint between the chirality constraint
     * and the volume of the tetrahedron from the positions. Since the target 
     * volume is calculated using the same distances and the two quantities are
     * therefore different only in sign, the deviation is either zero or two times
     * the volume.
     */
    const auto deviations = TemplateMagic::map(
      constraints,
      [&](const auto& constraint) -> double {
        return std::fabs(
          constraint.target - _getVolume(
            positions,
            constraint.indices[0],
            constraint.indices[1],
            constraint.indices[2],
            constraint.indices[3]
          )
        );
      }
    );

    // Return the summation of the deviations
    return TemplateMagic::sum(deviations);
  }
}

SymmetryFit::Fit::Fit(
  const Symmetry::Name& symmetryName,
  const unsigned& assignment,
  const std::vector<AtomIndexType> adjacentAtoms,
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

double SymmetryFit::Fit::_toRadians(const double& inDegrees) {
  return M_PI * inDegrees / 180;
}


double SymmetryFit::Fit::totalDeviation() const {
  return angleDeviation + oneThreeDeviation + chiralityDeviation;
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
  const Delib::PositionCollection& positions
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
    ) continue;

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

} // eo namespace

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
