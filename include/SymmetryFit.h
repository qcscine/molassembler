#ifndef INCLUDE_SYMMETRY_FIT_H
#define INCLUDE_SYMMETRY_FIT_H

#include <map>
#include "CNStereocenter.h"
#include "Types/PositionCollection.h"

// Forward-declarations so global-scope operator can be friended
namespace MoleculeManip { class SymmetryFit; }
std::ostream& operator << (std::ostream& os, const MoleculeManip::SymmetryFit& fit);


namespace MoleculeManip {

class SymmetryFit {
private:
/* Helper classes */
  class Fit {
  private:
    double _calculateAngleDeviation(
      const std::vector<AtomIndexType> adjacentAtoms,
      const Delib::PositionCollection& positions,
      std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
    );

    double _calculateOneThreeDeviation(
      const std::vector<AtomIndexType> adjacentAtoms,
      const Delib::PositionCollection& positions,
      std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
    );

  public:
    Symmetry::Name symmetryName;
    unsigned assignment;

    double angleDeviation;
    double oneThreeDeviation;
    double chiralityDeviation;

    Fit() = delete;
    Fit(
      const Symmetry::Name& symmetryName,
      const unsigned& assignment,
      const std::vector<AtomIndexType> adjacentAtoms,
      const Delib::PositionCollection& positions,
      std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr
    );

    static double _getAngle(
      const Delib::PositionCollection& positions,
      const AtomIndexType& i,
      const AtomIndexType& j,
      const AtomIndexType& k
    );

    static double _getDistance(
      const Delib::PositionCollection& positions,
      const AtomIndexType& i,
      const AtomIndexType& j
    );

    static inline double _toRadians(const double& inDegrees);

    double totalDeviation() const;

    bool operator < (const Fit& other) const;
  };
  
/* Private members */
  std::set<Fit> _fits;

/* Private member functions */
public:
  /* Public members */
  Symmetry::Name bestSymmetry;
  std::vector<unsigned> assignmentsWithLowestDeviation;
  double lowestDeviation;

  /* Constructor */
  SymmetryFit(
    std::shared_ptr<Stereocenters::CNStereocenter>& CNStereocenterPtr,
    const std::vector<AtomIndexType>& adjacentAtoms,
    const Delib::PositionCollection& positions
  ); 

  friend std::ostream& (::operator << )(std::ostream& os, const SymmetryFit& fit);
};

} // eo namespace

#endif
