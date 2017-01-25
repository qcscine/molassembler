#ifndef INCLUDE_CN4_STEREOCENTER_H
#define INCLUDE_CN4_STEREOCENTER_H


#include "ElementTypes.h"

#include <boost/optional.hpp>
#include <Eigen/Core>

#include "Stereocenter.h"
#include "steric_uniqueness/Assignment.h"

namespace MoleculeManip {

// forward-declare Molecule to avoid dependency
class Molecule;

namespace Stereocenters {

// Interestingly enough, the code is very general except for the constraint 
// collection, maybe some of that can be refactored into PermSymmetry...
class CN4Stereocenter : public Stereocenter {
private:
/* Typedefs */
  using AssignmentType = UniqueAssignments::Assignment;

/* Private data */
  const Molecule* _molPtr;
  const AtomIndexType _centerAtom;

  //! State of whether it is assigned, and if so, in which 
  boost::optional<unsigned> _assignment;

  //! Mapping between next neighbor atom index and symbolic ligand character
  std::map<
    AtomIndexType,
    char
  > _neighborCharMap;

  /*! 
   * Mapping between next neighbor AtomIndexTypes to Permutational Symmetry
   * positions
   */
  std::map<AtomIndexType, unsigned> _neighborSymmetryPositionMap;

  //! List of unique Assignments
  std::vector<AssignmentType> _uniqueAssignments;

  /*! 
   * If distanceConstraints is called before chiralityConstraints, then 
   * this option has a value
   */
  mutable boost::optional<
    Eigen::Matrix<double, 5, 5>
  > cayleyMengerOption;

/* Private member functions */
  /*!
   * Reduce substituent atoms at central atom to a mapping of their indices to
   * a symbolic ligand character. e.g.
   * map{4 -> 'C', 6 -> 'A', 9 -> 'A', 14 -> 'B'}
   */
  std::map<
    AtomIndexType,
    char
  > _reduceSubstituents() const;

  /*!
   * Reduce a mapping of atom indices to symbolic ligand characters to a 
   * character vector
   */
  std::vector<char> _reduceNeighborCharMap(
    const std::map<
      AtomIndexType,
      char
    >& neighborCharMap
  );

public:
/* Public member functions */
  /* Construction */
  CN4Stereocenter(
    const Molecule* molPtr,
    const AtomIndexType& center
  );

  /* Modification */
  /*!
   * Assign this feature
   */
  virtual void assign(const unsigned& assignment) override final; 

  /* Information */
  /*!
   * Return a string specifying the type of stereocenter
   */
  virtual std::string type() const override final {
    return "CN4Stereocenter";
  }

  /*!
   * Return a set of involved atom indices
   */
  virtual std::set<AtomIndexType> involvedAtoms() const override final {
    return {_centerAtom};
  }

  /*!
   * Return a list of distance constraints. Generates 1-2 and 1-3 distance
   * constraints for the central atom and its next neighbors. 
   */
  virtual std::vector<DistanceConstraint> distanceConstraints() const override final;

  /*!
   * Return a list of chirality constraints.  The target volume of the
   * chirality constraint created by the tetrahedron is calculated using
   * internal coordinates (the Cayley-Menger determinant), always leading to V
   * > 0, so depending on the current assignment, the sign of the result is
   * switched. The formula used later in chirality constraint calculation for
   * explicit coordinates is adjusted by V' = 6 V to avoid an unnecessary
   * factor, so we do that here too:
   *               
   *    288 V²  = |...|               | substitute V' = 6 V
   * -> 8 (V')² = |...|               
   * ->      V' = sqrt(|...| / 8)
   *
   * where the Cayley-Menger determinant |...| is square symmetric:
   *   
   *          |   0   1   1   1   1 |
   *          |       0 d12 d13 d14 |
   *  |...| = |           0 d23 d24 |
   *          |               0 d34 |
   *          |  ...              0 |
   *
   */
  virtual std::vector<ChiralityConstraint> chiralityConstraints() const override final;

  /*!
   * Return the number of possible assignments at this feature
   */
  virtual unsigned assignments() const override final {
    return _uniqueAssignments.size();
  }

  /*!
   * Return whether this feature has been assigned or not
   */
  virtual boost::optional<unsigned> assigned() const override final {
    return _assignment;
  }
};

} // eo namespace Stereocenters

}

#endif
