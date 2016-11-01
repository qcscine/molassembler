#ifndef INCLUDE_CN4_STEREOCENTER_H
#define INCLUDE_CN4_STEREOCENTER_H

#include <experimental/optional>

#include "ElementTypes.h"

#include "Stereocenter.h"
#include "Molecule.h"
#include "UniqueAssignments/Assignment.h"
#include "UniqueAssignments/SymmetryInformation.h"

namespace MoleculeManip {

namespace Stereocenters {

// Interestingly enough, the code is very general except for the constraint 
// collection, maybe some of that can be refactored into PermSymmetry...
class CN4Stereocenter : public Stereocenter {
private:
/* Typedefs */
  using AssignmentType = UniqueAssignments::Assignment<
    PermSymmetry::Tetrahedral
  >;

/* Private data */
  const Molecule* _molPtr;
  const AtomIndexType _centerAtom;

  //! State of whether it is assigned, and if so, in which 
  std::experimental::optional<unsigned> _assignment;

  //! Mapping between next neighbor atom index and symbolic ligand character
  std::map<
    AtomIndexType,
    char
  > _neighborCharMap;

  //! Mapping between next neighbor AtomIndexTypes to Permutational Symmetry positions
  std::map<AtomIndexType, uint8_t> _neighborSymmetryPositionMap;

  //! List of unique Assignments
  std::vector<AssignmentType> _uniqueAssignments;

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
   * Return a list of distance and chirality constraints. Generates 1-2 and 1-3
   * distance constraints for the central atom and its next neighbors. The 
   * target volume of the chirality constraint created by the tetrahedron is 
   * calculated using internal coordinates (the Cayley-Menger determinant),
   * always leading to V > 0, so depending on the current assignment, the sign 
   * of the result is switched. The formula used later in chirality constraint
   * calculation for explicit coordinates is adjusted by V' = 6 V to avoid an 
   * unnecessary factor, so we do that here too:
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
  virtual std::pair<
    std::vector<DistanceConstraint>,
    std::vector<ChiralityConstraint>
  > collectConstraints() const override final;

  /*!
   * Return the number of possible assignments at this feature
   */
  virtual unsigned assignments() const override final {
    return _uniqueAssignments.size();
  }

  /*!
   * Return whether this feature has been assigned or not
   */
  virtual std::experimental::optional<unsigned> assigned() const override final {
    return _assignment;
  }
};

} // eo namespace Stereocenters

}

#endif
