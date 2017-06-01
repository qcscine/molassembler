#ifndef INCLUDE_CN_STEREOCENTER_H
#define INCLUDE_CN_STEREOCENTER_H

#include "Stereocenter.h"
#include "geometry_assignment/Assignment.h"
#include "symmetry_information/Symmetries.h"

using namespace std::string_literals;

namespace MoleculeManip {

namespace Stereocenters {

class CNStereocenter : public Stereocenter {
private:
/* Typedefs */
  using AssignmentType = UniqueAssignments::Assignment;

/* Private State */
  // Mapping between next neighbor atom index and symbolic ligand character
  std::map<AtomIndexType, char> _neighborCharMap;
  // Mapping between next neighbor atom index to permutational symmetry position
  std::map<AtomIndexType, unsigned> _neighborSymmetryPositionMap;
  // Vector of rotationally unique Assignments
  std::vector<AssignmentType> _uniqueAssignments;
  
/* Private members */
  /*!
   * Reduce substituent atoms at central atom to a mapping of their indices to
   * a symbolic ligand character. e.g.
   * map{4 -> 'C', 6 -> 'A', 9 -> 'A', 14 -> 'B'}
   */
  std::map<
    AtomIndexType,
    char
  > _reduceSubstituents(
    const std::vector<AtomIndexType>& rankedSubstituentNextAtomIndices,
    const std::set<
      std::pair<
        AtomIndexType,
        AtomIndexType
      >
    >& equalSubstituentPairsSet
  ) const;
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
/* Public state */
  Symmetry::Name symmetry;
  // Central atom of the Stereocenter, const on assignment
  const AtomIndexType centerAtom; 
  // The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> assignment;

/* Constructors */
  CNStereocenter(
    // The symmetry of this Stereocenter
    const Symmetry::Name& symmetry,
    // The atom this Stereocenter is centered on
    const AtomIndexType& centerAtom,
    // A partially ordered list of substituents, low to high
    const std::vector<AtomIndexType> partiallySortedSubstituents,
    // A set of pairs denoting which substituents are equal priority
    const std::set<
      std::pair<AtomIndexType, AtomIndexType>
    > equalPairsSet
  );

/* Modification */
  void assign(const unsigned& assignment) final;

  void changeSymmetry(const Symmetry::Name& symmetryName);

  void fit(const Delib::PositionCollection& positions) final;

/* Information */
  double angle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const final;

  boost::optional<unsigned> assigned() const final;

  unsigned numAssignments() const final;

  std::vector<ChiralityConstraintPrototype> chiralityConstraints() const final;

  std::vector<DihedralLimits> dihedralLimits() const final;

  std::string info() const final;

  std::set<AtomIndexType> involvedAtoms() const final;

  Type type() const final;

/* Operators */
  bool operator == (const CNStereocenter& other) const;

  friend std::basic_ostream<char>& MoleculeManip::operator << (
    std::basic_ostream<char>& os,
    const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
  );
};

} // namespace Stereocenters

} // namespace MoleculeManip

#endif
