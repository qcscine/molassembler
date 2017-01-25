#ifndef LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H
#define LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H

#include <vector>
#include <map>
#include <set>

#include "symmetry_information/Symmetries.h"

/* TODO
 * - update documentation
 */

/* NOTES
 * - One path to implementation success of implementing permutations may be to
 *   implement iterators on Assignment that will then allow 
 *   std::next_permutation to run on Assignment! Positions, length, and when 
 *   characters are "different" must be clearly defined so as to mimic a 
 *   container 
 */

namespace UniqueAssignments {

/*!
 * This class represents a simplified model of a sterically unique assignment
 * of a set of ligands to a stereocenter. It exists to uniquely identify the 
 * steric configuration at this stereocenter, and provides methods to assist
 * a systematic generation of all possible configurations. It is generalized 
 * over a number of symmetries which are encoded elsewhere and that serve as 
 * template parameter to this class.
 * \tparam Symmetry The template class specifying which symmetry is to be 
 *  enforced.
 */
struct Assignment {
private:
/* Private member functions */
  /*!
   * Implementation of the generation of a set of all rotational equivalents of
   * this Assignment as defined by its symmetry template parameter. Takes an 
   * interrupt callback as an argument to which it passes *this and a new 
   * rotational structure every time one is found. If the callback returns
   * true, the generation of assignments is terminated and a pair containing
   * the set of generated assignments and a boolean with the value true is 
   * returned. If the generation is allowed to finish, the full set and the 
   * boolean false are returned.
   */
  std::pair<
    std::set<Assignment>,
    bool
  > _generateAllRotations(
    std::function<
      bool(const Assignment&, const Assignment&)
    > interruptCallbackOnNewAssignment
  ) const;

public:
/* Typedefs */
  using LinksSetType = std::set<
    std::pair<unsigned, unsigned>
  >;

/* Public members */
  const Symmetry::Name symmetryName;
  std::vector<char> characters;
  LinksSetType links;

  /* Constructors */
  // Do not default instantiate
  Assignment() = delete;
  /*!
   * Constructs an Assignment from a list of ligand characters.
   * \tparam Symmetry A SymmetryInformation derived class template.
   * \param characters A vector of chars signifying abstract ligands.
   */
  Assignment(
    const Symmetry::Name& passSymmetryName,
    const std::vector<char>& passCharacters
  );
  /*!
   * Construct an Assignment from a list of ligand characters and a list of 
   * bonded indices referencing the ligand characters.
   * \tparam Symmetry A SymmetryInformation derived class template.
   * \param characters A vector of chars signifying abstract ligands.
   * \param pairedIndices A vector of pairs. Describes which ligand characters 
   *  are bonded to one another.
   */
  Assignment(
    const Symmetry::Name& passSymmetryName,
    const std::vector<char>& passCharacters,
    const LinksSetType& passLinks
  );

  /* Modifiers ––––––––––––––––––––––––––––*/
  //! Swap two "columns"
  void columnSwap(const unsigned& a, const unsigned& b);

  //! Transform this Assignment into its lowest permutation.
  void lowestPermutation();

  //! Modify the "columns" to the previous permutation
  bool nextPermutation();

  //! Modify the "columns" to the previous permutation
  bool previousPermutation();

  //! Reverse a span of "columns"
  void reverseColumns(const unsigned& from, const unsigned& to);

  //! Rotate charactes according to template symmetry
  std::vector<char> rotateCharacters(
    const std::vector<char>& characters,
    const unsigned& rotationFunctionIndex
  );

  //! Rotate links according to template symmetry
  LinksSetType rotateLinks(
    const LinksSetType& links,
    const unsigned& rotationFunctionIndex
  );

  /*!
   * Applies a Symmetry rotation.
   */
  void applyRotation(const unsigned& index);

  /* Information –––––––––––––––––––––––– */
  /*!
   * Compares two "columns". An important note about how this works:
   * E.g. if we have: chars {A, A}, links {[0, 1]}, then columnSmaller(0, 1) is
   * false. These columns are considered equal size in order to avoid confusion
   * in the permutation code where the instruction columnSwap(0, 1) would have
   * no net effect.
   */
  bool columnSmaller(const unsigned& a, const unsigned& b) const;


  /*!
   * Generates a set of all rotational equivalents of this Assignment as 
   * defined by its symmetry template parameter.
   */
  std::set<Assignment> generateAllRotations() const;

  /*!
   * Gets a map of ligand symbol character to position in the permutational 
   * symmetry. 
   */
  std::map<
    char,
    std::vector<unsigned>
  > getCharMap() const;

  //! Returns whether the "columns" are sorted in ascending order
  bool isSortedAsc() const;

  /*!
   * Checks whether this Assignment is rotationally superimposable with
   * another.
   */
  bool isRotationallySuperimposable(const Assignment& other) const;

  /*!
   * Makes a set of a "column"'s connected indices.
   * If e.g. chars {A, A, A, A, A, A}, links {[0, 1], [1, 2]}, then
   * makeConnectedIndicesSet(1) = set {0, 2}.
   */
  std::set<unsigned> makeConnectedIndicesSet(const unsigned& index) const;

  //! Converts positionOccupations into a string for display
  std::string toString() const;

  /* Operators */
  bool operator < (const Assignment& other) const;
  bool operator > (const Assignment& other) const;
  bool operator == (const Assignment& other) const;
  bool operator != (const Assignment& other) const;
};

/*!
 * ostream operator for easier debugging
 */
std::ostream& operator << (
  std::ostream& os,
  const Assignment& a
);

} // eo namespace

#endif
