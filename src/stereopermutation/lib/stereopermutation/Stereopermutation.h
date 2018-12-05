// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_STEREOPERMUTATION_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_STEREOPERMUTATION_H

#include "chemical_symmetries/Names.h"
#include "temple/Preprocessor.h"

#include "boost/optional/optional_fwd.hpp"

#include <set>
#include <vector>
#include <map>
#include <functional>

/*! @file
 *
 * @brief Base class to describe substituent arrangements in symmetries
 *
 * Contains the base class employed for describing the particular manner in
 * which substituents are arranged in various symmetries.
 */

namespace stereopermutation {

/*!
 * This class represents a simplified model of a sterically unique assignment
 * of a set of ligands to a stereocenter. It exists to uniquely identify the
 * steric configuration at this stereocenter, and provides methods to assist
 * a systematic generation of all possible configurations. It is generalized
 * over a number of symmetries which are encoded in a separate library.
 */
struct Stereopermutation {
private:
/* Private member functions */
  /*!
   * Implementation of the generation of a set of all rotational equivalents of
   * this Stereopermutation as defined by its symmetry template parameter. Takes an
   * interrupt callback as an argument to which it passes *this and a new
   * rotational structure every time one is found. If the callback returns
   * true, the generation of assignments is terminated and a pair containing
   * the set of generated assignments and a boolean with the value true is
   * returned. If the generation is allowed to finish, the full set and the
   * boolean false are returned.
   */
  std::pair<std::set<Stereopermutation>, bool> _generateAllRotations(
    const std::function<
      bool(const Stereopermutation&, const Stereopermutation&)
    >& interruptCallbackOnNewStereopermutation,
    const Symmetry::Name& symmetryName
  ) const;

public:
//!@name Member types
//!@{
  using LinksSetType = std::set<
    std::pair<unsigned, unsigned>
  >;
//!@}

//!@name Static functions
//!@{
  //! Rotate charactes according to template symmetry
  static std::vector<char> rotateCharacters(
    const std::vector<char>& characters,
    const std::vector<unsigned>& rotationIndices
  );

  //! Rotate links according to template symmetry
  static LinksSetType rotateLinks(
    const LinksSetType& links,
    const std::vector<unsigned>& rotationIndices
  );
//!@}

//!@name Member data
//!@{
  std::vector<char> characters;
  LinksSetType links;
//!@}

//!@name Special member functions
//!@{
  // Do not default instantiate
  Stereopermutation() = delete;
  /*!
   * Constructs an Stereopermutation from a list of ligand characters.
   *
   * \param passSymmetryName The name of the employed symmetry.
   * \param passCharacters A vector of chars signifying abstract ligands.
   *
   * \throws If size of passSymmetryName does not match the number of
   *   characters in debug builds
   */
  Stereopermutation(
    Symmetry::Name passSymmetryName [[gnu::unused]],
    std::vector<char> passCharacters
  );
  /*!
   * Construct an Stereopermutation from a list of ligand characters and a list of
   * bonded indices referencing the ligand characters.
   *
   * \param passSymmetryName The name of the employed symmetry.
   * \param passCharacters A vector of chars signifying abstract ligands.
   * \param passLinks A vector of pairs. Describes which ligand characters
   *  are bonded to one another.
   *
   * \throws If size of passSymmetryName does not match the number of
   *   characters in debug builds
   */
  Stereopermutation(
    Symmetry::Name passSymmetryName [[gnu::unused]],
    std::vector<char> passCharacters,
    LinksSetType passLinks
  );
//!@}

//!@name Modifiers
//!@{
  //! Swap two "columns"
  void columnSwap(unsigned a, unsigned b);

  //! Transform this Stereopermutation into its lowest permutation.
  void lowestPermutation();

  //! Modify the "columns" to the previous permutation
  bool nextPermutation();

  //! Modify the "columns" to the previous permutation
  bool previousPermutation();

  //! Reverse a span of "columns"
  void reverseColumns(unsigned from, unsigned to);
//!@}

//!@name Information
//!@{
  //! Applies a Symmetry rotation.
  void applyRotation(const std::vector<unsigned>& rotationIndices);

  /* Information –––––––––––––––––––––––– */
  /*!
   * Compares two "columns". An important note about how this works:
   * E.g. if we have: chars {A, A}, links {[0, 1]}, then columnSmaller(0, 1) is
   * false. These columns are considered equal size in order to avoid confusion
   * in the permutation code where the instruction columnSwap(0, 1) would have
   * no net effect.
   */
  PURITY_WEAK bool columnSmaller(unsigned a, unsigned b) const;

  /*!
   * Generates a set of all rotational equivalents of this Stereopermutation as
   * defined by its symmetry template parameter.
   */
  std::set<Stereopermutation> generateAllRotations(const Symmetry::Name& symmetryName) const;

  /*!
   * Gets a map of ligand symbol character to position(s) in the permutational
   * symmetry.
   */
  std::map<
    char,
    std::vector<unsigned>
  > getCharMap() const;

  //! Returns whether the "columns" are sorted in ascending order
  PURITY_WEAK bool isSortedAsc() const;

  /*!
   * @brief Checks whether a stereopermutation is a mirror image of another
   *   within a particular symmetry
   *
   * @returns boost::none If the symmetry does not generate enantiomers
   * @returns true If the symmetry has enantiomers, and this is the enantiomeric
   *   stereopermutation to @p other
   * @return false If the symmetry has enantiomers, and this is not the
   *   enantiomeric to @p other
   */
  boost::optional<bool> isEnantiomer(
    const Stereopermutation& other,
    const Symmetry::Name& symmetryName
  ) const;

  /*!
   * Checks whether this Stereopermutation is rotationally superimposable with
   * another.
   */
  bool isRotationallySuperimposable(
    const Stereopermutation& other,
    const Symmetry::Name& symmetryName
  ) const;

  /*!
   * Makes a set of a "column"'s connected indices.
   * If e.g. chars {A, A, A, A, A, A}, links {[0, 1], [1, 2]}, then
   * makeConnectedIndicesSet(1) = set {0, 2}.
   */
  std::set<unsigned> makeConnectedIndicesSet(unsigned index) const;

  //! Converts positionOccupations into a string for display
  std::string toString() const;
//!@}

//!@name Operators
//!@{
  PURITY_WEAK bool operator < (const Stereopermutation& other) const;
  PURITY_WEAK bool operator > (const Stereopermutation& other) const;
  PURITY_WEAK bool operator == (const Stereopermutation& other) const;
  PURITY_WEAK bool operator != (const Stereopermutation& other) const;
//!@}
};

PURITY_WEAK std::size_t hash_value(const Stereopermutation& assignment);

/*!
 * ostream operator for easier debugging
 */
std::ostream& operator << (
  std::ostream& os,
  const Stereopermutation& a
);

} // namespace stereopermutation

#endif
