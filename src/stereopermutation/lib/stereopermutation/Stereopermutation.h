/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Base class to describe substituent arrangements in shapes
 *
 * Contains the base class employed for describing the particular manner in
 * which substituents are arranged in various shapes.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_STEREOPERMUTATION_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_STEREOPERMUTATION_H

#include "shapes/Shapes.h"
#include "temple/Preprocessor.h"
#include "temple/OperatorSuppliers.h"

#include "boost/optional/optional_fwd.hpp"

#include <set>
#include <vector>
#include <map>
#include <functional>

namespace Scine {

/**
 * @brief Permutations of possibly linked substituents around atom-centric
 *   shapes and combinations of two shapes
 */
namespace stereopermutation {

/*! @brief Represent abstract stereopermutation around atom center
 *
 * This class represents a simplified model of a sterically unique assignment
 * of a set of ligands to a stereocenter. It exists to uniquely identify the
 * steric configuration at this stereocenter, and provides methods to assist
 * a systematic generation of all possible configurations. It is generalized
 * over a number of shapes which are encoded in a separate library.
 */
class Stereopermutation : public temple::crtp::LexicographicComparable<Stereopermutation> {
public:
//!@name Member types
//!@{
  //! Type used to represent links between shape vertices
  using LinksSetType = std::set<
    std::pair<unsigned, unsigned>
  >;
//!@}

//!@name Static functions
//!@{
  /*! @brief Rotate characters
   *
   * @complexity{@math{\Theta(N)}}
   */
  static std::vector<char> rotateCharacters(
    const std::vector<char>& characters,
    const std::vector<unsigned>& rotation
  );

  /*! @brief Rotate links
   *
   * @complexity{@math{\Theta(L)}, not linear in shape size since it is small and constant}
   */
  static LinksSetType rotateLinks(
    const LinksSetType& links,
    const std::vector<unsigned>& rotation
  );
//!@}

//!@name Member data
//!@{
  //! Abstract representation of ranked substituents
  std::vector<char> characters;
  //! Links between characters
  LinksSetType links;
//!@}

//!@name Special member functions
//!@{
  // Do not default instantiate
  Stereopermutation() = delete;
  /*!
   * @brief Constructs an Stereopermutation from a list of ligand characters.
   *
   * \param passCharacters A vector of chars signifying abstract ligands.
   */
  Stereopermutation(std::vector<char> passCharacters);
  /*!
   * @brief Construct an Stereopermutation from a list of ligand characters and
   *   a list of bonded indices referencing the ligand characters.
   *
   * @param passCharacters A vector of chars signifying abstract ligands.
   * @param passLinks A vector of pairs. Describes which ligand characters
   *  are bonded to one another.
   */
  Stereopermutation(
    std::vector<char> passCharacters,
    LinksSetType passLinks
  );
//!@}

//!@name Modifiers
//!@{
  /*! @brief Applies a shape vertex rotation.
   *
   * @complexity{@math{O(N + L)}}
   */
  void applyRotation(const std::vector<unsigned>& rotation);
  /*! @brief Swap two "columns"
   *
   * @complexity{@math{\Theta(L)}}
   */
  void columnSwap(unsigned a, unsigned b);

  /*! @brief Transform this Stereopermutation into its lowest permutation.
   *
   * @complexity{@math{\Theta(1)} if already the lowest permutation,
   * @math{O(N!)} otherwise}
   */
  void lowestPermutation();

  /*! @brief Modify the "columns" to the previous permutation
   *
   * @complexity{@math{O(N / 2)}}
   */
  bool nextPermutation();

  /*! @brief Modify the "columns" to the previous permutation
   *
   * @complexity{@math{O(N / 2)}}
   */
  bool previousPermutation();

  /*! @brief Reverse a span of "columns"
   *
   * @complexity{@math{\Theta(N)}}
   */
  void reverseColumns(unsigned from, unsigned to);
//!@}

//!@name Information
//!@{
  /*! @brief Compare two "columns"
   *
   * Compares two "columns". An important note about how this works:
   * E.g. if we have: chars {A, A}, links {[0, 1]}, then columnSmaller(0, 1) is
   * false. These columns are considered equal size in order to avoid confusion
   * in the permutation code where the instruction columnSwap(0, 1) would have
   * no net effect.
   *
   * @complexity{@math{O(L \log L)}}
   */
  PURITY_WEAK bool columnSmaller(unsigned a, unsigned b) const;

  /*! @brief Generate all superimposable rotations of a stereopermutation
   *
   * Generates a set of all rotational equivalents of this Stereopermutation as
   * defined by its shape template parameter.
   *
   * @complexity{@math{O(\prod_i^Rm_i)} where @math{R} is the set of rotations and
   * @math{m_i} is the multiplicity of rotation @math{i}}
   */
  std::set<Stereopermutation> generateAllRotations(const Shapes::Shape& shape) const;

  /*!
   * Gets a map of ligand symbol character to position(s) in the permutational
   * shape.
   */
  std::map<
    char,
    std::vector<unsigned>
  > getCharMap() const;

  /*! @brief Returns whether the "columns" are sorted in ascending order
   *
   * @complexity{@math{O(N L \log L)}}
   */
  PURITY_WEAK bool isSortedAsc() const;

  /*!
   * @brief Checks whether a stereopermutation is a mirror image of another
   *   within a particular shape
   *
   * @complexity{As generateAllRotations}
   *
   * @returns boost::none If the shape does not generate enantiomers
   * @returns true If the shape has enantiomers, and this is the enantiomeric
   *   stereopermutation to @p other
   * @return false If the shape has enantiomers, and this is not the
   *   enantiomeric to @p other
   */
  boost::optional<bool> isEnantiomer(
    const Stereopermutation& other,
    const Shapes::Shape& shape
  ) const;

  /*!
   * Checks whether this Stereopermutation is rotationally superimposable with
   * another.
   *
   * @complexity{As generateAllRotations}
   */
  bool isRotationallySuperimposable(
    const Stereopermutation& other,
    const Shapes::Shape& shape
  ) const;

  /*!
   * Makes a set of a "column"'s connected indices.
   * If e.g. chars {A, A, A, A, A, A}, links {[0, 1], [1, 2]}, then
   * makeConnectedIndicesSet(1) = set {0, 2}.
   *
   * @complexity{@math{\Theta(L)}}
   */
  std::set<unsigned> makeConnectedIndicesSet(unsigned index) const;

  //! Converts positionOccupations into a string for display
  std::string toString() const;
//!@}

//!@name Operators
//!@{
  //! Yields members tied to tuple for crtp operator suppliers
  inline auto tie() const {
    return std::tie(characters, links);
  }
//!@}

private:
/* Private member functions */
  /*!
   * Implementation of the generation of a set of all rotational equivalents of
   * this Stereopermutation as defined by its shape template parameter. Takes an
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
    const Shapes::Shape shape
  ) const;
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

} // namespace Scine

#endif
