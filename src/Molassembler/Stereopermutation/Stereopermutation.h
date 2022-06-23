/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Base class to describe substituent arrangements in shapes
 *
 * Contains the base class employed for describing the particular manner in
 * which substituents are arranged in various shapes.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_STEREOPERMUTATION_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_STEREOPERMUTATION_H

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Temple/StrongIndexPermutation.h"

namespace Scine {
namespace Molassembler {

//! @brief Data classes for permutational spatial arrangement modeling
namespace Stereopermutations {

//! Helper tag to differentiate index types for site rank
struct rank_tag;
//! Ranked index of a site
using Rank = Temple::StrongIndex<rank_tag, unsigned>;

/*! @brief Represent abstract stereopermutation around atom center
 *
 * This class represents a simplified model of a sterically unique assignment
 * of a set of ligands to a stereocenter. It exists to uniquely identify the
 * steric configuration at this stereocenter, and provides methods to assist
 * a systematic generation of all possible configurations. It is generalized
 * over a number of shapes which are encoded in a separate library.
 */
class MASM_EXPORT Stereopermutation : public Temple::Crtp::LexicographicComparable<Stereopermutation> {
public:
//!@name Member types
//!@{
  //! Ranked site occupation
  using Occupation = Temple::StrongIndexPermutation<Shapes::Vertex, Rank>;

  //! Type used to represent a link between shape vertices
  using Link = std::pair<Shapes::Vertex, Shapes::Vertex>;

  //! Unordered type, but kept ordered by member functions
  using OrderedLinks = std::vector<Link>;
//!@}

//!@name Static functions
//!@{
  /*! @brief Rotate characters
   *
   * @complexity{@math{\Theta(N)}}
   */
  static Occupation permuteOccupation(
    const Occupation& occupation,
    const Shapes::Permutation& permutation
  );

  static Occupation occupationFromChars(const std::string& chars);

  /*! @brief Rotate links
   *
   * @complexity{@math{\Theta(L)}, not linear in shape size since it is small and constant}
   */
  static OrderedLinks permuteLinks(
    const OrderedLinks& links,
    const Shapes::Permutation& permutation
  );
//!@}

//!@name Member data
//!@{
  //! Occupation of shape vertices by site rank
  Occupation occupation;
  //! Links between shape vertices
  OrderedLinks links;
//!@}

//!@name Special member functions
//!@{
  //! Empty default-init (invalid state)
  Stereopermutation() = default;
  /*!
   * @brief Construct an Stereopermutation from an occupation and a
   *   a list of bonded shape vertices
   *
   * @param passOccupation An occupation of shape vertices by site rankings
   * @param passLinks A vector of shape vertex pairs denoting links
   */
  Stereopermutation(
    Occupation passOccupation,
    OrderedLinks passLinks = {}
  );
//!@}

//!@name Information
//!@{
  /*! @brief Applies a shape vertex permutation.
   *
   * @complexity{@math{O(N + L)}}
   */
  Stereopermutation applyPermutation(const Shapes::Permutation& permutation) const;

  //! A string for display
  std::string toString() const;
//!@}

//!@name Operators
//!@{
  //! Yields members tied to tuple for crtp operator suppliers
  inline auto tie() const {
    return std::tie(occupation, links);
  }
//!@}
};

PURITY_WEAK MASM_EXPORT std::size_t hash_value(const Stereopermutation& assignment);

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine

#endif
