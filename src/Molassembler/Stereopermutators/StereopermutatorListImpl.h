/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Private implementation of StereopermutatorList
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATOR_LIST_IMPL_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATOR_LIST_IMPL_H

#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/BondStereopermutator.h"

#include <unordered_map>

namespace Scine {
namespace Molassembler {

struct StereopermutatorList::Impl {
  using AtomMapType = std::unordered_map<AtomIndex, AtomStereopermutator>;
  using BondMapType = std::unordered_map<BondIndex, BondStereopermutator, boost::hash<BondIndex>>;

//!@name Modification
//!@{
  /*! @brief Access atom stereopermutator at an atom
   *
   * @throws std::out_of_range If there is no atom stereopermutator on the
   * passed atom
   *
   * @complexity{@math{\Theta(1)}}
   */
  AtomStereopermutator& at(AtomIndex index);

  /*! @brief Access bond stereopermutator at a bond
   *
   * @throws std::out_of_range If there is no bond stereopermutator on the
   * passed bond
   *
   * @complexity{@math{\Theta(1)}}
   */
  BondStereopermutator& at(const BondIndex& index);

  /*! @brief Add a new AtomStereopermutator to the list
   *
   * @complexity{@math{\Theta(1)} amortized}
   *
   * @returns An iterator pointing to the added stereopermutator
   */
  AtomStereopermutator& add(AtomStereopermutator stereopermutator);

  /*! @brief Add a new BondStereopermutator to the list
   *
   * @complexity{@math{\Theta(1)} amortized}
   *
   * @returns An iterator pointing to the added stereopermutator
   */
  BondStereopermutator& add(BondStereopermutator stereopermutator);

  /*! @brief Apply an index mapping to the list of stereopermutators
   *
   * Applies the permutation to its maps, transforming the keys (atom and bond
   * indices) and all stereopermutators.
   *
   * @complexity{@math{\Theta(A + B)}}
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*! @brief Remove all stereopermutators
   *
   * @complexity{@math{\Theta(1)}}
   */
  void clear();

  /*! @brief Remove all stereopermutators on bonds
   *
   * @complexity{@math{\Theta(1)}}
   */
  void clearBonds();

  /*! @brief Fetch a reference-option to an AtomStereopermutator
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<AtomStereopermutator&> option(AtomIndex index);

  /*! @brief Fetch a reference-option to a BondStereopermutator
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<BondStereopermutator&> option(const BondIndex& edge);

  /*! @brief Communicates the removal of a vertex index to all stereopermutators in the list
   *
   * Removing a vertex invalidates some vertex descriptors, which are used
   * liberally in the stereopermutator classes. This function ensures that
   * vertex descriptors are valid throughout.
   *
   * @complexity{@math{\Theta(A + B)}}
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  /*! @brief Removes the AtomStereopermutator on a specified index, if present
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool remove(AtomIndex index);

  /*! @brief Removes the BondStereopermutator on a specified edge, if present
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool remove(const BondIndex& edge);
//!@}

//!@name Information
//!@{
  /*! @brief Access atom stereopermutator at an atom
   *
   * @throws std::out_of_range If there is no atom stereopermutator on the
   * passed atom
   *
   * @complexity{@math{\Theta(1)}}
   */
  const AtomStereopermutator& at(AtomIndex index) const;

  /*! @brief Access bond stereopermutator at a bond
   *
   * @throws std::out_of_range If there is no bond stereopermutator on the
   * passed bond
   *
   * @complexity{@math{\Theta(1)}}
   */
  const BondStereopermutator& at(const BondIndex& index) const;

  /*! @brief Modular comparison with another StereopermutatorList using a bitmask
   *
   * @complexity{@math{O(A + B)}}
   */
  bool compare(
    const Impl& other,
    AtomEnvironmentComponents componentBitmask
  ) const;

  /*! @brief Returns true if there are no stereopermutators
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool empty() const;

  /*! @brief Returns true if there are any stereopermutators with zero possible assignments
   *
   * @complexity{@math{O(A + B)}}
   */
  bool hasZeroAssignmentStereopermutators() const;

  /*! @brief Returns true if there are unassigned stereopermutators
   *
   * @complexity{@math{O(A + B)}}
   */
  bool hasUnassignedStereopermutators() const;

  /*! @brief Fetch a const ref-option to an AtomStereopermutator, if present
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<const AtomStereopermutator&> option(AtomIndex index) const;

  /*! @brief Fetch a const ref-option to a BondStereopermutator, if present
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<const BondStereopermutator&> option(const BondIndex& edge) const;

  /*! @brief Returns the number of AtomStereopermutators
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned A() const;

  /*! @brief Returns the number of BondStereopermutators
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned B() const;

  /*! @brief Combined size of atom and bond-stereopermutator lists
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned size() const;
//!@}

//!@name Operators
//!@{
  /*! @brief Strict equality comparison
   *
   * @complexity{@math{O(A + B)}}
   */
  bool operator == (const Impl& other) const;
//!@}

  //! The underlying storage for atom stereopermutators
  AtomMapType atomStereopermutators;

  //! The underlying storage for bond stereopermutators
  BondMapType bondStereopermutators;
};

} // namespace Molassembler
} // namespace Scine
#endif
