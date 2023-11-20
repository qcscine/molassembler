/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Owning class storing all stereopermutators in a molecule
 *
 * Contains the declaration for a class that stores a list of all stereopermutators
 * in a molecule.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATOR_LIST_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATOR_LIST_H

#include "Molassembler/Types.h"
#include "Molassembler/IteratorRange.h"

#include "boost/optional/optional_fwd.hpp"
#include <memory>
#include <vector>

namespace Scine {
namespace Molassembler {

class AtomStereopermutator;
class BondStereopermutator;

/**
 * @brief Manages all stereopermutators that are part of a Molecule
 */
class MASM_EXPORT StereopermutatorList {
private:
  struct Impl;

public:
//!@name Public types
//!@{
  template<typename Permutator>
  class iterator {
  public:
    static_assert(
      std::is_same<std::decay_t<Permutator>, AtomStereopermutator>::value
      || std::is_same<std::decay_t<Permutator>, BondStereopermutator>::value,
      "This type may not be instantiated for any type other than atom and bond stereopermutators"
    );

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Permutator;
    using pointer = Permutator*;
    using reference = Permutator&;

    iterator(iterator&& other) noexcept;
    iterator& operator = (iterator&& other) noexcept;
    iterator(const iterator& other);
    iterator& operator = (const iterator& other);
    ~iterator();

    iterator();

    using ImplType = std::conditional_t<
      std::is_const<Permutator>::value,
      const StereopermutatorList::Impl&,
      StereopermutatorList::Impl&
    >;
    iterator(ImplType impl, bool begin);

    iterator& operator ++ ();
    iterator operator ++ (int);
    reference operator * () const;

    bool operator == (const iterator& other) const;
    bool operator != (const iterator& other) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };

  using AtomStereopermutatorIterator = iterator<AtomStereopermutator>;
  using AtomStereopermutatorConstIterator = iterator<const AtomStereopermutator>;
  using BondStereopermutatorIterator = iterator<BondStereopermutator>;
  using BondStereopermutatorConstIterator = iterator<const BondStereopermutator>;
//!@}

//!@name Special member functions
//!@{
  StereopermutatorList();
  StereopermutatorList(StereopermutatorList&& other) noexcept;
  StereopermutatorList& operator = (StereopermutatorList&& other) noexcept;
  StereopermutatorList(const StereopermutatorList& other);
  StereopermutatorList& operator = (const StereopermutatorList& other);
  ~StereopermutatorList();
//!@}

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
   * @returns A reference to the added stereopermutator
   */
  AtomStereopermutator& add(AtomStereopermutator stereopermutator);

  /*! @brief Add a new BondStereopermutator to the list
   *
   * @complexity{@math{\Theta(1)} amortized}
   *
   * @returns A reference to the added stereopermutator
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
   * @returns Whether a permutator was deleted
   */
  bool remove(AtomIndex index);

  /*! @brief Removes the BondStereopermutator on a specified edge, if present
   *
   * @complexity{@math{\Theta(1)}}
   * @returns Whether a permutator was deleted
   */
  bool remove(const BondIndex& edge);

  /*! @brief Removes the AtomStereopermutator on a specified index, if present
   *
   * @complexity{@math{\Theta(1)}}
   */
  [[deprecated("Prefer remove")]]
  void try_remove(AtomIndex index);

  /*! @brief Removes the BondStereopermutator on a specified edge, if present
   *
   * @complexity{@math{\Theta(1)}}
   */
  [[deprecated("Prefer remove")]]
  void try_remove(const BondIndex& edge);
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
    const StereopermutatorList& other,
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

//!@name Ranges (not thread-safe)
//!@{
  /*! @brief Returns an iterable object with modifiable atom stereopermutator references
   *
   * @complexity{@math{\Theta(1)}}
   */
  IteratorRange<AtomStereopermutatorIterator> atomStereopermutators();

  /*! @brief Returns an iterable object with unmodifiable atom stereopermutator references
   *
   * @complexity{@math{\Theta(1)}}
   */
  IteratorRange<AtomStereopermutatorConstIterator> atomStereopermutators() const;

  /*! @brief Returns an iterable object with modifiable bond stereopermutator references
   *
   * @complexity{@math{\Theta(1)}}
   */
  IteratorRange<BondStereopermutatorIterator> bondStereopermutators();

  /*! @brief Returns an iterable object with unmodifiable bond stereopermutator references
   *
   * @complexity{@math{\Theta(1)}}
   */
  IteratorRange<BondStereopermutatorConstIterator> bondStereopermutators() const;
//!@}

//!@name Operators
//!@{
  /*! @brief Strict equality comparison
   *
   * @complexity{@math{O(A + B)}}
   */
  bool operator == (const StereopermutatorList& other) const;
  //! Inverts @p operator ==
  bool operator != (const StereopermutatorList& other) const;
//!@}

private:
  std::unique_ptr<Impl> impl_;
};

} // namespace Molassembler
} // namespace Scine
#endif
