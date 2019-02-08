/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Owning class storing all stereopermutators in a molecule
 *
 * Contains the declaration for a class that stores a list of all stereopermutators
 * in a molecule.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATOR_LIST_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATOR_LIST_H

#include "boost/range/adaptor/map.hpp"
#include "boost/functional/hash.hpp"

#include "temple/constexpr/Bitmask.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"

#include <unordered_map>

namespace Scine {

namespace molassembler {

/**
 * @brief Manages all stereopermutators that are part of a Molecule
 */
class StereopermutatorList {
public:
/* Typedefs */
  using AtomMapType = std::unordered_map<AtomIndex, AtomStereopermutator>;
  using BondMapType = std::unordered_map<BondIndex, BondStereopermutator, boost::hash<BondIndex>>;

/* Modification */
  //! Add a new AtomStereopermutator to the list
  AtomMapType::iterator add(AtomStereopermutator stereopermutator);

  //! Add a new BondStereopermutator to the list
  BondMapType::iterator add(BondStereopermutator stereopermutator);

  //! Apply an index mapping to the list of stereopermutators
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  //! Remove all stereopermutators
  void clear();

  //! Remove all stereopermutators on bonds
  void clearBonds();

  //! Fetch a reference-option to an AtomStereopermutator, if present
  boost::optional<AtomStereopermutator&> option(AtomIndex index);

  //! Fetch a reference-option to a BondStereopermutator, if present
  boost::optional<BondStereopermutator&> option(const BondIndex& edge);

  /*! Communicates the removal of a vertex index to all stereopermutators in the list
   *
   * Removing a vertex invalidates some vertex descriptors, which are used
   * liberally in the stereopermutator classes. This function ensures that
   * vertex descriptors are valid throughout.
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  //! Removes the AtomStereopermutator on a specified index
  void remove(AtomIndex index);

  //! Removes the BondStereopermutator on a specified edge
  void remove(const BondIndex& edge);

  //! Removes the AtomStereopermutator on a specified index, if present
  void try_remove(AtomIndex index);

  //! Removes the BondStereopermutator on a specified edge, if present
  void try_remove(const BondIndex& edge);

/* Information */
  //! Modular comparison with another StereopermutatorList using a bitmask
  bool compare(
    const StereopermutatorList& other,
    AtomEnvironmentComponents componentBitmask
  ) const;

  //! Returns true if there are no stereopermutators
  bool empty() const;

  //! Returns true if there are any stereopermutators with zero possible assignments
  bool hasZeroAssignmentStereopermutators() const;

  //! Returns true if there are unassigned stereopermutators
  bool hasUnassignedStereopermutators() const;

  //! Fetch a const ref-option to an AtomStereopermutator, if present
  boost::optional<const AtomStereopermutator&> option(AtomIndex index) const;

  //! Fetch a const ref-option to a BondStereopermutator, if present
  boost::optional<const BondStereopermutator&> option(const BondIndex& edge) const;

  //! Returns the number of AtomStereopermutators
  unsigned A() const;

  //! Returns the number of BondStereopermutators
  unsigned B() const;

  //! Combined size of atom and bond-stereopermutator lists
  unsigned size() const;

/* Iterators */
  //! Returns an iterable object with modifiable atom stereopermutator references
  boost::range_detail::select_second_mutable_range<AtomMapType> atomStereopermutators();

  //! Returns an iterable object with unmodifiable atom stereopermutator references
  boost::range_detail::select_second_const_range<AtomMapType> atomStereopermutators() const;

  //! Returns an iterable object with modifiable bond stereopermutator references
  boost::range_detail::select_second_mutable_range<BondMapType> bondStereopermutators();

  //! Returns an iterable object with unmodifiable bond stereopermutator references
  boost::range_detail::select_second_const_range<BondMapType> bondStereopermutators() const;

/* Operators */
  //! Strict equality comparison
  bool operator == (const StereopermutatorList& other) const;
  //! Inverts @p operator ==
  bool operator != (const StereopermutatorList& other) const;

private:
  //! The underlying storage for atom stereopermutators
  AtomMapType _atomStereopermutators;

  //! The underlying storage for bond stereopermutators
  BondMapType _bondStereopermutators;
};

} // namespace molassembler

} // namespace Scine
#endif
