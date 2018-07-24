#ifndef INCLUDE_MOLASSEMBLER_STEREOCENTER_LIST_H
#define INCLUDE_MOLASSEMBLER_STEREOCENTER_LIST_H

#include "boost/range/adaptor/map.hpp"

#include "temple/constexpr/Bitmask.h"

#include "AtomStereocenter.h"
#include "BondStereocenter.h"

#include <unordered_map>
#include <map>

/*! @file
 *
 * Contains the declaration for a class that stores a list of all stereocenters
 * in a molecule.
 */

/* TODO
 * - Make a hash for GraphType::edge_descriptor so we can use unordered_map
 * - bond stereocenter state propagation
 */

namespace molassembler {

class StereocenterList {
public:
/* Typedefs */
  using AtomMapType = std::unordered_map<AtomIndexType, AtomStereocenter>;
  using BondMapType = std::map<GraphType::edge_descriptor, BondStereocenter>;

/* Modification */
  //! Add a new AtomStereocenter to the list
  void add(AtomIndexType i, AtomStereocenter stereocenter);

  //! Add a new BondStereocenter to the list
  void add(GraphType::edge_descriptor edge, BondStereocenter stereocenter);

  //! Remove all stereocenters
  void clear();

  //! Fetch a reference-option to an AtomStereocenter, if present
  boost::optional<AtomStereocenter&> option(const AtomIndexType index);

  //! Fetch a reference-option to an BondStereocenter, if present
  boost::optional<BondStereocenter&> option(const GraphType::edge_descriptor edge);

  /*! Communicates the removal of a vertex index to all stereocenters in the list
   *
   * Removing a vertex invalidates some vertex descriptors, which are used
   * liberally in the stereocenter classes. This function ensures that
   * vertex descriptors are valid throughout.
   */
  void propagateVertexRemoval(const AtomIndexType removedIndex);

  //! Removes the AtomStereocenter on a specified index
  void remove(const AtomIndexType index);

  //! Removes the BondStereocenter on a specified edge
  void remove(const GraphType::edge_descriptor edge);

  //! Removes the AtomStereocenter on a specified index, if present
  void try_remove(const AtomIndexType index);

  //! Removes the BondStereocenter on a specified edge, if present
  void try_remove(const GraphType::edge_descriptor edge);

/* Information */
  //! Modular comparison with another StereocenterList using a bitmask
  bool compare(
    const StereocenterList& other,
    const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
  ) const;

  //! Returns true if there are no stereocenters
  bool empty() const;

  //! Returns true if there are any stereocenters with zero possible assignments
  bool hasZeroAssignmentStereocenters() const;

  //! Returns true if there are unassigned stereocenters
  bool hasUnassignedStereocenters() const;

  //! Fetch a const ref-option to an AtomStereocenter, if present
  boost::optional<const AtomStereocenter&> option(const AtomIndexType index) const;

  //! Fetch a const ref-option to an BondStereocenter, if present
  boost::optional<const BondStereocenter&> option(const GraphType::edge_descriptor edge) const;

  //! Combined size of atom and bond-stereocenter lists
  unsigned size() const;

/* Iterators */
  //! Returns an iterable object with modifiable atom stereocenter references
  boost::range_detail::select_second_mutable_range<AtomMapType> atomStereocenters();

  //! Returns an iterable object with unmodifiable atom stereocenter references
  boost::range_detail::select_second_const_range<AtomMapType> atomStereocenters() const;

  //! Returns an iterable object with modifiable bond stereocenter references
  boost::range_detail::select_second_mutable_range<BondMapType> bondStereocenters();

  //! Returns an iterable object with unmodifiable bond stereocenter references
  boost::range_detail::select_second_const_range<BondMapType> bondStereocenters() const;

/* Operators */
  //! Strict equality comparison
  bool operator == (const StereocenterList& other) const;
  bool operator != (const StereocenterList& other) const;

private:
  //! The underlying storage for atom stereocenters
  AtomMapType _atomStereocenters;

  //! The underlying storage for bond stereocenters
  BondMapType _bondStereocenters;
};

}

#endif
