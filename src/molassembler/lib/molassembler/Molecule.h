#ifndef INCLUDE_MOLASSEMBLER_MOLECULE_H
#define INCLUDE_MOLASSEMBLER_MOLECULE_H

#include "Delib/ElementTypes.h"
#include "boost/optional.hpp"
#include "chemical_symmetries/Names.h"
#include "temple/constexpr/Bitmask.h"

#include "molassembler/Containers/AngstromWrapper.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

#include <memory>

/* TODO
 * - Dynamism of Eta bond type is not implemented. Under molecule edits, bonds
 *   may become or cease to be eta bond types.
 */

/*! @file
 *
 * Contains the Molecule class declaration, which is the central class of the
 * library.
 */

// External forward declarations
namespace Delib {
class ElementTypeCollection;
} // namespace Delib

namespace molassembler {

// Forward declarations
class OuterGraph;
class StereocenterList;
struct RankingInformation;

//! Central class of the library, modeling a molecular graph with all state.
class Molecule {
public:
//!@name Special member functions
//!@{
  /* Rule of five members */
  Molecule(Molecule&& other) noexcept;
  Molecule& operator = (Molecule&& rhs) noexcept;
  Molecule(const Molecule& other);
  Molecule& operator = (const Molecule& rhs);
  ~Molecule();

  //! Default-constructor creates a hydrogen molecule.
  Molecule() noexcept;

  //! Construct a minimal molecule from two element types and a mutual bond type
  Molecule(
    Delib::ElementType a,
    Delib::ElementType b,
    BondType bondType
  ) noexcept;

  //! Constructs from connectivity alone, inferring the stereocenters from graph
  explicit Molecule(OuterGraph graph);

  /*! Construct from connectivity and positions
   *
   * \param bondStereocenterCandidatesOptional If boost::none, all bonds are
   *   candidates for BondStereocenter. Otherwise, only the specified bonds are
   *   checked for BondStereocenters.
   */
  Molecule(
    OuterGraph graph,
    const AngstromWrapper& positions,
    const boost::optional<
      std::vector<BondIndex>
    >& bondStereocenterCandidatesOptional = boost::none
  );

  //! Construct a molecule from the underlying data fragments
  Molecule(
    OuterGraph graph,
    StereocenterList stereocenters
  );
//!@}

//!@name Modifiers
//!@{
  //! Adds an atom by attaching it to an existing atom.
  AtomIndex addAtom(
    Delib::ElementType elementType,
    AtomIndex adjacentTo,
    BondType bondType
  );

  //! Adds a bond between existing atoms.
  void addBond(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  /*! Sets the stereocenter assignment at a particular atom
   *
   * This sets the stereocenter assignment at a specific atom index. For this,
   * a stereocenter must be instantiated and contained in the StereocenterList
   * returned by stereocenters(). The supplied assignment must be either
   * boost::none or smaller than stereocenterPtr->numAssignments().
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenter(
    AtomIndex a,
    const boost::optional<unsigned>& assignment
  );

  /*! Assigns a stereocenter stereopermutation at random
   *
   * This sets the stereocetner assignment at a specific index, taking relative
   * statistical occurence weights of each stereopermutation into account.
   *
   * \pre There must be an AtomStereocenter at the passed index
   *
   * \throws If no AtomStereocenter exists at the passed index
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators to its members invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenterRandomly(AtomIndex a);

  /*! Assigns a bond stereocenter to a random assignment
   *
   * \pre There must be a BondStereocenter at the passed edge
   * \throws If no BondStereocenter exists at the passed edge
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators to its members invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenterRandomly(const BondIndex& e);

  /*! Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * \throws if the supplied index is invalid or isSafeToRemoveAtom returns false.
   */
  void removeAtom(AtomIndex a);

  /*!
   * Removes an atom after checking if removing that bond is safe, i.e. does not
   * disconnect the graph. An example of bonds that can always be removed are
   * ring-closing bonds, since they never disconnect the molecular graph.
   *
   * \throws if isSafeToRemoveBond returns false.
   *
   * \note It is not safe to remove a bond just because one of the involved
   * atoms is terminal, since that atom would then be disconnected from the
   * rest of the molecule. This function merely removes a bond from the graph.
   * It is, however, considered safe to remove the terminal vertex, which
   * involves removing the bond to it.
   */
  void removeBond(AtomIndex a, AtomIndex b);

  //! Changes a bond type. Returns whether the bond already existed
  bool setBondType(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  //! Changes an existing atom's element type
  void setElementType(
    AtomIndex a,
    Delib::ElementType elementType
  );

  /*! Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a AtomStereocenter instantiated at this atom index, its
   * underlying symmetry is altered. If there is no AtomStereocenter at
   * this index, one is instantiated. In all cases, new or modified
   * stereocenters are default-assigned if there is only one possible
   * assignment.
   *
   * \throws if
   *   - the supplied atomic index is invalid
   *   - there is an BondStereocenter at that index
   *   - or the provided symmetry is a different size than that of an existing
   *     AtomStereocenter or the expected symmetry
   */
  void setGeometryAtAtom(
    AtomIndex a,
    Symmetry::Name symmetryName
  );
//!@}

//!@name Information
//!@{
  /*! Determines what the local geometry at a non-terminal atom ought to be
   *
   * Returns the expected symmetry name at a non-terminal atom.
   * \throws if the supplied atomic index is invalid
   */
  Symmetry::Name determineLocalGeometry(
    AtomIndex index,
    const RankingInformation& ranking
  ) const;

  //! Returns a graphivz string representation of the molecule
  std::string dumpGraphviz() const;

  //! Provides read-only access to the graph member
  const OuterGraph& graph() const;

  //! Provides read-only access to the list of stereocenters
  const StereocenterList& stereocenters() const;

  /*! Generates stereocenters from connectivity and positional information
   *
   * Positions are an important source of information for stereocenters as they
   * will alleviate graph-based symmetry-determination errors and allow for the
   * determination of stereocenter assignments through spatial fitting.
   *
   * \param angstromWrapper Wrapped positions in angstrom length units
   * \param explicitBondStereocenterCandidatesOption Permits the specification
   *   of a limited set of bonds on which BondStereocenter instantiation is
   *   attempted. In Interpret.h, for instance, you can choose not to
   *   instantiate BondStereocenters below a fractional bond order threshold to
   *   avoid spurious frozen dihedrals. By default, all bonds are candidates.
   *
   * \throws std::domain_error if a BondIndex in
   *   explicitBondStereocenterCandidatesOption does not reference an existing
   *   bond (irrelevant if left default).
   */
  StereocenterList inferStereocentersFromPositions(
    const AngstromWrapper& angstromWrapper,
    const boost::optional<
      std::vector<BondIndex>
    >& explicitBondStereocenterCandidatesOption = boost::none
  ) const;

  /*! Modular comparison of this Molecule with another.
   *
   * This permits detailed specification of which elements of the molecular
   * information you want to use in the comparison.
   *
   * Equality comparison is performed in several stages: First, at each atom
   * position, a hash is computed that encompasses all local information that
   * is specified to be used in the comparisonBitmask. This hash is then used
   * during graph isomorphism calculation to avoid finding an isomorphism that
   * does not consider the specified factors.
   *
   * If an isomorphism is found, it is then validated. Bond orders and
   * stereocenters across both molecules are compared using the found
   * isomorphism as an index map.
   *
   * \note The number of stereopermutations that a stereocenter has is
   * considered part of the Symmetry ComparisonOptions.
   *
   * \note If you choose to discard bond order checking, this merely
   * deactivates bond order hashing and a post-isomorphism-search bond order
   * re-check. Bond order information - if present in the molecule prior to
   * calling this function - is also present in stereocenter ranking information
   * and hence can influence the number of stereopermutations and the currently
   * set stereopermutation index. This can lead to unexpected but logically
   * consistent comparison behavior.
   */
  bool modularCompare(
    const Molecule& other,
    const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
  ) const;

  RankingInformation rankPriority(
    AtomIndex a,
    const std::vector<AtomIndex>& excludeAdjacent = {},
    const boost::optional<AngstromWrapper>& positionsOption = boost::none
  ) const;
//!@}

//!@name Operators
//!@{
  //! Equality operator, performs most strict equality comparison
  bool operator == (const Molecule& other) const;
  bool operator != (const Molecule& other) const;
//!@}

private:
  struct Impl;

#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
  std::experimental::propagate_const<
    std::unique_ptr<Impl>
  > _pImpl;
#else
  std::unique_ptr<Impl> _pImpl;
#endif
};

} // namespace molassembler

std::ostream& operator << (
  std::ostream& os,
  const molassembler::Molecule& molecule
);

#endif
