/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Molecule class interface
 *
 * Contains the Molecule class declaration, which is the central class of the
 * library.
 *
 * @todo
 * - Dynamism of Eta bond type is not implemented. Under molecule edits, bonds
 *   may become or cease to be eta bond types, and this is not checked for
 *   and appropriately altered in the graph.
 */

#ifndef INCLUDE_MOLASSEMBLER_MOLECULE_H
#define INCLUDE_MOLASSEMBLER_MOLECULE_H

#include "Utils/ElementTypes.h"
#include "boost/optional.hpp"
#include "chemical_symmetries/Names.h"
#include "temple/constexpr/Bitmask.h"

#include "molassembler/AngstromWrapper.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

#include <memory>

// External forward declarations
namespace Scine {
namespace Utils {
using ElementTypeCollection = std::vector<ElementType>;
} // namespace Utils
} // namespace Scine

/**
 * @brief Master namespace to which the molassembler project belongs
 */
namespace Scine {

/**
 * @brief Molecular graph and stereocenter manipulation
 */
namespace molassembler {

// Forward declarations
class OuterGraph;
class StereopermutatorList;
struct RankingInformation;

/*!
 * @brief Models a molecule as a graph (connectivity of atoms) and a list of
 *   stereopermutators.
 *
 * This class models chemical molecules as a combination of a graph
 * and a list of stereopermutators.
 *
 * In the graph, vertices additionally store an element type (and thus
 * represent atoms) and edges store a discretized bond type (and therefore
 * represent bonds).
 *
 * Stereopermutators represent the absolute configuration at an atom or bond in an
 * abstract fashion and do not store coordinate information.
 *
 * @note You may be surprised to see that all basic editing of Molecules, even
 * if it seems to concern the graph only, happens in this class interface. That
 * is because a graph edit may affect rankings at any stereopermutators in the
 * molecule due to the algorithm by which substituents are ranked. This means
 * that for every tiny edit, all stereopermutator substituents are re-ranked and
 * chiral state, if present, is propagated through a possible ranking change.
 * For that, the list of stereopermutators is required, which is accessible only
 * in this class.
 */
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
    Utils::ElementType a,
    Utils::ElementType b,
    BondType bondType
  ) noexcept;

  /*!
   * @brief Constructs from connectivity alone, inferring the stereopermutators from graph
   *
   * Constructs a molecule from connectivity alone. Local symmetries and
   * stereopermutators are inferred from the graph alone.
   *
   * @throws std::logic_error If the supplied graph has multiple connected
   *   components or there are less than 2 atoms
   */
  explicit Molecule(OuterGraph graph);

  /*!
   * @brief Construct from connectivity and positions
   *
   * Construct an instance from a constituting graph and positional information.
   * Local symmetries are deduced from positional information. Stereopermutators
   * are inferred from the graph and assigned using the supplied positional
   * information.
   *
   * @param bondStereopermutatorCandidatesOptional If boost::none, all bonds are
   *   candidates for BondStereopermutator. Otherwise, only the specified bonds are
   *   checked for BondStereopermutators.
   *
   * @throws std::logic_error If the supplied graph has multiple connected
   *   components or there are less than 2 atoms
   */
  Molecule(
    OuterGraph graph,
    const AngstromWrapper& positions,
    const boost::optional<
      std::vector<BondIndex>
    >& bondStereopermutatorCandidatesOptional = boost::none
  );

  /*!
   * @brief Construct a molecule from the underlying data fragments
   *
   * @throws std::logic_error If the supplied graph has multiple connected
   *   components or there are less than 2 atoms
   *
   * @warning This function is not intended for library consumers. It is used
   * internally in implementation details.
   */
  Molecule(
    OuterGraph graph,
    StereopermutatorList stereopermutators
  );
//!@}

//!@name Modifiers
//!@{
  /*!
   * @brief Adds an atom by attaching it to an existing atom.
   *
   * Adds a new atom, attaching it to an existing atom by a specified bond type.
   *
   * @param elementType The element type of the new atom
   * @param adjacentTo The atom to which the new atom is to be attached
   * @param bondType The bond type with which the new atom is to be attached
   *
   * @throws std::out_of_range If adjacentTo is invalid, i.e. >= N()
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  AtomIndex addAtom(
    Scine::Utils::ElementType elementType,
    AtomIndex adjacentTo,
    BondType bondType
  );

  /*!
   * @brief Adds a bond between two existing atoms
   *
   * Adds a bond between two already-existing atoms.
   *
   * @param a The first atom index
   * @param b The second atom index
   * @param bondType The bond type with which to connect a and b.
   *
   * @throws std::out_of_range If either atom index is invalid, i.e. >= N()
   * @throws std::logic_error If the atom indices match or the edge already
   *   exists.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void addBond(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  /*!
   * @brief Sets the stereopermutator assignment at a particular atom
   *
   * This sets the stereopermutator assignment at a specific atom index. For this,
   * a stereopermutator must be instantiated and contained in the StereopermutatorList
   * returned by stereopermutators(). The supplied assignment must be either
   * boost::none or smaller than stereopermutatorPtr->numAssignments().
   *
   * @param a The atom index at which a stereopermutator is to be assigned.
   * @param assignmentOption The new assignment. The special value boost::none
   *   makes the stereopermutator indeterminate. Any indeterminate atom
   *   stereopermutators in a molecule at conformation-generation will be assigned
   *   with a probability according to their sterepermutation's relative
   *   statistical occurence.
   *
   * @throws std::out_of_range If the atom index is invalid, i.e. >= N(),
   *   there is no stereopermutator at this position or the assignment index is
   *   invalid.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void assignStereopermutator(
    AtomIndex a,
    const boost::optional<unsigned>& assignmentOption
  );

  /*!
   * @brief Sets the stereopermutator assignment on a bond
   *
   * This sets the stereopermutator assignment at a specific bond index. For this,
   * a stereopermutator must be instantiated and contained in the StereopermutatorList
   * returned by stereopermutators(). The supplied assignment must be either
   * boost::none or smaller than stereopermutatorPtr->numAssignments().
   *
   * @param edge The edge at which a stereopermutator is to be assigned.
   * @param assignmentOption The new assignment. The special value boost::none
   *   makes the stereopermutator indeterminate. Any indeterminate bond
   *   stereopermutators in a molecule at conformation-generation will be assigned
   *   randomly.
   *
   * @throws std::out_of_range If the BondIndex is invalid (i.e. either atom
   *   index >= N()), there is no bond stereopermutator at the supplied edge
   *   or the assignment index is inalid.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void assignStereopermutator(
    const BondIndex& edge,
    const boost::optional<unsigned>& assignmentOption
  );

  /*!
   * @brief Assigns a stereopermutator stereopermutation at random
   *
   * This sets the stereopermutator assignment at a specific index, taking relative
   * statistical occurence weights of each stereopermutation into account.
   *
   * @param a The atom index at which the atom stereopermutator is to be assigned
   *   randomly
   *
   * @throws std::out_of_range If the atom index is invalid (i.e. is >= N()) or
   *   there is no atom stereopermutator at this bond index.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void assignStereopermutatorRandomly(AtomIndex a);

  /*!
   * @brief Assigns a bond stereopermutator to a random assignment
   *
   * @throws std::out_of_range If the bond index is invalid (i.e. either atom
   *   index is >= N()) or there is no bond stereopermutator at this bond index.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void assignStereopermutatorRandomly(const BondIndex& e);

  /*!
   * @brief Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * @throws std::out_of_range If the supplied index is invalid, i.e. >= N()
   * @throws std::logic_error If removing the atom disconnects the graph.
   *
   * @warning Invalidates **all** atom indices due to renumbering
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void removeAtom(AtomIndex a);

  /*!
   * @brief Removes a bond from the graph
   *
   * Removes an atom after checking if removing that bond is safe, i.e. does not
   * disconnect the graph. An example of bonds that can always be removed are
   * ring-closing bonds, since they never disconnect the molecular graph.
   *
   * @throws std::out_of_range If the supplied bond index is invalid, i.e. either
   *   atom index >= N() or the specified bond does not exist.
   * @throws std::logic_error If graph().canRemove() returns false.
   *
   * @note It is not safe to remove a bond just because one of the involved
   *   atoms is terminal, since that atom would then be disconnected from the
   *   rest of the molecule. This function merely removes a bond from the graph.
   *   It is, however, considered safe to remove the terminal vertex, which
   *   involves removing the bond to it.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void removeBond(AtomIndex a, AtomIndex b);

  /*!
   * @brief Changes a bond type. Returns whether the bond already existed
   *
   * Changes the bond type between two atom indices. If the bond does not exist
   * yet, adds the bond.
   *
   * @param a The first index of the bond whose type should be changed
   * @param b The second index of the bond whose type should be changed
   * @param bondType The new bond type
   *
   * @return If the bond type already existed.
   *
   * @throws out_of_range If a or b are invalid, i.e. >= N()
   * @throws std::logic_error If bondType is specified as BondType::Eta. The
   *   representation of bonding to haptic ligands and its dynamism is handled
   *   internally.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  bool setBondType(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  /*!
   * @brief Changes an existing atom's element type
   *
   * Changes the element type of an existing atom.
   *
   * @param a The atom index of the atom whose element type is to be changed
   * @param The new element type
   *
   * @throws std::out_of_range If a is invalid >= N()
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void setElementType(
    AtomIndex a,
    Utils::ElementType elementType
  );

  /*!
   * @brief Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a AtomStereopermutator instantiated at this atom index, its
   * underlying symmetry is altered. If there is no AtomStereopermutator at
   * this index, one is instantiated. In all cases, new or modified
   * stereopermutators are default-assigned if there is only one possible
   * assignment.
   *
   * @throws if
   *   - the supplied atomic index is invalid
   *   - the provided symmetry is a different size than that of an existing
   *     AtomStereopermutator
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void setGeometryAtAtom(
    AtomIndex a,
    Symmetry::Name symmetryName
  );
//!@}

//!@name Information
//!@{
  /*!
   * @brief Determines what the local geometry at a non-terminal atom ought to be
   *
   * Returns the expected symmetry name at a non-terminal atom by inference
   * from graph information only.
   *
   * Can be extended to apply various levels of theory to determine the
   * geometry based on various information contained in the graph, such as
   * ligand field theory.
   *
   * @param index The atom index where the local geometry should be determined
   * @param ranking The ranking of all substituents at the supplied index. See
   *   rankPriority()
   *
   * @returns A local geometry (of the appropriate size to fit the number of
   *   ligand sites) which may or may not be the symmetry with the lowest
   *   energy.
   *
   * @throws std::out_of_range If the supplied atomic index is >= N()
   * @throws std::logic_error If the supplied atom is terminal or if the amount
   *   of ligand sites exceeds the size of the largest defined symmetry.
   *
   * @note Currently applies VSEPR if the element type at the supplied index is
   *   not a transition metal, and returns the first symmetry of appropriate
   *   size otherwise.
   */
  Symmetry::Name determineLocalGeometry(
    AtomIndex index,
    const RankingInformation& ranking
  ) const;

  /*!
   * @brief Returns a graphivz string representation of the molecule
   *
   * Creates a graphviz representation of a molecule that can be written into a
   * dotfile and processed with graphviz's `dot` binary to create an image of
   * the molecular graph.
   *
   * Includes tooltip information on stereopermutators when hovering over
   * individual nodes.
   */
  std::string dumpGraphviz() const;

  //! Provides read-only access to the graph representation
  const OuterGraph& graph() const;

  //! Provides read-only access to the list of stereopermutators
  const StereopermutatorList& stereopermutators() const;

  /*!
   * @brief Generates stereopermutators from connectivity and positional information
   *
   * Positions are an important source of information for stereopermutators as they
   * will alleviate graph-based symmetry-determination errors and allow for the
   * determination of stereopermutator assignments through spatial fitting.
   *
   * @param angstromWrapper Wrapped positions in angstrom length units
   * @param explicitBondStereopermutatorCandidatesOption Permits the specification
   *   of a limited set of bonds on which BondStereopermutator instantiation is
   *   attempted. In Interpret.h, for instance, you can choose not to
   *   instantiate BondStereopermutators below a fractional bond order threshold to
   *   avoid spurious frozen dihedrals. By default, all bonds are candidates.
   *
   * @throws std::out_of_range if a BondIndex in
   *   explicitBondStereopermutatorCandidatesOption does not reference an existing
   *   bond (irrelevant if left default).
   */
  StereopermutatorList inferStereopermutatorsFromPositions(
    const AngstromWrapper& angstromWrapper,
    const boost::optional<
      std::vector<BondIndex>
    >& explicitBondStereopermutatorCandidatesOption = boost::none
  ) const;

  /*!
   * @brief Modular comparison of this Molecule with another.
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
   * stereopermutators across both molecules are compared using the found
   * isomorphism as an index map.
   *
   * @note The number of stereopermutations that a stereopermutator has is
   * considered part of the Symmetry ComparisonOptions.
   *
   * @note If you choose to discard bond order checking, this merely
   * deactivates bond order hashing and a post-isomorphism-search bond order
   * re-check. Bond order information - if present in the molecule prior to
   * calling this function - is also present in stereopermutator ranking information
   * and hence can influence the number of stereopermutations and the currently
   * set stereopermutation index. This can lead to unexpected but logically
   * consistent comparison behavior.
   */
  bool modularCompare(
    const Molecule& other,
    const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
  ) const;

  /*!
   * @brief Rank substituents of an atom
   *
   * Performs a ranking algorithm that attempts to differentiate branches
   * extending at each substituent atom (haptic ligands are not considered a
   * single unit, rather their component atoms are all individual substituents).
   *
   * Groups substituents into ligands (these are not the same since haptic
   * ligands exist) and ranks those too.
   *
   * @param a The atom whose substituents are to be ranked
   * @param excludeAdjacent A list of substituent atom indices that should be
   *   excluded from ranking
   * @param positionsOption Positional information can be used to determine
   *   auxiliary stereopermutator assignments that arise in the ranking algorithm
   *   and may have different sub-rankings than at the same position in the
   *   molecule considered on its own. It is preferable to supply this if
   *   positional information is present!
   *
   * @throws std::out_of_range If the supplied atom index is i.e. >= N()
   *
   * @returns a RankingInformation instance that contains all gathered
   * information.
   */
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
  //! Inverts @see operator ==
  bool operator != (const Molecule& other) const;
//!@}

private:
  //! Private implementation member
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

} // namespace Scine

std::ostream& operator << (
  std::ostream& os,
  const Scine::molassembler::Molecule& molecule
);

#endif
