/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Private implementation of Molecule
 */

#ifndef INCLUDE_MOLASSEMBLER_MOLECULE_IMPL_H
#define INCLUDE_MOLASSEMBLER_MOLECULE_IMPL_H

#include "molassembler/Molecule.h"

#include "molassembler/Graph.h"
#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/StereopermutatorList.h"
#include "Utils/Geometry/AtomCollection.h"

namespace Scine {
namespace Molassembler {

struct Molecule::Impl {
  static Utils::AtomCollection applyCanonicalizationMap(
    const std::vector<AtomIndex>& canonicalizationIndexMap,
    const Utils::AtomCollection& atomCollection
  );

  Graph adjacencies_;
  StereopermutatorList stereopermutators_;
  boost::optional<AtomEnvironmentComponents> canonicalComponentsOption_;

/* "Private" helpers */
  void tryAddAtomStereopermutator_(
    AtomIndex candidateIndex,
    StereopermutatorList& stereopermutators
  ) const;

  void tryAddBondStereopermutator_(
    const BondIndex& bond,
    StereopermutatorList& stereopermutators
  ) const;

  //! Generates a list of stereopermutators based on graph properties alone
  StereopermutatorList detectStereopermutators_() const;

  //! Ensures basic expectations about what constitutes a Molecule are met
  void ensureModelInvariants_() const;

  //! Returns whether the specified index is valid or not
  bool isValidIndex_(AtomIndex index) const;

  //! Returns whether an edge is double, triple or higher bond order
  bool isGraphBasedBondStereopermutatorCandidate_(BondType bondType) const;

  //! Updates the molecule's StereopermutatorList after a graph modification
  void propagateGraphChange_();


//!@name Constructors
//!@{
  //! Default constructor
  Impl() noexcept;

  //! Mono-atomic constructor
  Impl(Utils::ElementType element) noexcept;

  //! Diatomic constructor
  Impl(
    Utils::ElementType a,
    Utils::ElementType b,
    BondType bondType
  ) noexcept;

  //! Graph-only constructor
  explicit Impl(Graph graph);

  //! Graph and positions constructor
  Impl(
    Graph graph,
    const AngstromPositions& positions,
    const boost::optional<
      std::vector<BondIndex>
    >& bondStereopermutatorCandidatesOptional = boost::none
  );

  //! Graph and stereopermutators constructor
  Impl(
    Graph graph,
    StereopermutatorList stereopermutators,
    boost::optional<AtomEnvironmentComponents> canonicalComponentsOption
  );
//!@}

//!@name Modifiers
//!@{
  //! Adds an atom by attaching it to an existing atom.
  AtomIndex addAtom(
    Utils::ElementType elementType,
    AtomIndex adjacentTo,
    BondType bondType
  );

  //! Adds a bond between existing atoms.
  BondIndex addBond(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  //! Applies an index permutation to all member state
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*! Sets the stereopermutator assignment at a particular atom
   *
   * This sets the stereopermutator assignment at a specific atom index. For this,
   * a stereopermutator must be instantiated and contained in the StereopermutatorList
   * returned by stereopermutators(). The supplied assignment must be either
   * boost::none or smaller than stereopermutatorPtr->numAssignments().
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereopermutatorList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereopermutators or the removal of old ones.
   */
  void assignStereopermutator(
    AtomIndex a,
    const boost::optional<unsigned>& assignment
  );

  /*! Sets the stereopermutator assignment on a bond
   *
   * This sets the stereopermutator assignment at a specific bond index. For this,
   * a stereopermutator must be instantiated and contained in the StereopermutatorList
   * returned by stereopermutators(). The supplied assignment must be either
   * boost::none or smaller than stereopermutatorPtr->numAssignments().
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereopermutatorList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereopermutators or the removal of old ones.
   */
  void assignStereopermutator(
    const BondIndex& edge,
    const boost::optional<unsigned>& assignment
  );

  /*! Assigns a stereopermutator stereopermutation at random
   *
   * This sets the stereocetner assignment at a specific index, taking relative
   * statistical occurence weights of each stereopermutation into account.
   *
   * \pre There must be an AtomStereopermutator at the passed index
   *
   * \throws If no AtomStereopermutator exists at the passed index
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereopermutatorList you have accessed prior to calling this function and
   * particularly any iterators to its members invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereopermutators or the removal of old ones.
   */
  void assignStereopermutatorRandomly(AtomIndex a, Random::Engine& engine);

  /*! Assigns a bond stereopermutator to a random assignment
   *
   * \pre There must be a BondStereopermutator at the passed edge
   * \throws If no BondStereopermutator exists at the passed edge
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereopermutatorList you have accessed prior to calling this function and
   * particularly any iterators to its members invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereopermutators or the removal of old ones.
   */
  void assignStereopermutatorRandomly(const BondIndex& e, Random::Engine& engine);

  /**
   * @brief Canonicalizes the graph, invalidating all atom and bond indices.
   *
   * @param componentBitmask The components of the molecular graph to include in the
   *   canonicalization procedure.
   *
   * @warning Any comparisons made on canonical graphs must be made with a less
   *   or equally strict @p components bitmask. If you choose a stricter
   *   bitmask than you have used in the canonicalization procedure, you risk
   *   false positives and false negatives.
   *
   * @warning This invalidates all atom indices and bond indices and any
   *   references to constituting members of the molecule.
   *
   * @note Use @see canonicalCompare to compare instances of canonicalized
   *   molecules.
   *
   * @return Permutation mapping from old indices to new:
   * @code{.cpp}
   * auto indexMapping = mol.canonicalize();
   * AtomIndex newIndex = indexMapping.at(oldIndex);
   * @endcode
   * You can use this to update invalidated indices.
   */
  std::vector<AtomIndex> canonicalize(
    AtomEnvironmentComponents componentBitmask
  );

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

  //! Changes an existing bond's type
  bool setBondType(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  //! Changes an existing atom's element type
  void setElementType(
    AtomIndex a,
    Utils::ElementType elementType
  );

  /*! Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a AtomStereopermutator instantiated at this atom index, its
   * underlying shape is altered. If there is no AtomStereopermutator at
   * this index, one is instantiated. In all cases, new or modified
   * stereopermutators are default-assigned if there is only one possible
   * assignment.
   * \throws if
   *   - the supplied atomic index is invalid
   *   - there is an BondStereopermutator at that index
   *   - or the provided shape is a different size than that of an existing
   *     AtomStereopermutator or the expected shape
   */
  void setShapeAtAtom(
    AtomIndex a,
    Shapes::Shape shape
  );
//!@}

//!@name Information
//!@{
  //! Yield which components were used in canonicalization
  boost::optional<AtomEnvironmentComponents> canonicalComponents() const;

  /*! Determines what the local geometry at a non-terminal atom ought to be
   *
   * Returns the expected shape name at a non-terminal atom.
   * \throws if the supplied atomic index is invalid
   */
  boost::optional<Shapes::Shape> inferShape(
    AtomIndex index,
    const RankingInformation& ranking
  ) const;

  //! Returns a graphivz string representation of the molecule
  std::string dumpGraphviz() const;

  //! Provides read-only access to the graph member
  const Graph& graph() const;

  //! Convolutional hash
  std::size_t hash() const;

  //! Provides read-only access to the list of stereopermutators
  const StereopermutatorList& stereopermutators() const;

  StereopermutatorList inferStereopermutatorsFromPositions(
    const AngstromPositions& angstromWrapper,
    const boost::optional<
      std::vector<BondIndex>
    >& explicitBondStereopermutatorCandidatesOption = boost::none
  ) const;

  //! Compares two canonical instances with one another
  bool canonicalCompare(
    const Impl& other,
    AtomEnvironmentComponents componentBitmask
  ) const;

  //! Modular comparison of this Impl with another.
  boost::optional<std::vector<AtomIndex>> modularIsomorphism(
    const Impl& other,
    AtomEnvironmentComponents componentBitmask
  ) const;

  //! Returns a command-line interface information string
  std::string str() const;

  RankingInformation rankPriority(
    AtomIndex a,
    const std::vector<AtomIndex>& excludeAdjacent = {},
    const boost::optional<AngstromPositions>& positionsOption = boost::none
  ) const;
//!@}

//!@name Operators
//!@{
  //! Equality operator, performs most strict equality comparison
  bool operator == (const Impl& other) const;
  //! Negates @see operator ==
  bool operator != (const Impl& other) const;
//!@}
};

} // namespace Molassembler
} // namespace Scine

#endif
