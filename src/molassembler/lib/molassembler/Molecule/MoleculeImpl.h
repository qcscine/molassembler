#ifndef INCLUDE_MOLASSEMBLER_MOLECULE_IMPL_H
#define INCLUDE_MOLASSEMBLER_MOLECULE_IMPL_H

#include "molassembler/Molecule.h"

#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/StereopermutatorList.h"

namespace molassembler {

struct Molecule::Impl {
  OuterGraph _adjacencies;
  StereopermutatorList _stereopermutators;

/* "Private" helpers */
  //! Generates a list of stereopermutators based on graph properties alone
  StereopermutatorList _detectStereopermutators() const;

  //! Ensures basic expectations about what constitutes a Molecule are met
  void _ensureModelInvariants() const;

  //! Returns whether the specified index is valid or not
  bool _isValidIndex(AtomIndex index) const;

  //! Returns whether an edge is double, aromtic, triple or higher bond order
  bool _isGraphBasedBondStereopermutatorCandidate(const InnerGraph::Edge& e) const;

  //! Updates the molecule's StereopermutatorList after a graph modification
  void _propagateGraphChange();


//!@name Constructors
//!@{
  //! Default constructor
  Impl() noexcept;

  //! Diatomic constructor
  Impl(
    Delib::ElementType a,
    Delib::ElementType b,
    BondType bondType
  ) noexcept;

  //! Graph-only constructor
  explicit Impl(OuterGraph graph);

  //! Graph and positions constructor
  Impl(
    OuterGraph graph,
    const AngstromWrapper& positions,
    const boost::optional<
      std::vector<BondIndex>
    >& bondStereopermutatorCandidatesOptional = boost::none
  );

  //! Graph and stereopermutators constructor
  Impl(
    OuterGraph graph,
    StereopermutatorList stereopermutators
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
  void assignStereopermutatorRandomly(AtomIndex a);

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
  void assignStereopermutatorRandomly(const BondIndex& e);

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
    Delib::ElementType elementType
  );

  /*! Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a AtomStereopermutator instantiated at this atom index, its
   * underlying symmetry is altered. If there is no AtomStereopermutator at
   * this index, one is instantiated. In all cases, new or modified
   * stereopermutators are default-assigned if there is only one possible
   * assignment.
   * \throws if
   *   - the supplied atomic index is invalid
   *   - there is an BondStereopermutator at that index
   *   - or the provided symmetry is a different size than that of an existing
   *     AtomStereopermutator or the expected symmetry
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

  //! Provides read-only access to the list of stereopermutators
  const StereopermutatorList& stereopermutators() const;

  StereopermutatorList inferStereopermutatorsFromPositions(
    const AngstromWrapper& angstromWrapper,
    const boost::optional<
      std::vector<BondIndex>
    >& explicitBondStereopermutatorCandidatesOption = boost::none
  ) const;

  //! Modular comparison of this Impl with another.
  bool modularCompare(
    const Impl& other,
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
  bool operator == (const Impl& other) const;
  bool operator != (const Impl& other) const;
//!@}
};

} // namespace molassembler

#endif
