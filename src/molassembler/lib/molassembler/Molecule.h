#ifndef INCLUDE_MOLASSEMBLER_MOLECULE_H
#define INCLUDE_MOLASSEMBLER_MOLECULE_H

#include "Delib/ElementTypes.h"

#include "molassembler/detail/RangeForTemporary.h"
#include "molassembler/detail/AngstromWrapper.h"
#include "molassembler/AtomEnvironmentHash.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

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
class Cycles;
class RankingInformation;

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

  //! Construct a minimal molecule from two element types and a shared bond type
  Molecule(
    Delib::ElementType a,
    Delib::ElementType b,
    BondType bondType
  ) noexcept;

  //! Constructs a molecule from connectivity alone, inferring the stereocenters
  explicit Molecule(GraphType graph);

  /*! Construct a molecule from connectivity and 3D information.
   *
   * \note Assumes that the provided position collection is in Angstrom units.
   */
  Molecule(
    GraphType graph,
    const AngstromWrapper& positions
  );

  //! Construct a molecule from the underlying data fragments
  Molecule(
    GraphType graph,
    StereocenterList stereocenters
  );
//!@}

//!@name Modifiers
//!@{
  //! Adds an atom by attaching it to an existing atom.
  AtomIndexType addAtom(
    Delib::ElementType elementType,
    AtomIndexType adjacentTo,
    BondType bondType
  );

  //! Adds a bond between existing atoms.
  void addBond(
    AtomIndexType a,
    AtomIndexType b,
    BondType bondType
  );

  /*! Sets the stereocenter assignment at a particular atom
   *
   * This sets the stereocenter assignment at a specific atom index. For this,
   * a stereocenter must be instantiated and contained in the StereocenterList
   * returned by getStereocenterList(). The supplied assignment must be either
   * boost::none or smaller than stereocenterPtr->numAssignments().
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenter(
    AtomIndexType a,
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
  void assignStereocenterRandomly(AtomIndexType a);

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
  void assignStereocenterRandomly(const GraphType::edge_descriptor& e);

  /*! Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * \throws if the supplied index is invalid or isSafeToRemoveAtom returns false.
   */
  void removeAtom(AtomIndexType a);

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
  void removeBond(AtomIndexType a, AtomIndexType b);

  //! Changes an existing bond's type
  bool setBondType(
    AtomIndexType a,
    AtomIndexType b,
    BondType bondType
  );

  //! Changes an existing atom's element type
  void setElementType(
    AtomIndexType a,
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
    AtomIndexType a,
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
    AtomIndexType index,
    const RankingInformation& ranking
  ) const;

  //! Returns a graphivz string representation of the molecule
  std::string dumpGraphviz() const;

  //! Get an edge descriptor for a pair of indices
  GraphType::edge_descriptor edge(
    AtomIndexType a,
    AtomIndexType b
  ) const;

  /*! Fetches the atomic indices of vertices adjacent to a particular index
   *
   * Fetches the atomic indices of vertices adjacent to a particular index.
   * \throws if the supplied atomic index is invalid
   */
  std::vector<AtomIndexType> getAdjacencies(AtomIndexType a) const;

  //! Fetches the optional bond type between two atom indices
  boost::optional<BondType> getBondType(
    AtomIndexType a,
    AtomIndexType b
  ) const;

  BondType getBondType(const GraphType::edge_descriptor& edge) const;

  Cycles getCycleData() const;

  /*! Returns the element type of an atomic index
   *
   * Returns the element type of an atomic index
   * \throws if the atomic index is invalid
   */
  Delib::ElementType getElementType(AtomIndexType index) const;

  //! Returns a collection detailing all element types
  Delib::ElementTypeCollection getElementCollection() const;

  //! Provides read-only access to the graph member
  const GraphType& getGraph() const;

  //! Provides read-only access to the list of stereocenters
  const StereocenterList& getStereocenterList() const;

  //! Returns the number of adjacencies of an atomic position
  unsigned getNumAdjacencies(AtomIndexType a) const;

  StereocenterList inferStereocentersFromPositions(const AngstromWrapper& angstromWrapper) const;

  //! Checks if two atomic indices are connected by a bond
  bool isAdjacent(
    AtomIndexType a,
    AtomIndexType b
  ) const;

  //! An atom is considered removable if it isn't an articulation vertex
  bool isSafeToRemoveAtom(AtomIndexType a) const;

  //! A bond is considered removable if it isn't a bridge edge
  bool isSafeToRemoveBond(AtomIndexType a, AtomIndexType b) const;

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

  //! Fetch the number of atoms
  unsigned numAtoms() const;

  //! Fetch the number of bonds
  unsigned numBonds() const;

  RankingInformation rankPriority(
    AtomIndexType a,
    const std::set<AtomIndexType>& excludeAdjacent = {},
    const boost::optional<AngstromWrapper>& positionsOption = boost::none
  ) const;

  //! Get the vertex indices on both ends of a graph edge
  std::array<AtomIndexType, 2> vertices(const GraphType::edge_descriptor& edge) const;
//!@}

//!@name Iterators
//!@{
  //! Returns a range-for temporary object iterating through all atom indices
  RangeForTemporary<GraphType::vertex_iterator> iterateAtoms() const;

  /*! Returns a range-for temporary object allowing c++11 style for loop
   * iteration through an atom's adjacencies
   */
  RangeForTemporary<GraphType::adjacency_iterator> iterateAdjacencies(
    AtomIndexType a
  ) const;

  /*! Returns a range-for temporary object allowing c++11-style for loop
   * iteration through edges
   */
  RangeForTemporary<GraphType::edge_iterator> iterateEdges() const;

  /*! Returns a range-for temporary object allowing c++11-style for loop
   * iteration through edges around a specific atom
   */
  RangeForTemporary<GraphType::out_edge_iterator> iterateEdges(AtomIndexType a) const;
//!@}


//!@name Operators
//!@{
  //! Returns the adjacencies of the specified atom index
  RangeForTemporary<GraphType::adjacency_iterator> operator [] (
    AtomIndexType a
  ) const;

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
