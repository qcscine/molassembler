#ifndef INCLUDE_MOLECULE_MANIP_MOLECULE_H
#define INCLUDE_MOLECULE_MANIP_MOLECULE_H

#include "Edges.h"
#include "StereocenterList.h"
#include "CycleData.h"
#include "VSEPR.h"
// #include "MemFnCache.h"

/*! @file
 *
 * Contains the Molecule class declaration, which is the central class of the
 * library.
 */

/* TODO
 * - Modifiers must update StereocenterList, invalidate Cache
 * - the results of detectStereocenters and inferStereocenters are both
 *   dependent on the sequence of candidates!  If e.g.::
 *   
 *          E   F
 *           % :
 *       B    2    B
 *         \ / \ /
 *      A ▶ 1   3 ◀ A
 *          △   △
 *          C   C
 *
 *    And if the sequence is 1-2-3, we get ABCD, ABCD, ABCD, which is correct.
 *    Each of 1-3 is a stereocenter.  If the sequence is 2-1-3, we get ABCC,
 *    ABCD, ABCD, which is incorrect. Needs to consider that rankings can
 *    differ since stereocenters may be as yet unassigned. One way would be to
 *    make a graph of dependencies and then resolve that in a fashion that they
 *    can be sequentially assigned. That may not be possible however!
 *
 *    What happens if::
 *
 *         A   B
 *          % /
 *           1
 *      A.  / \  .A
 *        :2 - 3:
 *      B°       °B
 *
 *   Dependencies are circular, yet there are clearly two ways of arranging the
 *   substituents relative to one another: all on one side or one on one side,
 *   two on the other. Can this even be achieved with permuting stereocenters?
 *
 * - How should changing a stereocenter's assignment work now? Via the public 
 *   member is no longer an option now that we have the cache, since there is no
 *   way to guarantee that the cache is invalidated and the list of stereocenters
 *   updated
 *
 * - Cannot leave it to the user to notify Molecule of changes in a 
 *   stereocenter, it ought to update all other stereocenters itself if that
 *   happens. Need a proxy object that notifies the StereocenterList that
 *   something was changed and updates all other stereocenters. for-range
 *   iteration through the map then has to be discouraged unleass I can modify
 *   Stereocenters in-place to accomodate ranking changes and am not forced to
 *   remove and re-insert, which will invalidate iterators to the affected
 *   stereocenter
 */


namespace MoleculeManip {

/*! 
 * Central class of the library, modeling a molecular graph with all state.
 */
class Molecule {
public:
  /*enum class CacheKeys {
    RemovalSafetyData
  };*/

private:
/* State */
  GraphType _adjacencies;
  StereocenterList _stereocenters;

  // Properties cache
  /*using CacheType = MemFnCache<CacheKeys, Molecule>;
  mutable CacheType _cache;*/

/* Members */
/* Private members */
  /*! 
   * Adds an vertex to the graph and sets it's element type property, returning
   * the new index.
   */
  AtomIndexType _addAtom(const Delib::ElementType& elementType);

  //! Returns the list of cache generators needed to initialize the cache
  /*std::initializer_list<
    std::pair<CacheKeys, CacheType::LambdaType>
  > _cacheGenerators() const;*/

  //! Returns whether the specified index is valid or not
  bool _isValidIndex(const AtomIndexType& index) const;

  //! Returns a list of atom indices that have at least 3 bonded neighbors
  std::vector<AtomIndexType> _getCNStereocenterCandidates() const;

  /*! 
   * Returns a list of edge indices where each endpoint has 1 or two additional
   * substituent besides the edge neighbor
   */
  std::vector<EdgeIndexType> _getEZStereocenterCandidates() const;

  //!  Reduces an atom's neighbors to ligand types
  std::vector<LocalGeometry::LigandType> _reduceToLigandTypes(
    const AtomIndexType& index
  ) const;

  //!  Updates the molecule's StereocenterList after a graph modification
  void _updateStereocenterList();

public:
/* Typedefs */
  using ExplicitEdge = std::pair<
    Edges::MapType::key_type, // pair<AtomIndexType, AtomIndexType>
    Edges::MapType::mapped_type // BondType
  >;

/* Constructors */
  /*! 
   * Default-constructor is deleted since the minimal molecule we represent is
   * two bonded atoms.
   */
  Molecule() = delete;

  //! Construct a minimal molecule from two element types and a shared bond type
  Molecule(
    const Delib::ElementType& a,
    const Delib::ElementType& b,
    const BondType& bondType
  );

  /*! 
   * Shorthand construction of a molecule using a list of element types and
   * a set of explicit edges
   */
  Molecule(
    const Delib::ElementTypeCollection& elements,
    const Edges& edges
  );

  //! Constructs a molecule from connectivity alone, inferring the stereocenters
  Molecule(const GraphType& graph);

  //! Construct a molecule from connectivity and 3D information
  Molecule(
    const GraphType& graph,
    const Delib::PositionCollection& positions
  );

/* Modification */
  //! Adds an atom by attaching it to an existing atom.
  AtomIndexType addAtom(
    const Delib::ElementType& elementType,
    const AtomIndexType& adjacentTo,
    const BondType& bondType
  );

  //! Adds a bond between existing atoms.
  void addBond(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bondType
  );

  //! Changes an existing atom's element type
  void changeElementType(
    const AtomIndexType& a,
    const Delib::ElementType& elementType
  );

  /*! 
   * Removes an atom after checking if removing that atom is safe, i.e. does
   * not disconnect the graph.
   *
   * Throws if isSafeToRemoveAtom returns false.
   */
  void removeAtom(const AtomIndexType& a);

  /*!
   * Removes an atom after checking if removing that bond is safe, i.e. does not
   * disconnect the graph. An example of bonds that can always be removed are
   * ring-closing bonds, since they never disconnect the molecular graph.
   *
   * Throws if isSafeToRemoveBond returns false. Note that it is not safe to
   * remove a bond just because one of the involved atoms is terminal, since
   * that atom would then be disconnected from the rest of the molecule. It is
   * however considered safe to remove the terminal vertex, which involves
   * removing the bond to it.
   */
  void removeBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  );

/* Information */
  StereocenterList detectStereocenters() const;

  Symmetry::Name determineLocalGeometry(
    const AtomIndexType& index
  ) const;

  std::string dumpGraphviz() const;

  std::vector<AtomIndexType> getAdjacencies(const AtomIndexType& a) const;

  boost::optional<BondType> getBondType(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

  CycleData getCycleData() const;

  // Creates a copy of the contained data suitable for the Edges class
  std::vector<ExplicitEdge> getEdges() const;

  Delib::ElementType getElementType(const AtomIndexType& index) const;

  const GraphType& getGraph() const;

  const StereocenterList& getStereocenterList() const;

  unsigned getNumAdjacencies(const AtomIndexType& a) const;

  StereocenterList inferStereocentersFromPositions(
    const Delib::PositionCollection& positions
  ) const;

  bool isAdjacent(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

  bool isSafeToRemoveAtom(const AtomIndexType& a) const;

  bool isSafeToRemoveBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

  /*! Returns a range-for temporary object allowing c++11 style for loop 
   * iteration through an atom's adjacencies
   */
  RangeForTemporary<GraphType::adjacency_iterator> iterateAdjacencies(
    const AtomIndexType& a
  ) const;

  /*! Returns a range-for temporary object allowing c++11-style for loop 
   * iteration through edges
   */
  RangeForTemporary<GraphType::edge_iterator> iterateEdges() const;

  /*! Returns a range-for temporary object allowing c++11-style for loop 
   * iteration through edges around a specific atom
   */
  RangeForTemporary<GraphType::out_edge_iterator> iterateEdges(
    const AtomIndexType& a
  ) const;

  RankingInformation rankPriority(
    const AtomIndexType& a,
    const std::set<AtomIndexType>& excludeAdjacent = {}
  ) const;

  unsigned numAtoms() const;

  unsigned numBonds() const;

/* Operators */
  //! Returns the adjacencies of the specified atom index
  RangeForTemporary<GraphType::adjacency_iterator> operator [] (
    const AtomIndexType& a
  ) const;

/* Friends */
  friend struct MoleculeValidator;
};

} // namespace MoleculeManip

std::ostream& operator << (
  std::ostream& os,
  const MoleculeManip::Molecule& molecule
);

#endif
