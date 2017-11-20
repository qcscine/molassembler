#ifndef INCLUDE_MOLECULE_MANIP_MOLECULE_H
#define INCLUDE_MOLECULE_MANIP_MOLECULE_H

#include "Edges.h"
#include "StereocenterList.h"
#include "CycleData.h"
#include "VSEPR.h"

/*! @file
 *
 * Contains the Molecule class declaration, which is the central class of the
 * library.
 */

/* TODO
 */

namespace MoleculeManip {

/*! 
 * Central class of the library, modeling a molecular graph with all state.
 */
class Molecule {
public:
  enum class TemperatureRegimeOption {
    LowTemperature,
    HighTemperature
  };

private:
/* State */
  GraphType _adjacencies;
  StereocenterList _stereocenters;

/* Members */
/* Private members */
  /*! 
   * Adds an vertex to the graph and sets it's element type property, returning
   * the new index.
   */
  AtomIndexType _addAtom(const Delib::ElementType& elementType);

  StereocenterList _detectStereocenters() const;

  /*! Returns if an atom can be a CNStereocenter with multiple assignments
   *
   * Criteria applied are:
   * - Minimum of three adjacent indices
   * - If the high-temperature approximation is invoked, trivalent nitrogen
   *   inverts too rapidly to carry stereoinformation (unless part of a cycle
   *   of size 4 or smaller, where strain hinders inversion)
   */
  bool _isCNStereocenterCandidate(
    const AtomIndexType& atomIndex,
    const TemperatureRegimeOption& temperatureRegime = TemperatureRegimeOption::HighTemperature
  ) const;

  //! Returns whether the specified index is valid or not
  bool _isValidIndex(const AtomIndexType& index) const;

  /*! 
   * Returns a list of edge indices where each endpoint has 1 or two additional
   * substituent besides the edge neighbor
   */
  std::vector<EdgeIndexType> _getEZStereocenterCandidates() const;

  void _pickyFitStereocenter(
    Stereocenters::CNStereocenter& stereocenter,
    const Symmetry::Name& expectedSymmetry,
    const Delib::PositionCollection& positions
  ) const;

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
  explicit Molecule(const GraphType& graph);

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

  /*! Sets the stereocenter assignment at a particular atom
   *
   * This sets the stereocenter assignment at a specific atom index. For this,
   * a stereocenter must be instantiated and contained in the StereocenterList
   * returned by getStereocenterList(). The supplied assignment must be either
   * boost::none or smaller than stereocenterPtr->numAssignments().
   *
   * NOTE: Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and 
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenterAtAtom(
    const AtomIndexType& a,
    const boost::optional<unsigned>& assignment
  );

  //! Changes an existing atom's element type
  void changeElementType(
    const AtomIndexType& a,
    const Delib::ElementType& elementType
  );

  //! TODO Temporary function prior to proper editing state correctness
  void refreshStereocenters();

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

  // TODO implement, with trigger to keep valid state (can change ranking)
  /*! Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently: In case the symmetry does
   * not fit the number of substituents at that atom, this function throws. If
   * there is already a CNStereocenter instantiated at this atom index, its
   * underlying symmetry is altered. If there is an EZStereocenter
   * instantiated on this atom, this function will throw. If there is no
   * CNStereocenter at this index, one is instantiated. In all cases, new or
   * modified stereocenters are default-assigned if there is only one possible
   * assignment.
   */
  void setGeometryAtAtom(
    const AtomIndexType& a,
    const Symmetry::Name& symmetryName
  );

/* Information */

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

  unsigned numAtoms() const;

  unsigned numBonds() const;

  RankingInformation rankPriority(
    const AtomIndexType& a,
    const std::set<AtomIndexType>& excludeAdjacent = {},
    const boost::optional<Delib::PositionCollection>& positionsOption = boost::none
  ) const;

/* Operators */
  //! Returns the adjacencies of the specified atom index
  RangeForTemporary<GraphType::adjacency_iterator> operator [] (
    const AtomIndexType& a
  ) const;

/* Friends */
  friend struct MoleculeValidator;
  friend class RankingTree;
};

} // namespace MoleculeManip

std::ostream& operator << (
  std::ostream& os,
  const MoleculeManip::Molecule& molecule
);

#endif
