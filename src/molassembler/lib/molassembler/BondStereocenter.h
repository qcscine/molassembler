#ifndef INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H
#define INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H

#include "stereopermutation/Composites.h"

#include "detail/AngstromWrapper.h"
#include "detail/SharedTypes.h"

/*! @file
 *
 * Contains the BondStereocenter class declaration, which models E/Z double bond
 * stereocenters in molecules.
 */

/* TODO
 * - Various TODO in implementation
 * - Consistency in 0 = E, 1 = Z (or other) must be established somehow (using
 *   ligandsRanking)
 * - State propagation is going to be a nasty piece of work
 *   Will need to store more data, using AtomState, for both, but probably also
 *   including the symmetryPositionMap and the assignment
 */

namespace molassembler {

// Forward-declarations
class AtomStereocenter;

namespace DistanceGeometry {
class SpatialModel;
struct ChiralityConstraint;
} // namespace DistanceGeometry

class BondStereocenter {
public:
  static constexpr double chiralityConstraintTolerance = 0.1;

  struct AtomState {
    AtomIndexType index;
    Symmetry::Name symmetry;

    bool operator == (const AtomState& other) const;
    bool operator != (const AtomState& other) const;
  };

private:
  stereopermutation::Composite _composite;
  GraphType::edge_descriptor _edge;
  boost::optional<unsigned> _assignment;
  AtomState _left, _right;

  static std::vector<char> _charifyRankedLigands(
    const std::vector<std::vector<unsigned>> ligandsRanking
  );

public:
  BondStereocenter() = delete;
  /*!
   * @warning No haptic ligands allowed.
   * @note The left stereocenter must have the smaller
   */
  BondStereocenter(
    const AtomStereocenter& leftStereocenter,
    const AtomStereocenter& rightStereocenter,
    const GraphType::edge_descriptor edge
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom();

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereocenter& left,
    const AtomStereocenter& right
  );

/* Information */
  boost::optional<unsigned> assigned() const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    const double looseningMultiplier,
    const AtomStereocenter& left,
    const AtomStereocenter& right
  ) const;

  std::string info() const;

  std::string rankInfo() const;

  GraphType::edge_descriptor edge() const;

  const AtomState& left() const {
    return _left;
  }

  const AtomState& right() const {
    return _right;
  }

  void setModelInformation(
    DistanceGeometry::SpatialModel& model,
    const AtomStereocenter& left,
    const AtomStereocenter& right,
    const double looseningMultiplier
  ) const;

/* Operators */
  bool operator == (const BondStereocenter& other) const;
  bool operator != (const BondStereocenter& other) const;
};

} // namespace molassembler

#endif
