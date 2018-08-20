#ifndef INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H
#define INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H

#include "detail/SharedTypes.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

/*! @file
 *
 * Contains the BondStereocenter class declaration, which models E/Z double bond
 * stereocenters in molecules.
 */

/* TODO
 * - Various TODO in implementation
 * - Consistency in 0 = E, 1 = Z (or other) must be established somehow (using
 *   ligandsRanking hopefully)
 * - State propagation is going to be a nasty piece of work
 */

namespace molassembler {

struct AngstromWrapper;

// Forward-declarations
class AtomStereocenter;

namespace DistanceGeometry {
class SpatialModel;
struct ChiralityConstraint;
} // namespace DistanceGeometry

class BondStereocenter {
public:
  static constexpr double chiralityConstraintTolerance = 0.1;
  static constexpr double assignmentAcceptanceDihedralThreshold = M_PI * 3.0 / 180.0;

public:
  /* Rule of five members */
  BondStereocenter(BondStereocenter&& other) noexcept;
  BondStereocenter& operator = (BondStereocenter&& other) noexcept;
  BondStereocenter(const BondStereocenter& other);
  BondStereocenter& operator = (const BondStereocenter& other);
  ~BondStereocenter();

  BondStereocenter() = delete;
  /*!
   * @warning No haptic ligands allowed.
   * @note The left stereocenter must have the smaller
   */
  BondStereocenter(
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB,
    const GraphType::edge_descriptor edge
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom();

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB
  );

/* Information */
  boost::optional<unsigned> assigned() const;

  bool hasSameCompositeOrientation(const BondStereocenter& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    double looseningMultiplier,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB
  ) const;

  std::string info() const;

  std::string rankInfo() const;

  GraphType::edge_descriptor edge() const;

  void setModelInformation(
    DistanceGeometry::SpatialModel& model,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB,
    double looseningMultiplier
  ) const;

/* Operators */
  bool operator == (const BondStereocenter& other) const;
  bool operator != (const BondStereocenter& other) const;

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

#endif
