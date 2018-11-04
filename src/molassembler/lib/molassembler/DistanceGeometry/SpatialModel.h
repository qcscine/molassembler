// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_SPATIAL_MODEL_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_SPATIAL_MODEL_H

#include "molassembler/Conformers.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"

#include <cmath>

/*! @file
 *
 * @brief Convert a Molecule to atom-pair distance bounds and chiral constraints
 *
 * The molecular graph with all its conformational specifications via
 * stereopermutator assignments must be transformed into a spatial model that
 * describes its internal degrees of freedom in a manner translatable to
 * pairwise distance bounds for the distance geometry algorithm. This file
 * contains the class declaration for that spatial model.
 */

namespace molassembler {

namespace DistanceGeometry {

template<int size>
std::array<AtomIndex, size> orderedIndexSequence(
  const std::array<AtomIndex, size>& source
) {
  if(source.front() > source.back()) {
    auto copy = source;
    std::reverse(
      copy.begin(),
      copy.end()
    );
    return copy;
  }

  return source;
}

template<int size>
std::array<AtomIndex, size> orderedIndexSequence(
  std::initializer_list<AtomIndex>& initializer
) {
  assert(size == initializer.size());
  return orderedIndexSequence(
    std::array<AtomIndex, size>(initializer)
  );
}

template<typename ... Inds>
auto orderedSequence(Inds ... inds) {
  std::array<AtomIndex, sizeof...(inds)> indices {{
    static_cast<AtomIndex>(inds)...
  }};

  if(indices.front() > indices.back()) {
    std::reverse(indices.begin(), indices.end());
  }

  return indices;
}

/*!
 * Keeps a record of the internal coordinate bounds that a molecular graph is
 * interpreted as and generates a list of atom-pairwise bounds from which
 * classes performing smoothing can be constructed.
 */
class SpatialModel {
public:
//!@name Member types
//!@{
  struct ModelGraphWriter;

  using BoundsList = std::map<
    std::array<AtomIndex, 2>,
    ValueBounds
  >;
//!@}

//!@name Static members
//!@{
  /*! Relative bond distance variance, 0.0x meaning x% variance. Must fulfill
   * 0 < x << 1
   */
  static constexpr double bondRelativeVariance = 0.01;
  //! Absolute angle variance in radians. Must fulfill 0 < x << M_PI
  static constexpr double angleAbsoluteVariance = M_PI / 90; // ~ 2°
  //! Absolute dihedral angle variance in radians.
  static constexpr double dihedralAbsoluteVariance = M_PI / 90; // ~ 2°

  static constexpr ValueBounds angleClampBounds {0.0, M_PI};
  static ValueBounds defaultDihedralBounds;

/* Static functions */
  static boost::optional<ValueBounds> coneAngle(
    const std::vector<AtomIndex>& baseConstituents,
    const ValueBounds& coneHeightBounds,
    double bondRelativeVariance,
    const OuterGraph& graph,
    const Cycles& etaLessCycles
  );

  static double spiroCrossAngle(double alpha, double beta);

  static ValueBounds ligandDistanceFromCenter(
    const std::vector<AtomIndex>& ligandIndices,
    AtomIndex centralIndex,
    double bondRelativeVariance,
    const OuterGraph& graph
  );

  /**
   * @brief Yields central value plus/minus some absolute variance bounds
   *
   * @param centralValue The central value
   * @param absoluteVariance The value to subtract and add to yield new bounds
   *
   * @return bounds with median central value and width 2 * absoluteVariance
   */
  static ValueBounds makeBoundsFromCentralValue(
    double centralValue,
    double absoluteVariance
  );

  /**
   * @brief Analogous to C++17's clamp, reduces ValueBounds to the maximal
   *   extent specified by some clamping bounds
   *
   * @param bounds The bounds to be clamped
   * @param clampBounds The range in which bounds may exist
   *
   * @return New bounds that are within the clamp bounds
   */
  static ValueBounds clamp(
    const ValueBounds& bounds,
    const ValueBounds& clampBounds
  );

  /**
   * @brief Ensures that the preconditions for fixed positions in Distance
   *   Geometry detailed in the Configuration object are met
   *
   * @param molecule The molecule for which the Configuration is valid
   * @param configuration The Configuration object detailing fixed positions
   *
   * @throws std::runtime_error If the preconditions are not met
   */
  static void checkFixedPositionsPreconditions(
    const Molecule& molecule,
    const Configuration& configuration
  );
//!@}

//!@name Special member functions
//!@{
  /**
   * @brief Model a molecule into internal coordinate bounds stored internally
   *
   * Main space where modeling is performed. This is where a molecule's
   * graph and list of stereopermutators are transformed into explicit bounds
   * on 1-2 (bonds), 1-3 (angles) and 1-4 (dihedral) internal coordinates. This
   * determines which conformations are accessible.
   *
   * @param molecule The molecule that is to be modeled. This may not contain
   *   stereopermutators with zero assignments or unassigned stereopermutators.
   * @param configuration The Distance Geometry configuration object. Relevant
   *   for this stage of the process are fixed positions, if set.
   * @param looseningMultiplier A multiplier that loosens all bounds set upon
   *   internal coordinates. Affects the entire conformation.
   */
  SpatialModel(
    const Molecule& molecule,
    const Configuration& configuration,
    double looseningMultiplier = 1.0
  );
//!@}

//!@name Modifiers
//!@{
  //! Sets the bond bounds to the model.
  void setBondBoundsIfEmpty(
    const std::array<AtomIndex, 2>& bondIndices,
    double centralValue
  );

  //! Sets bond bounds to exact value bounds.
  void setBondBoundsIfEmpty(
    const std::array<AtomIndex, 2>& bondIndices,
    const ValueBounds& bounds
  );

  /*!
   * Adds the angle bounds to the model, but only if the information for that
   * set of indices does not exist yet.
   */
  void setAngleBoundsIfEmpty(
    const std::array<AtomIndex, 3>& angleIndices,
    const ValueBounds& bounds
  );

  /*!
   * Adds the dihedral bounds to the model, but only if the information for that
   * set of indices does not exist yet.
   */
  void setDihedralBoundsIfEmpty(
    const std::array<AtomIndex, 4>& dihedralIndices,
    const ValueBounds& bounds
  );

  void addAtomStereopermutatorInformation(
    const AtomStereopermutator& permutator,
    const std::function<double(const AtomIndex)>& cycleMultiplierForIndex,
    double looseningMultiplier,
    const std::unordered_map<AtomIndex, Delib::Position>& fixedAngstromPositions
  );

  void addBondStereopermutatorInformation(
    const BondStereopermutator& permutator,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    double looseningMultiplier,
    const std::unordered_map<AtomIndex, Delib::Position>& fixedAngstromPositions
  );

  //! Adds [0, π] default angle bounds for all bonded atom triples
  void addDefaultAngles();

  /*! @brief Adds [0, π] default dihedrals to the model (of connected sequences)
   *
   * Adds [0, π] default dihedrals to the model. Avoids treating 1-4 pairs as
   * nonbonded.
   */
  void addDefaultDihedrals();
//!@}

//!@name Information
//!@{
  boost::optional<ValueBounds> coneAngle(
    const std::vector<AtomIndex>& ligandIndices,
    const ValueBounds& coneHeightBounds
  ) const;

  void dumpDebugInfo() const;

  ValueBounds ligandDistance(
    const std::vector<AtomIndex>& ligandIndices,
    AtomIndex centralIndex
  ) const;

  std::vector<ChiralityConstraint> getChiralityConstraints() const;

  BoundsList makeBoundsList() const;

  std::string dumpGraphviz() const;

  void writeGraphviz(const std::string& filename) const;
//!@}

private:
  // Molecule closure
  const Molecule& _molecule;

  // Mutable state
  double _looseningMultiplier;
  StereopermutatorList _stereopermutators;

  //! Constraints by fixed positions
  std::map<
    std::array<AtomIndex, 2>,
    ValueBounds
  > _constraints;

  //!@name Bounds on internal angles
  //!@{
  std::map<
    std::array<AtomIndex, 2>,
    ValueBounds
  > _bondBounds;
  std::map<
    std::array<AtomIndex, 3>,
    ValueBounds
  > _angleBounds;
  std::map<
    std::array<AtomIndex, 4>,
    ValueBounds
  > _dihedralBounds;
  //!@}

  //! Chiral constraints
  std::vector<ChiralityConstraint> _chiralConstraints;
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
