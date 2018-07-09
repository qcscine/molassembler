#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_SPATIAL_MODEL_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_SPATIAL_MODEL_H

#include "Molecule.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"

/*! @file
 *
 * The molecular graph with all its conformational specifications via
 * stereocenter assignments must be transformed into a spatial model that
 * describes its internal degrees of freedom in a manner translatable to
 * pairwise distance bounds for the distance geometry algorithm. This file
 * contains the class declaration for that spatial model.
 */

namespace molassembler {

namespace DistanceGeometry {

template<int size>
std::array<AtomIndexType, size> orderedIndexSequence(
  const std::array<AtomIndexType, size>& source
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
std::array<AtomIndexType, size> orderedIndexSequence(
  std::initializer_list<AtomIndexType>& initializer
) {
  assert(size == initializer.size());
  return orderedIndexSequence(
    std::array<AtomIndexType, size>(initializer)
  );
}

template<typename ... Inds>
auto orderedSequence(Inds ... inds) {
  std::array<AtomIndexType, sizeof...(inds)> indices {{
    static_cast<AtomIndexType>(inds)...
  }};

  if(indices.front() > indices.back()) {
    std::reverse(indices.begin(), indices.end());
  }

  return indices;
}

/*!
 * Keeps a record of the internal dimension bounds that a molecular graph is
 * interpreted as and permits the generation of a distance bounds matrix.
 */
class SpatialModel {
public:
/* Typedefs */
  struct ModelGraphWriter;

private:
  // Closures
  const Molecule& _molecule;

  // Mutable state

  double _looseningMultiplier;

  std::map<
    std::array<AtomIndexType, 2>,
    ValueBounds
  > _bondBounds;
  std::map<
    std::array<AtomIndexType, 3>,
    ValueBounds
  > _angleBounds;
  std::map<
    std::array<AtomIndexType, 4>,
    ValueBounds
  > _dihedralBounds;

  std::map<
    AtomIndexType,
    std::shared_ptr<Stereocenters::Stereocenter>
  > _stereocenterMap;

public:
/* Static constants */
  /*! Relative bond distance variance, 0.0x meaning x% variance. Must fulfill
   * 0 < x << 1
   */
  static constexpr double bondRelativeVariance = 0.01;
  //! Absolute angle variance in radians. Must fulfill 0 < x << M_PI
  static constexpr double angleAbsoluteVariance = M_PI / 90; // ~ 2°
  //! Absolute dihedral angle variance in radians.
  static constexpr double dihedralAbsoluteVariance = M_PI / 36; // ~ 5°

  static const ValueBounds angleClampBounds;
  static const ValueBounds dihedralClampBounds;

/* Static functions */
  static boost::optional<ValueBounds> coneAngle(
    const std::vector<AtomIndexType>& baseConstituents,
    const ValueBounds& coneHeightBounds,
    const double bondRelativeVariance,
    const GraphType& graph,
    const Cycles& etaLessCycles
  );

  static double spiroCrossAngle(const double alpha, const double beta);

  static ValueBounds ligandDistanceFromCenter(
    const std::vector<AtomIndexType>& ligandIndices,
    const AtomIndexType centralIndex,
    const double bondRelativeVariance,
    const GraphType& graph
  );

  static ValueBounds makeBoundsFromCentralValue(
    const double centralValue,
    const double absoluteVariance
  );

  static ValueBounds clamp(
    const ValueBounds& bounds,
    const ValueBounds& clampBounds
  );

/* Constructor */
  SpatialModel(
    const Molecule& molecule,
    const double looseningMultiplier = 1.0
  );

/* Modification */
  //! Sets the bond bounds to the model.
  void setBondBoundsIfEmpty(
    const std::array<AtomIndexType, 2>& bondIndices,
    const double centralValue
  );

  //! Sets bond bounds to exact value bounds.
  void setBondBoundsIfEmpty(
    const std::array<AtomIndexType, 2>& bondIndices,
    const ValueBounds& bounds
  );

  /*!
   * Adds the angle bounds to the model, but only if the information for that
   * set of indices does not exist yet.
   */
  void setAngleBoundsIfEmpty(
    const std::array<AtomIndexType, 3>& angleIndices,
    const ValueBounds& bounds
  );

  /*!
   * Adds the dihedral bounds to the model, but only if the information for that
   * set of indices does not exist yet.
   */
  void setDihedralBoundsIfEmpty(
    const std::array<AtomIndexType, 4>& dihedralIndices,
    const ValueBounds& bounds
  );

  //! Adds [0, 2π] default angle bounds for all bonded atom triples
  void addDefaultAngles();

  /*!
   * Adds [0, 2π] default dihedrals to the model. Use immediately before
   * calling makeDistanceBounds if you want default dihedrals modeled in the
   * distance bounds as well. In principle, the default dihedral distances are
   * inferable from the existing information using bound smoothing, but this
   * fashion is probably significantly faster.
   */
  void addDefaultDihedrals();

  using BoundsList = std::map<
    std::array<AtomIndexType, 2>,
    ValueBounds
  >;

/* Information */

  boost::optional<ValueBounds> coneAngle(
    const std::vector<AtomIndexType>& ligandIndices,
    const ValueBounds& coneHeightBounds
  ) const;

  void dumpDebugInfo() const;

  ValueBounds ligandDistance(
    const std::vector<AtomIndexType>& ligandIndices,
    const AtomIndexType centralIndex
  ) const;


  std::vector<ChiralityConstraint> getChiralityConstraints() const;

  [[deprecated]]
  DistanceBoundsMatrix makeBounds() const;

  BoundsList makeBoundsList() const;

  std::string dumpGraphviz() const;

  void writeGraphviz(const std::string& filename) const;
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
