/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Convert a Molecule to atom-pair distance bounds and chiral constraints
 *
 * The molecular graph with all its conformational specifications via
 * stereopermutator assignments must be transformed into a spatial model that
 * describes its internal degrees of freedom in a manner translatable to
 * pairwise distance bounds for the distance geometry algorithm. This file
 * contains the class declaration for that spatial model.
 *
 * @todo With the change in default representation of a Molecule in that
 * AtomStereopermutators basically must exist on any non-terminal atom, it
 * may be possible to remove the StereopermutatorList member of SpatialModel.
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_SPATIAL_MODEL_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_SPATIAL_MODEL_H

#include "molassembler/Conformers.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"

#include <cmath>

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

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
  /*!
   * @brief Relative bond distance variance
   * 0.0x mean x% variance. Must fulfill 0 < x << 1
   */
  static constexpr double bondRelativeVariance = 0.01;
  //! Absolute angle variance in radians. Must fulfill 0 < x << M_PI
  static constexpr double angleAbsoluteVariance = M_PI / 90; // ~ 2°
  //! Absolute dihedral angle variance in radians.
  static constexpr double dihedralAbsoluteVariance = M_PI / 90; // ~ 2°

  //! Defines clamping bounds on angles
  static constexpr ValueBounds angleClampBounds {0.0, M_PI};
  //! Defines the range (-π,π] using std::nextafter (not constexpr)
  static ValueBounds defaultDihedralBounds;

/* Static functions */
  /**
   * @brief Tries to calculate the cone angle for a possibly haptic ligand site
   *
   * @param baseConstituents Atom indices constituting the ligand site
   * @param coneHeightBounds Bounds on the distance of the ligand site to
   *   the central atom
   * @param graph The molecular graph to model
   * @param etaLessCycles Reference to a cycle interpretation object that
   *   ignored eta bonds
   *
   * @return If calculable, bounds on the cone angle spanned by the ligand site
   */
  static boost::optional<ValueBounds> coneAngle(
    const std::vector<AtomIndex>& baseConstituents,
    const ValueBounds& coneHeightBounds,
    const OuterGraph& graph,
    const Cycles& etaLessCycles
  );

  /**
   * @brief Calculates the cross angle between opposite cycle atoms in a spirocenter
   *
   * @param alpha The left cycle-internal angle
   * @param beta The right cycle-internal angle
   *
   * @return The cross angle between cycle atoms in the spirocenter
   */
  static double spiroCrossAngle(double alpha, double beta);

  /**
   * @brief Calculates bounds on the distance of a (possibly haptic) ligand
   *   site to a central index
   *
   * @param siteAtomList Atom indices constituting the ligand site
   * @param centralIndex The central index to which the ligand is bound
   * @param graph The molecular graph to model
   *
   * @return Bounds on the distance of the ligand site to the central index
   */
  static ValueBounds siteDistanceFromCenter(
    const std::vector<AtomIndex>& siteAtomList,
    AtomIndex centralIndex,
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
    ValueBounds bounds,
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
   *   for this stage of the process are the loosening multiplier and fixed
   *   positions, if set.
   */
  SpatialModel(
    const Molecule& molecule,
    const Configuration& configuration
  );
//!@}

//!@name Modifiers
//!@{
  /*!
   * @brief Sets bond bounds to exact value bounds if unset
   *
   * @param bondIndices The atom pair on which to place bond distance bounds
   * @param bounds The distance bounds to enforce
   * @pre @p bondIndices must be ordered (i.e. front() < back())
   */
  void setBondBoundsIfEmpty(
    std::array<AtomIndex, 2> bondIndices,
    ValueBounds bounds
  );

  /*!
   * @brief Adds the angle bounds to the model, but only if the information for that
   *   set of indices does not exist yet.
   *
   * @param angleIndices The atom index sequence specifying an angle
   * @param bounds The angle bounds to enforce
   * @pre @p angleIndices must be ordered (i.e. front() < back())
   *
   * @note Passed angle bounds are clamped against the angleClampBounds
   */
  void setAngleBoundsIfEmpty(
    std::array<AtomIndex, 3> angleIndices,
    ValueBounds bounds
  );

  /*!
   * @brief Adds the dihedral bounds to the model, but only if the information for that
   *   set of indices does not exist yet.
   *
   * @param dihedralIndices The atom index sequence specifying a dihedral
   * @param bounds The dihedral bounds to enforce
   * @pre @p dihedralIndices must be ordered (i.e. front() < back())
   */
  void setDihedralBoundsIfEmpty(
    std::array<AtomIndex, 4> dihedralIndices,
    ValueBounds bounds
  );

  /**
   * @brief Adds angle information to the internal coordinate bounds and
   *   collects chiral constraints
   *
   * @param permutator The AtomStereopermutator to collect information from
   * @param cycleMultiplierForIndex A function yielding factors with which to
   *   multiply angular variances depending on the smalles cycle an atom is a
   *   member of
   * @param looseningMultiplier A loosening factor for the overall model
   * @param fixedAngstromPositions A mapping between atom indices and fixed
   *   spatial positions
   */
  void addAtomStereopermutatorInformation(
    const AtomStereopermutator& permutator,
    const std::function<double(const AtomIndex)>& cycleMultiplierForIndex,
    double looseningMultiplier,
    const std::unordered_map<AtomIndex, Scine::Utils::Position>& fixedAngstromPositions
  );

  /**
   * @brief Adds dihedral information to the internal coordinate bounds and
   *   collects dihedral constraints
   *
   * @param permutator The BondStereopermutator to collect information from
   * @param stereopermutatorA One AtomStereopermutator constituting the BondStereopermutator
   * @param stereopermutatorB The other AtomStereopermutator constituting the BondStereopermutator
   * @param looseningMultiplier A loosening factor for the overall model
   * @param fixedAngstromPositions A mapping between atom indices and fixed
   *   spatial positions
   */
  void addBondStereopermutatorInformation(
    const BondStereopermutator& permutator,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    double looseningMultiplier,
    const std::unordered_map<AtomIndex, Scine::Utils::Position>& fixedAngstromPositions
  );

  //! Adds [0, π] default angle bounds for all bonded atom triples
  void addDefaultAngles();

  /*!
   * @brief Adds [0, π] default dihedrals to the model (of connected sequences)
   *
   * Adds [0, π] default dihedrals to the model. Avoids treating 1-4 pairs as
   * nonbonded.
   */
  void addDefaultDihedrals();
//!@}

//!@name Information
//!@{
  //! Yields all collected chiral constraints
  std::vector<ChiralityConstraint> getChiralityConstraints() const;

  //! Yields all collected dihedral constraints
  std::vector<DihedralConstraint> getDihedralConstraints() const;

  /**
   * @brief Generates a list of atom-pairwise distance bounds from the internal
   *   coordinate bounds and fixed positions
   *
   * @return A list of atom-pairwise distance bounds
   */
  BoundsList makeBoundsList() const;

  /**
   * @brief Generates a string graphviz representation of the modeled molecule
   *
   * The graph contains basic connectivity, stereopermutator information
   * and internal coordinate bounds.
   *
   * @return A string that can be converted into an image using graphviz of
   *   the molecule being modeled.
   */
  std::string dumpGraphviz() const;

  /**
   * @brief Writes a graphviz representation of a modeled molecule to a file
   *
   * The graph contains basic connectivity, stereopermutator information
   * and internal coordinate bounds.
   *
   * @param filename The filename to which to write the graphviz representation
   */
  void writeGraphviz(const std::string& filename) const;
//!@}

private:
  // Molecule closure
  const Molecule& _molecule;

  // Mutable state
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
  std::vector<DihedralConstraint> _dihedralConstraints;
};

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
