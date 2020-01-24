/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/Molecule.h"
#include "molassembler/Graph.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/AtomStereopermutator.h"

#include "boost/functional/hash.hpp"
#include <unordered_map>

namespace Scine {
namespace molassembler {
namespace distance_geometry {

/*! @brief Class performing spatial modeling of molecules
 *
 * Keeps a record of the internal coordinate bounds that a molecular graph is
 * interpreted as and generates a list of atom-pairwise bounds from which
 * classes performing smoothing can be constructed.
 */
class SpatialModel {
public:
//!@name Member types
//!@{
  //! Type that can be used to write a graphviz representation of the model
  struct ModelGraphWriter;

  //! ValueBounds container for various internal coordinates
  template<unsigned size>
  using BoundsMapType = std::unordered_map<
    std::array<AtomIndex, size>,
    ValueBounds,
    boost::hash<
      std::array<AtomIndex, size>
    >
  >;

  //! Type used to store fixed positions in angstrom
  using FixedPositionsMapType = std::unordered_map<AtomIndex, Utils::Position>;

  /*!
   * @brief Type used to represent atom-pairwise distance bounds constructed
   *   from bounds on internal coordinates
   *
   * Diagonal entries are undefined. Strict lower triangle contains lower
   * bounds, strict upper triangle contains upper bounds.
   */
  using BoundsMatrix = Eigen::MatrixXd;
//!@}

//!@name Static members
//!@{
  /*!
   * @brief Relative bond distance variance
   * 0.0x mean x% variance. Must fulfill 0 < x << 1
   */
  static constexpr double bondRelativeVariance = 0.01;
  /*! @brief Relative angle variance
   *
   * 0.0x mean x& variance. Must fulfill 0 < x << 1
   */
  static constexpr double angleRelativeVariance = 0.02;
  //! Absolute dihedral angle variance in radians.
  static constexpr double dihedralAbsoluteVariance = M_PI / 90; // ~ 2°

  //! Defines clamping bounds on angles
  static constexpr ValueBounds angleClampBounds {0.0, M_PI};
  //! Defines the range (-π,π] using std::nextafter (not constexpr)
  static const ValueBounds defaultDihedralBounds;

/* Static functions */
  /*! @brief Models the equilibrium distance between two bonded atoms
   *
   * @complexity{@math{\Theta(1)}}
   * @pre @p i and @p j are bonded
   */
  static double modelDistance(AtomIndex i, AtomIndex j, const PrivateGraph& graph);
  //! @overload
  static double modelDistance(const BondIndex& bond, const PrivateGraph& graph);

  /*! @brief Tries to find a cycle that consists of a particular set of atoms
   *
   * @complexity{@math{\Theta(C)} where @math{C} is the number of cycles in the
   * molecule}
   *
   * @param atoms The exact set of atoms that constitute the cycle.
   * @param graph The graph context
   *
   * @returns A vector of bond indices of exactly the matching cycle if found.
   * An empty vector is returned otherwise.
   */
  static std::vector<BondIndex> cycleConsistingOfExactly(
    const std::vector<AtomIndex>& atoms,
    const PrivateGraph& graph
  );

  /** @brief Tries to calculate the cone angle for a possibly haptic ligand site
   *
   * @complexity{@math{O(C)} where @math{C} is the number of cycles in the
   * molecule}
   *
   * @param baseConstituents Atom indices constituting the ligand site
   * @param coneHeightBounds Bounds on the distance of the ligand site to
   *   the central atom
   * @param graph The molecular graph to model
   *
   * @return If calculable, bounds on the cone angle spanned by the ligand site
   */
  static boost::optional<ValueBounds> coneAngle(
    const std::vector<AtomIndex>& baseConstituents,
    const ValueBounds& coneHeightBounds,
    const Graph& graph
  );

  /** @brief Calculates the cross angle between opposite cycle atoms in a spirocenter
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param alpha The left cycle-internal angle
   * @param beta The right cycle-internal angle
   *
   * @return The cross angle between cycle atoms in the spirocenter
   */
  static double spiroCrossAngle(double alpha, double beta);

  /** @brief Calculates bounds on the distance of a (possibly haptic) ligand
   *   site to a central index
   *
   * @complexity{@math{\Theta(A)} where @math{A} is the number of atoms comprising
   * the site}
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
    const Graph& graph
  );

  /** @brief Yields central value plus/minus some absolute variance bounds
   *
   * @complexity{@math{\Theta(1)}}
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

  /*! @brief Idealization strictness loosening multiplier due to small cycles
   *
   * @todo complexity
   *
   * @returns A multiplier intended for the absolute angle variance for a
   *   particular atom. If that atom is part of a cycle of size < 6, the
   *   multiplier is > 1.
   */
  static double smallestCycleDistortionMultiplier(
    const AtomIndex i,
    const Cycles& cycles
  );

  /*
   * @brief Generates a matrix of atom-pairwise distance bounds
   *
   * @complexity{@math{O(P_2 + P_3 + P_4)} where @math{P_i} is the number of
   * distinct paths of length @math{i} in the graph. That should scale at least
   * linearly in the number of vertices.}
   *
   * @param N size of the system
   * @param fixedPositionbounds Distances to enforce due to fixed positions
   * @param bondBounds Distances to enforce for bonds
   * @param angleBounds Angle bounds to enforce on atom index triples
   * @param dihedralBounds Dihedral bounds to enforce on atom index quadruplet
   *
   * @return A matrix of atom-pairwise distance bounds. Lower bounds are in the
   * strict lower triangle of the matrix, upper bounds in the strict upper
   * triangle.
   */
  static BoundsMatrix makePairwiseBounds(
    unsigned N,
    const BoundsMapType<2>& fixedPositionBounds,
    const BoundsMapType<2>& bondBounds,
    const BoundsMapType<3>& angleBounds,
    const BoundsMapType<4>& dihedralBounds
  );

  /** @brief Determines the central value of the angle between
   *   AtomStereopermutator sites
   *
   * @complexity{Varies. For most cases, @math{\Omega(1)}}
   *
   * @param centralIndex The atom index of the central index of the angle
   * @param shape The local shape at @p centralIndex
   * @param ranking The ranking of substituents at @p centralIndex
   * @param shapeVertexMap The mapping from site indices to shape vertices
   * @param sites The two sites between which the angle is to be determined
   * @param inner Graph instance being modeled
   *
   * @return The central value of the angle between the sites
   */
  static double siteCentralAngle(
    AtomIndex centralIndex,
    const shapes::Shape& shape,
    const RankingInformation& ranking,
    const std::vector<unsigned>& shapeVertexMap,
    const std::pair<unsigned, unsigned>& sites,
    const PrivateGraph& inner
  );

  //! @overload
  static double siteCentralAngle(
    const AtomStereopermutator& permutator,
    const std::pair<unsigned, unsigned>& sites,
    const PrivateGraph& inner
  );

  /** @brief Models bounds on the angle between AtomStereopermutator sites
   *
   * @complexity{Varies. For most cases, @math{\Omega(1)}}
   *
   * @param permutator The permutator being modeled
   * @param sites The two sites between which the angle is to be determined
   * @param looseningMultiplier Loosening of that particular atom
   * @param inner Graph instance being modeled
   *
   * @return Bounds on the angle between two AtomStereopermutator sites
   */
  static ValueBounds modelSiteAngleBounds(
    const AtomStereopermutator& permutator,
    const std::pair<unsigned, unsigned>& sites,
    const double looseningMultiplier,
    const PrivateGraph& inner
  );

  /** @brief Creates a volume-specific chiral constraint from a list of
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param minimalConstraint Site index sequence defining a chiral constraint
   * @param permutator The AtomStereopermutator from which the minimalConstraint
   *   has been generated
   * @param looseningMultiplier Describes the spatial modeling loosening of
   *   bounds
   *
   * @return A volume constraint on the site index sequence in a final
   *   conformation
   */
  static ChiralConstraint makeChiralConstraint(
    const AtomStereopermutator::MinimalChiralConstraint& minimalConstraint,
    const AtomStereopermutator& permutator,
    const double looseningMultiplier
  );

  /** @brief Analogous to C++17's clamp, reduces ValueBounds to the maximal
   *   extent specified by some clamping bounds
   *
   * @complexity{@math{\Theta(1)}}
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
   * @complexity{At least @math{O(P_2 + P_3 + P_4)} where @math{P_i} is the
   * number of distinct paths of length @math{i} in the graph. That should
   * scale at least linearly in the number of vertices.}
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
  /*! @brief Sets bond bounds to exact value bounds if unset
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param bondIndices The atom pair on which to place bond distance bounds
   * @param bounds The distance bounds to enforce
   * @pre @p bondIndices must be ordered (i.e. front() < back())
   */
  void setBondBoundsIfEmpty(
    std::array<AtomIndex, 2> bondIndices,
    ValueBounds bounds
  );

  /*! @brief Adds angle bounds to the model if unset
   *
   * @complexity{@math{\Theta(1)}}
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
   * @brief Adds the dihedral bounds to the model if unset
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param dihedralIndices The atom index sequence specifying a dihedral
   * @param bounds The dihedral bounds to enforce
   * @pre @p dihedralIndices must be ordered (i.e. front() < back())
   */
  void setDihedralBoundsIfEmpty(
    std::array<AtomIndex, 4> dihedralIndices,
    ValueBounds bounds
  );

  /** @brief Adds angle information to the internal coordinate bounds and
   *   collects chiral constraints
   *
   * @complexity{@math{O(S^2)} where @math{S} is the size of the modeled
   * shape}
   *
   * @param permutator The AtomStereopermutator to collect information from
   * @param graph The graph context
   * @param looseningMultiplier A loosening factor for the overall model
   * @param fixedAngstromPositions A mapping between atom indices and fixed
   *   spatial positions
   * @param forceChiralConstraintEmission If set, the stereopermutator
   *   emits chiral constraints even if the stereopermutator does not have chiral
   *   state, i.e. numAssignments() <= 1, as long as it is assigned.
   */
  void addAtomStereopermutatorInformation(
    const AtomStereopermutator& permutator,
    const PrivateGraph& graph,
    double looseningMultiplier,
    const std::unordered_map<AtomIndex, Utils::Position>& fixedAngstromPositions,
    bool forceChiralConstraintEmission
  );

  /** @brief Adds dihedral information to the internal coordinate bounds and
   *   collects dihedral constraints
   *
   * @complexity{@math{O(S^2)} where @math{S} is the size of the larger modeled
   * shape}
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
    const std::unordered_map<AtomIndex, Utils::Position>& fixedAngstromPositions
  );
//!@}

//!@name Information
//!@{
  /*! @brief Yields all collected chiral constraints
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::vector<ChiralConstraint> getChiralConstraints() const;

  /*! @brief Yields all collected dihedral constraints
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::vector<DihedralConstraint> getDihedralConstraints() const;

  /** @brief Generate unsmoothed atom-pairwise distance bounds matrix
   *
   * Generates a matrix of atom-pairwise distance bounds from the
   * internal coordinate bounds and fixed positions from which this
   * was constructed.
   *
   * @complexity{@math{O(P_2 + P_3 + P_4)} where @math{P_i} is the number of
   * distinct paths of length @math{i} in the graph. That should scale at least
   * linearly in the number of vertices.}
   *
   * @return A matrix of atom-pairwise distance bounds. Lower bounds are in the
   * strict lower triangle of the matrix, upper bounds in the strict upper
   * triangle.
   */
  BoundsMatrix makePairwiseBounds() const;

  /** @brief Generates a string graphviz representation of the modeled molecule
   *
   * The graph contains basic connectivity, stereopermutator information
   * and internal coordinate bounds.
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @return A string that can be converted into an image using graphviz of
   *   the molecule being modeled.
   */
  std::string dumpGraphviz() const;

  /** @brief Writes a graphviz representation of a modeled molecule to a file
   *
   * The graph contains basic connectivity, stereopermutator information
   * and internal coordinate bounds.
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @param filename The filename to which to write the graphviz representation
   */
  void writeGraphviz(const std::string& filename) const;
//!@}

private:
  // Molecule closure
  const Molecule& _molecule;

  //! Constraints by fixed positions
  BoundsMapType<2> _constraints;
  //! Bond distance value bounds on index pairs
  BoundsMapType<2> _bondBounds;
  //! Angle value bounds on index triples
  BoundsMapType<3> _angleBounds;
  //! Dihedral value bounds on index quadruplets
  BoundsMapType<4> _dihedralBounds;

  //! Chiral constraints
  std::vector<ChiralConstraint> _chiralConstraints;
  std::vector<DihedralConstraint> _dihedralConstraints;

//!@name Private methods
//!@{
  //! Adds [0, π] default angle bounds for all bonded atom triples
  void _addDefaultAngles();

  /*!
   * @brief Adds [0, π] default dihedrals to the model (of connected sequences)
   *
   * Adds [0, π] default dihedrals to the model. Avoids treating 1-4 pairs as
   * nonbonded.
   */
  void _addDefaultDihedrals();

  /*!
   * It is permissible in some circumstances not to have AtomStereopermutators
   * even on non-terminal atoms. These have to be re-added in order for us
   * to be able to model everywhere.
   */
  void _instantiateMissingAtomStereopermutators(random::Engine& engine);

  //! Add bond distances to the underlying model
  void _modelBondDistances(
    const FixedPositionsMapType& fixedAngstromPositions,
    double looseningFactor
  );

  //! Add precise information to the model for cycles that are flat
  void _modelFlatCycles(
    const FixedPositionsMapType& fixedAngstromPositions,
    double looseningFactor
  );

  //! Adds spirocenter modelling to the data set
  void _modelSpirocenters(
    const FixedPositionsMapType& fixedAngstromPositions
  );
//!@}
};

} // namespace distance_geometry
} // namespace molassembler
} // namespace Scine

#endif
