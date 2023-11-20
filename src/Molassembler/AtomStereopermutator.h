/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Handle arrangements of substituents at corners of an atom-centered
 *   shape
 *
 * Coordinative stereopermutator class header file. Permits the storage of
 * particular arrangements of bonded atoms around a central atom and their
 * manipulation.
 *
 * Handles the stereopermutation issue, allowing users to cycle through
 * non-mutually-superimposable arrangements of substituents, here called
 * 'assignments'.
 */

#ifndef INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_H
#define INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_H

#include "Molassembler/RankingInformation.h"
#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Temple/StrongIndexPermutation.h"

#include <array>
#include <memory>

namespace Scine {
namespace Molassembler {

/* Forward declarations */
struct RankingInformation;
class AngstromPositions;
class Graph;

namespace Stereopermutators {
struct Abstract;
struct Feasible;
} // namespace Stereopermutators

namespace Random {
class Engine;
} // namespace Random

namespace DistanceGeometry {

class SpatialModel;
struct ChiralConstraint;

} // namespace DistanceGeometry

/**
 * @brief Handles the steric permutation of substituents of a non-terminal
 *   central atom
 *
 * This class handles the permutation of ranked ligands around a central atom.
 * It models its haptic ligands' binding sites and bridges in multidentate
 * ligands in order to decide which stereopermutations are feasible. A
 * stereopermutation may be infeasible, i.e. not realizable in
 * three-dimensional space, if either haptic ligands would intersect due to
 * too close ligand angles for their spatial heft, or if a multidentate
 * ligand's bridge length between binding sites were too short to match the
 * angle. The list of stereopermutations reduced by infeasible
 * stereopermutations is then re-indexed and those indices into the list are
 * called assignments.
 *
 * A Stereopermutator can be unassigned, i.e. the distinct stereopermutation
 * that the substituents are can be indeterminate. If you choose to generate
 * conformers for a molecule that includes unassigned stereopermutators, every
 * conformer will choose an assignment from the pool of feasible assignments
 * randomly, but consistent with relative statistical occurrence weights.
 *
 * E.g. a square planar AABC ligand set will have an A-A cis stereopermutation
 * that occurs twice as often as the A-A trans stereopermutation.
 *
 * @note An instance of this class on a given central atom does not indicate
 *   that that atom is a stereocenter. That is only the case if there are
 *   multiple stereopermutations of the ranked substituents / ligands.
 */
class MASM_EXPORT AtomStereopermutator {
public:
  //! Mapping between site indices of a ranking and vertices of a shape
  using ShapeMap = Temple::StrongIndexPermutation<SiteIndex, Shapes::Vertex>;

  //! Spatial positions of the centroids of sites
  using SiteCentroids = Eigen::Matrix<double, 3, Eigen::Dynamic>;

  //! Predicate to decide whether an atom stereopermutator is thermalized
  using ThermalizationPredicate = std::function<bool(AtomIndex, Shapes::Shape, const RankingInformation&)>;

  //! Generator function to list spatially feasible permutation indices
  using FeasiblesGenerator = std::function<
    std::vector<unsigned>(
      const Stereopermutators::Abstract& abstract,
      Shapes::Shape shape,
      AtomIndex placement,
      const RankingInformation& ranking,
      std::vector<std::vector<SiteIndex>> siteGroups
    )
  >;

  //! Old state dumped upon propagation
  using PropagatedState = std::tuple<RankingInformation, ShapeMap>;

  /*!
   * @brief Site index sequence defining a chiral constraint. If a site index
   *   is None, then it denotes the position of the central index
   */
  using MinimalChiralConstraint = std::array<boost::optional<SiteIndex>, 4>;

//!@name Special member functions
//!@{
  AtomStereopermutator(AtomStereopermutator&& other) noexcept;
  AtomStereopermutator& operator = (AtomStereopermutator&& other) noexcept;
  AtomStereopermutator(const AtomStereopermutator& other);
  AtomStereopermutator& operator = (const AtomStereopermutator& other);
  ~AtomStereopermutator();

  /*! @brief Construct an AtomStereopermutator
   *
   * @param graph The molecule's graph. This information is needed to model
   *   haptic ligands.
   * @param shape The local idealized shape to model. Typically the
   *   result of Molecule's inferShape.
   * @param centerAtom The atom index within the molecule that is the center of
   *   the local idealized shape
   * @param ranking The ranking of the central atom's substituents and ligand
   *   sites. Typically the result of Molecule's rankPriority.
   *
   * @complexity{@math{L\cdot S!} where @math{L} is the number of links and
   * @math{S} is the size of @p shape}
   *
   * @post The stereopermutator is unassigned
   */
  [[deprecated("Prefer constructor without graph argument")]]
  AtomStereopermutator(
    const Graph& graph,
    Shapes::Shape shape,
    AtomIndex centerAtom,
    RankingInformation ranking
  );

  /*! @brief Construct an AtomStereopermutator
   *
   * @param centerAtom The atom index within the molecule that is the center of
   *   the local idealized shape
   * @param shape The local idealized shape to model. Typically the
   *   result of Molecule's inferShape.
   * @param ranking The ranking of the central atom's substituents and ligand
   *   sites. Typically the result of Molecule's rankPriority.
   * @param feasibility A functor for determining which abstract permutations
   *   are feasible.
   * @param thermalization A functor for deciding whether the stereopermutator
   *   is thermalized and all stereopermutations coalesce.
   * @param siteGroups The position groups for the binding sites.
   *
   * @complexity{@math{L\cdot S!} where @math{L} is the number of links and
   * @math{S} is the size of @p shape}
   *
   * @post The stereopermutator is unassigned
   */
  AtomStereopermutator(
    AtomIndex centerAtom,
    Shapes::Shape shape,
    RankingInformation ranking,
    const FeasiblesGenerator& feasibility = {},
    const ThermalizationPredicate& thermalization = {},
    const std::vector<std::vector<SiteIndex>>& siteGroups = {}
  );
//!@}

//!@name Shape picking
//!@{
  /*! @brief Picks a shape retaining as much chiral state as possible on a
   *   shape size increase
   *
   * @complexity{@math{O(S!)} if uncached, @math{\Theta(1)} otherwise}
   * @throws std::logic_error If there are no larger shapes
   */
  static Shapes::Shape up(Shapes::Shape shape);

  /*!
   * @brief Picks a shape retaining as much chiral state as possible on a
   *   shape size decrease
   *
   * @complexity{@math{O(S!)} if uncached, @math{\Theta(1)} otherwise}
   * @throws std::logic_error If there are no smaller shapes
   */
  static Shapes::Shape down(Shapes::Shape shape, Shapes::Vertex removedVertex);
//!@}

//!@name Statics
//!@{
  /*! @brief Decide whether stereopermutations are thermalized based on the graph
   *
   * If Options::Thermalization::pyramidalInversion is set, considers
   * asymmetric pyramidal nitrogen atom centers a stereocenter only if part of a
   * cycle of size 4 or smaller.
   *
   * If Options::Thermalization::berryPseudorotation is set, the shape is a
   * trigonal bipyramid and there are no links between sites,
   * stereopermutations are thermalized.
   *
   * If Options::Thermalization::bartellMechanism is set, the shape is a
   * pentagonal bipyramid and there are no links between sites,
   * stereopermutations are thermalized.
   */
  static bool thermalized(
    AtomIndex centerAtom,
    Shapes::Shape shape,
    const RankingInformation& ranking,
    const Graph& graph
  );

  //! Binds a graph instance into a functor invoking thermalized
  static auto thermalizationFunctor(const Graph& g) {
    return [&g](AtomIndex centerAtom, Shapes::Shape shape, const RankingInformation& ranking) {
      return AtomStereopermutator::thermalized(centerAtom, shape, ranking, g);
    };
  }
//!@}

//!@name Modifiers
//!@{
  /*! @brief Changes the assignment of the stereopermutator
   *
   * @param assignment The new assignment of the stereopermutator. May be
   *   @p boost::none, which sets the chiral state as indeterminate. Must be
   *   less than the number of assignments if not None.
   * @param siteGroups The shape position groups of the sites.
   *
   * @complexity{@math{\Theta(1)} if @p assignment is @c boost::none.
   * @math{\Theta(S)} otherwise}
   */
  void assign(boost::optional<unsigned> assignment,
              const std::vector<std::vector<SiteIndex>>& siteGroups = {});

  /*! @brief Assign the Stereopermutator randomly using relative statistical weights
   *
   * Stereopermutations are generated with relative statistical occurrence
   * weights. The assignment is then chosen from the possible stereopermutations
   * with a discrete distribution whose weights are the corresponding relative
   * statistical occurrences.
   *
   * @complexity{@math{\Theta(1)} if @p assignment is @c boost::none.
   * @math{\Theta(S)} otherwise}
   *
   * @parblock @note If the stereocenter is already assigned, it is reassigned.
   * @endparblock
   *
   * @parblock @note The state of the passed PRNG is advanced.
   * @endparblock
   */
  void assignRandom(Random::Engine& engine);

  /** @brief Applies an atom index permutation
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param permutation The permutation to apply
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*! @brief Determine the shape and assignment realized in positions
   *
   * The shape and assignment are determined from Cartesian coordinates by the
   * probability that the continuous shape measure with respect to an idealized
   * shape is least likely to be drawn from a sample of randomly generated
   * vertex positions. For shapes with many vertices, generating enough samples
   * for a statistical argument is too expensive, and the shape whose
   * continuous shape measure is lowest is chosen instead.
   *
   * @param centroids Spatial centroids of each site in the contained ranking,
   *   and the position of the central index last
   * @returns A mapping of site indices to shape vertices if a
   *   stereopermutation could be found, None otherwise
   *
   * @complexity{@math{\Theta(S!)}}
   */
  boost::optional<ShapeMap> fit(
    const SiteCentroids& centroids,
    const FeasiblesGenerator& feasibility = {},
    const ThermalizationPredicate& thermalization = {}
  );

  /*! @brief Determine the shape and assignment realized in positions
   *
   * The shape and assignment are determined from Cartesian coordinates by the
   * probability that the continuous shape measure with respect to an idealized
   * shape is least likely to be drawn from a sample of randomly generated
   * vertex positions. For shapes with many vertices, generating enough samples
   * for a statistical argument is too expensive, and the shape whose
   * continuous shape measure is lowest is chosen instead.
   *
   * @param angstromWrapper The wrapped positions
   * @returns A mapping of site indices to shape vertices if a
   *   stereopermutation could be found, None otherwise
   *
   * @complexity{@math{\Theta(S!)}}
   */
  boost::optional<ShapeMap> fit(
    const AngstromPositions& wrapper,
    const FeasiblesGenerator& feasibility = {},
    const ThermalizationPredicate& thermalization = {}
  );

  [[deprecated("Prefer graph-less alternative parameters")]]
  boost::optional<ShapeMap> fit(
    const Graph& graph,
    const AngstromPositions& angstromWrapper
  );

  /*! @brief Propagate the stereocenter state through a possible ranking change
   *
   * In case a graph modification changes the ranking of this stereopermutator's
   * substituents, it must be redetermined whether the new configuration is a
   * stereopermutator and if so, which assignment corresponds to the previous one.
   *
   * @complexity{@math{L\cdot S!} where @math{L} is the number of links and
   * @math{S} is the size of @p shape}
   */
  MASM_NO_EXPORT boost::optional<PropagatedState> propagate(
    RankingInformation newRanking,
    boost::optional<Shapes::Shape> shapeOption,
    const FeasiblesGenerator& feasibility = {},
    const ThermalizationPredicate& thermalization = {}
  );

  /*! @brief Adapts atom indices in the internal state to the removal of an atom
   *
   * Atom indices are adapted to a graph-level removal of an atom. The removed
   * index is changed to a placeholder value.
   *
   * @complexity{@math{\Theta(1)}}
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  /*! @brief Change the underlying shape of the permutator
   *
   * @complexity{@math{L\cdot S!} where @math{L} is the number of links and
   * @math{S} is the size of @p shape}
   *
   * @todo Consider trying to propagate within same shape size
   * @post The permutator is unassigned (chiral state is discarded)
   */
  void setShape(
    Shapes::Shape shape,
    const FeasiblesGenerator& feasibility = {},
    const ThermalizationPredicate& thermalization = {},
    const std::vector<std::vector<SiteIndex>>& siteGroups = {}
  );

  [[deprecated("Prefer graph-less alternative parameters")]]
  void setShape(
    Shapes::Shape shape,
    const Graph& graph,
    const std::vector<std::vector<SiteIndex>>& siteGroups = {}
  );

  //! Unconditionally thermalize the stereopermutations
  void thermalize(bool thermalization = true);
//!@}

//!@name Observers
//!@{
  /*! @brief Fetches angle between binding sites in the idealized shape
   *
   * @pre The stereopermutator must be assigned
   *
   * @param i Site index one
   * @param j Site index two
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @throws std::runtime_error If the stereopermutator is unassigned
   * @throws std::out_of_range If any site index is out of range
   *
   * @sa getRanking()
   */
  double angle(SiteIndex i, SiteIndex j) const;

  /*! @brief Returns the permutation index within the set of feasible permutations, if set
   *
   * Returns the information of whether the stereopermutator is assigned
   * or not, and if so, which assignment it is.
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<unsigned> assigned() const;

  /*! @brief Returns the atom this permutator is placed on
   *
   * @complexity{@math{\Theta(1)}}
   */
  AtomIndex placement() const;

  /*! @brief Returns IOP within the set of symbolic ligand permutations
   *
   * This is different to the assignment. The assignment denotes the index
   * within the set of possible (more specifically, not obviously infeasible)
   * stereopermutations.
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<unsigned> indexOfPermutation() const;

  /*! @brief Returns a minimal representation of chiral constraints
   *
   * Every minimal representation consists only of site indices. If no site
   * index is present, this position is the location of the central atom.
   *
   * The minimal representation assumes that all shape tetrahedra are
   * defined as Positive targets, which is checked in the
   * shapes tests.
   *
   * @complexity{@math{\Theta(T)} where @math{T} is the number of tetrahedra
   * defined for the modeled shape}
   *
   * @param enforce Emit minimal representations of chiral constraints even if
   * the stereopermutator does not have any chiral state, i.e.
   * numStereopermutators() <= 1, as long as it is assigned.
   */
  std::vector<MinimalChiralConstraint> minimalChiralConstraints(bool enforce = false) const;

  /*! @brief Returns an information string for diagnostic purposes
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::string info() const;

  /*! @brief Returns an information string for ranking equality checking purposes
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::string rankInfo() const;

  /*! @brief Returns site indices grouped by rotational interconversion
   *
   * @pre The stereopermutator must be assigned
   *
   * @complexity{@math{\Theta(S^2)}}
   *
   * @throws std::runtime_error If the stereopermutator is unassigned
   */
  std::vector<std::vector<SiteIndex>> siteGroups() const;

  //! Generate site centroid positions from a whole-molecule set of positions
  SiteCentroids sitePositions(
    const AngstromPositions& wrapper,
    const std::vector<std::pair<AtomIndex, AtomIndex>>& substitutions = {}
  ) const;

  //! Returns whether the stereopermutations are thermalized
  bool thermalized() const;

  /*! @brief Returns the underlying feasible stereopermutations object
   *
   * @complexity{@math{\Theta(1)}}
   * @note This is library-internal and not part of the public API
   */
  MASM_NO_EXPORT const Stereopermutators::Abstract& getAbstract() const;

  /*!  @brief Returns the feasible stereopermutation indices
   *
   * @complexity{@math{\Theta(1)}}
   * @note This is library-internal and not part of the public API
   */
  MASM_NO_EXPORT const std::vector<unsigned>& getFeasible() const;

  /*! @brief Returns the underlying ranking
   *
   * @complexity{@math{\Theta(1)}}
   */
  const RankingInformation& getRanking() const;

  /*! @brief Returns the underlying shape
   *
   * @complexity{@math{\Theta(1)}}
   */
  Shapes::Shape getShape() const;

  /*! @brief Yields the mapping from site indices to shape positions
   *
   * @complexity{@math{\Theta(1)}}
   * @throws std::logic_error if the stereopermutator is unassigned.
   */
  const ShapeMap& getShapePositionMap() const;

  /*! @brief Returns the number of possible assignments
   *
   * The number of possible assignments is the number of non-superposable
   * arrangements of the abstract ligand case reduced by trans-arranged
   * multidentate pairs where the bridge length is too short or overlapping
   * haptic cones.
   *
   * For instance, if octahedral M[(A-A)3], there are four abstract arrangements
   * - trans-trans-trans
   * - trans-cis-cis
   * - 2x cis-cis-cis (Δ and Λ isomers, ship propeller-like)
   *
   * However, the number of stereopermutations for a concrete case in which the
   * bridges are too short to allow trans bonding is reduced by all
   * arrangements containing a trans-bonded bidentate ligand, i.e. only Δ and Λ
   * remain. The number of assignments is then only two.
   *
   * This is the upper exclusive bound on Some-type arguments to assign().
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned numAssignments() const;

  /*!
   * @brief Returns the number of possible stereopermutations
   *
   * The number of possible stereopermutations is the number of
   * non-superposable arrangements of the abstract ligand case without removing
   * trans-arranged multidentate pairs or overlapping haptic cones.
   *
   * For instance, if octahedral M[(A-A)3], there are four abstract arrangements
   * - trans-trans-trans
   * - trans-cis-cis
   * - 2x cis-cis-cis (Δ and Λ isomers, ship propeller-like chirality)
   *
   * However, the number of assignments for a concrete case in which the bridges
   * are too short to allow trans binding is reduced by all arrangements
   * containing a trans-bonded bidentate ligand, i.e. only Δ and Λ remain.
   *
   * Fetches the number of permutations determined by symbolic ligand
   * calculation, not considering linking or haptic ligand cones.
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned numStereopermutations() const;
//!@}


//!@name Operators
//!@{
  /** @brief Checks whether the underlying shape, central atom index, number
   *   of stereopermutations and current assignment match
   *
   * @param other The other atom stereopermutator to compare against
   *
   * @return Whether the underlying shape, central atom index, number of
   *   stereopermutations and current assignment match
   */
  bool operator == (const AtomStereopermutator& other) const;
  //! Inverts operator ==
  bool operator != (const AtomStereopermutator& other) const;

  /**
   * @brief Lexicographically compares the central atom index, the shape,
   *   the number of stereopermutations, and the current assignment
   *
   * @param other The other atom stereopermutator to compare against
   *
   * @return Which atom stereopermutator is lexicographically smaller
   */
  bool operator < (const AtomStereopermutator& other) const;
//!@}

private:
  class Impl;
  std::unique_ptr<Impl> pImpl_;
};

} // namespace Molassembler
} // namespace Scine

#endif
