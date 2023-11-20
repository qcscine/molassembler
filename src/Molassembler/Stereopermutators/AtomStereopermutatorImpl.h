/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Private implementation of AtomStereopermutator
 */

#ifndef INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_IMPL_H
#define INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_IMPL_H

#include "Molassembler/AtomStereopermutator.h"

#include "Molassembler/DistanceGeometry/DistanceGeometry.h"
#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"

#include "boost/optional.hpp"

namespace Scine {
namespace Molassembler {

namespace Shapes {
namespace Properties {
struct ShapeTransitionGroup;
} // namespace Properties
} // namespace Shapes

class AtomStereopermutator::Impl : public Temple::Crtp::LexicographicComparable<Impl> {
public:
//!@name Static functions
//!@{
  /*!
   * @brief Picks a shape retaining as much chiral state as possible on a
   *   shape size increase
   * @throws std::logic_error If there are no larger shapes
   * @note Behavior is dependent on ChiralStatePreservation option
   */
  static Shapes::Shape up(Shapes::Shape shape);

  /*!
   * @brief Picks a shape retaining as much chiral state as possible on a
   *   shape size decrease
   * @throws std::logic_error If there are no smaller shapes
   * @note Behavior is dependent on ChiralStatePreservation option
   */
  static Shapes::Shape down(Shapes::Shape shape, Shapes::Vertex removedVertex);

  /*!
   * @brief Selects a shape vertex mapping from a shape transition group
   */
  static boost::optional<std::vector<Shapes::Vertex>> selectTransitionMapping(
    const Shapes::Properties::ShapeTransitionGroup& mappingsGroup,
    const ChiralStatePreservation& preservationOption
  );

  //! @brief Whether the stereopermutations interconvert rapidly at selected temp
  static bool thermalized(
    AtomIndex centerAtom,
    Shapes::Shape shape,
    const RankingInformation& ranking,
    const Graph& graph
  );
//!@}

/* Constructors */
  /**
   * @brief Constructor
   *
   * @param graph Base graph
   * @param shape Shape the stereopermutator represents
   * @param centerAtom Placement
   * @param ranking Ranking of its substituents
   */
  Impl(
    const Graph& graph,
    Shapes::Shape shape,
    AtomIndex centerAtom,
    RankingInformation ranking
  );

  Impl(
    AtomIndex centerAtom,
    Shapes::Shape shape,
    RankingInformation ranking,
    const FeasiblesGenerator& feasibility,
    const ThermalizationPredicate& thermalization,
    const std::vector<std::vector<SiteIndex>>& siteGroups = {});

/* Modification */
  //! Changes the assignment of the stereopermutator
  void assign(boost::optional<unsigned> assignment,
              const std::vector<std::vector<SiteIndex>>& siteGroups = {});

  /*! @brief Assign with a vertex mapping
   *
   * @pre Requires that this mapping was generated with centroid last!
   */
  void assign(std::vector<Shapes::Vertex> vertexMapping);

  //! Assigns the Stereopermutator randomly using relative stereopermutation weights
  void assignRandom(Random::Engine& engine);

  /**
   * @brief Applies an atom index permutation
   *
   * @param permutation The permutation to apply
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*!
   * The shape and assignment are determined based on three-dimensional
   * positions using angle and chiral distortions from the respective idealized
   * shapes.
   */
  boost::optional<ShapeMap> fit(
    const SiteCentroids& centroids,
    const FeasiblesGenerator& feasibility,
    const ThermalizationPredicate& thermalization
  );

  /*!
   * In case a graph modification changes the ranking of this stereopermutator's
   * substituents, it must be redetermined whether the new configuration is a
   * stereopermutator and if so, which assignment corresponds to the previous one.
   */
  boost::optional<PropagatedState> propagate(
    RankingInformation newRanking,
    boost::optional<Shapes::Shape> shapeOption,
    const FeasiblesGenerator& feasibility,
    const ThermalizationPredicate& thermalization
  );

  /*!
   * Prepares for the removal of an atom on the graph level, which involves
   * the generation of new atom indices.
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  //! If the shape is changed, we must adapt
  void setShape(
    Shapes::Shape shape,
    const FeasiblesGenerator& feasibility,
    const ThermalizationPredicate& thermalization,
    const std::vector<std::vector<SiteIndex>>& siteGroups = {}
  );

  //! Unconditionally alter the thermalization of stereopermutations
  void thermalize(bool thermalization);

/* Information */
  //! Returns the angle between two site indices in the idealized shape
  double angle(SiteIndex i, SiteIndex j) const;

  /*!
   * @brief Returns the permutation index within the set of possible
   *   permutations, if set
   *
   * Returns the (public) information of whether the stereopermutator is assigned
   * or not, and if so, which assignment it is.
   */
  boost::optional<unsigned> assigned() const;

  //! Returns a single-element vector containing the central atom
  AtomIndex placement() const;

  /*!
   * @brief Returns IOP within the set of symbolic ligand permutations
   *
   * This is different to the assignment. The assignment denotes the index
   * within the set of possible (more specifically, not obviously infeasible)
   * stereopermutations.
   */
  boost::optional<unsigned> indexOfPermutation() const;

  /*!
   * @brief Returns a minimal representation of chiral constraints
   *
   * Every minimal representation consists only of ligand indices.
   *
   * The minimal representation assumes that all shape tetrahedron
   * definitions are defined to be Positive targets, which is checked in
   * the shapes tests.
   *
   * @param enforce Emit minimal representations of chiral constraints even if
   * the stereopermutator does not have any chiral state, i.e.
   * numStereopermutators() <= 1, as long as it is assigned.
   */
  std::vector<MinimalChiralConstraint> minimalChiralConstraints(bool enforce = false) const;

  //! Returns an information string for diagnostic purposes
  std::string info() const;

  //! Returns an information string for ranking equality checking purposes
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

  inline bool thermalized() const {
    return thermalized_;
  }

  /*!
   * @brief Returns the underlying AbstractStereopermutation
   * @note This is library-internal and not part of the public API
   */
  const Stereopermutators::Abstract& getAbstract() const;

  //! Returns the list of feasible abstract permutations
  const std::vector<unsigned>& getFeasible() const;

  //! Returns the underlying ranking
  const RankingInformation& getRanking() const;

  //! Returns the underlying shape
  Shapes::Shape getShape() const;

  /*!
   * @brief Yields the mapping from site indices to shape positions
   * @throws std::logic_error if the stereopermutator is unassigned.
   */
  const ShapeMap& getShapePositionMap() const;

  /*!
   * @brief Returns the number of possible permutations
   *
   * Fetches the number of different assignments possible with the current
   * substituent ranking and connectivity information. This is also the upper
   * exclusive bound on the assignment indices that can be used to change the
   * arrangement of substituents.
   */
  unsigned numAssignments() const;

  /*!
   * @brief Returns the number of symbolic ligand permutations
   *
   * Fetches the number of permutations determined by symbolic ligand
   * calculation, not considering linking or haptic ligand cones.
   */
  unsigned numStereopermutations() const;

/* Operators */
  inline auto tie() const {
    return std::make_tuple(shape_, centerAtom_, numStereopermutations(), assigned());
  }

private:
/* State */
  //! Central atom of the Stereopermutator
  AtomIndex centerAtom_;

  //! The shape the stereopermutator represents
  Shapes::Shape shape_;

  //! Ranking information of substituents
  RankingInformation ranking_;

  //! Abstract stereopermutations and intermediate state
  Stereopermutators::Abstract abstract_;

  //! List of abstract permutation indices that are feasible in three dimensions
  std::vector<unsigned> feasibles_;

  //! The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> assignmentOption_;

  //! Derived property of @p assignmentOption_
  ShapeMap shapePositionMap_;

  //! Whether all feasible assignments interconvert thermally
  bool thermalized_;
};

} // namespace Molassembler
} // namespace Scine

#endif
