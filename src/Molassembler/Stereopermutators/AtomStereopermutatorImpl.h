/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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
/* Typedefs */
  using StereopermutationType = Stereopermutations::Stereopermutation;

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
    const Graph& graph,
    AtomIndex centerAtom,
    const Shapes::Shape shape,
    const RankingInformation& ranking,
    const TemperatureRegime temperature
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

/* Modification */
  //! Changes the assignment of the stereopermutator
  void assign(boost::optional<unsigned> assignment);

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
  void fit(
    const Graph& graph,
    const AngstromPositions& angstromWrapper
  );

  /*!
   * In case a graph modification changes the ranking of this stereopermutator's
   * substituents, it must be redetermined whether the new configuration is a
   * stereopermutator and if so, which assignment corresponds to the previous one.
   */
  boost::optional<PropagatedState> propagate(
    const Graph& graph,
    RankingInformation newRanking,
    boost::optional<Shapes::Shape> shapeOption
  );

  /*!
   * Prepares for the removal of an atom on the graph level, which involves
   * the generation of new atom indices.
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  //! If the shape is changed, we must adapt
  void setShape(
    Shapes::Shape shape,
    const Graph& graph
  );

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

  /*!
   * @brief Returns the underlying AbstractStereopermutation
   * @note This is library-internal and not part of the public API
   */
  const Stereopermutators::Abstract& getAbstract() const;

  /*!
   * @brief Returns the underlying FeasibleStereopermutation
   * @note This is library-internal and not part of the public API
   */
  const Stereopermutators::Feasible& getFeasible() const;

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

  //! Models abstract stereopermutations and decides three-dimensional feasibility
  Stereopermutators::Feasible feasible_;

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
