/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Handle rotational arrangements of adjacent atom-centered shapes
 *
 * Contains the BondStereopermutator class declaration, which models E/Z double bond
 * stereopermutators in molecules.
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_H
#define INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_H

#include "Molassembler/AtomStereopermutator.h"

#include <string>

namespace Scine {
namespace Molassembler {

/* Forward-declarations */
namespace Stereopermutations {
class Composite;
} // namespace Stereopermutations

namespace Random {
class Engine;
} // namespace Random

class AngstromPositions;
class PrivateGraph;
class StereopermutatorList;

namespace DistanceGeometry {
class SpatialModel;
struct ChiralConstraint;
} // namespace DistanceGeometry

/**
 * @brief Handles specific relative arrangements of two atom stereopermutators
 *   joined by a bond
 *
 * This class exists to model rotational barriers along bonds joining an
 * arbitrary pair of idealized shapes.
 */
class MASM_EXPORT BondStereopermutator {
public:
//!@name Public types
//!@{
  /**
   * @brief How dihedrals are aligned in the generation of stereopermutations
   */
  enum class Alignment {
    //! At least two shape vertices eclipse one another along the axis
    Eclipsed,
    //! At least one pair of substituents are staggered along the axis
    Staggered,
    //! Both eclipsed and staggered alignments
    EclipsedAndStaggered,
    //! Offset exactly halfway between eclipsed and staggered alignments
    BetweenEclipsedAndStaggered
  };

  /**
   * @brief Differentiates how viable assignments are chosen during fitting.
   */
  enum class FittingMode {
    //! Positions must be close to the idealized assignment geometry
    Thresholded,
    //! Choose whichever assignment best represents the geometry directly
    Nearest
  };

  //! Spatial centroids of atom stereopermutator site atoms
  using SitePositions = Eigen::Matrix<double, 3, Eigen::Dynamic>;
  //! Spatial centroids of two atom stereopermutators' site atoms
  using SitePositionsPair = std::pair<SitePositions, SitePositions>;
//!@}

//!@name Static members
//!@{
  /*! @brief Assignment acceptance is dependent on this parameter
   *
   * An assignment is accepted if the sum of absolute angular deviations of all
   * dihedrals involved in a stereopermutation divided by the number of
   * dihedrals is less than this factor times two pi divided by the larger
   * number of involved substituents at the ends of the bond.
   *
   * For the classic E/Z double bond case, there are two involved substituents
   * at both sides. The tolerance per dihedral is then 0.1 * 360° / 2 = 18°.
   */
  static constexpr double assignmentAcceptanceParameter = 0.1;
//!@}

//!@name Static functions
//!@{
  /*! @brief Determine which stereopermutations aren't obviously infeasible
   *
   * @complexity{@math{\Theta(P)}, where @math{P} is the number of permutations
   * in the underlying composite}
   */
  static std::vector<unsigned> notObviouslyInfeasibleStereopermutations(
    const PrivateGraph& graph,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    const Stereopermutations::Composite& composite
  );
//!@}

//!@name Special member functions
//!@{
  BondStereopermutator(BondStereopermutator&& other) noexcept;
  BondStereopermutator& operator = (BondStereopermutator&& other) noexcept;
  BondStereopermutator(const BondStereopermutator& other);
  BondStereopermutator& operator = (const BondStereopermutator& other);
  ~BondStereopermutator();
//!@}

//!@name Constructors
//!@{
  BondStereopermutator() = delete;
  /*!
   * @brief Constructs a bond stereopermutator on two atom stereopermutators
   *   without checking whether stereopermutations are feasible or not
   *
   * @complexity{@math{O(S!)} where @math{S} is the size of the larger involved
   * shape}
   */
  BondStereopermutator(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    const BondIndex& edge,
    Alignment alignment = Alignment::Eclipsed
  );

  /*! @brief Constructs a bond stereopermutator on two atom stereopermutators,
   *   removing obviously infeasible stereopermutations
   *
   * @complexity{@math{O(S!)} where @math{S} is the size of the larger involved
   * shape}
   */
  BondStereopermutator(
    const PrivateGraph& graph,
    const StereopermutatorList& stereopermutators,
    const BondIndex& edge,
    Alignment alignment = Alignment::Eclipsed
  );
//!@}

//!@name Modification
//!@{
  /*!
   * @brief Changes the assignment of the stereopermutator
   *
   * @param assignment The new assignment of the stereopermutator. May be
   *   @p boost::none, which sets the chiral state as indeterminate. Must be
   *   less than the number of assignments if not None.
   *
   * @complexity{@math{\Theta(1)}}
   */
  void assign(boost::optional<unsigned> assignment);

  /*!
   * @brief Assign the Stereopermutator at random
   *
   * @complexity{@math{\Theta(1)}}
   * @parblock @note If the stereocenter is already assigned, it is reassigned.
   * @endparblock
   * @parblock @note The state of the passed PRNG is advanced.
   * @endparblock
   */
  void assignRandom(Random::Engine& engine);

  /** @brief Applies an atom index permutation
   *
   * @complexity{@math{\Theta(1)}}
   * @param permutation The permutation to apply
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  struct FittingReferences {
    const AtomStereopermutator& stereopermutator;
    const AtomStereopermutator::ShapeMap& shapeMap;
  };

  /*! @brief Determines the assignment the permutator is in from positional
   *   information
   *
   * The assignment of this permutator is determined from three-dimensional
   * positions using penalties to modeled dihedral angles. An assignment is
   * considered matched if each dihedral that is part of the permutation is
   * within @p assignmentAcceptanceDihedralThreshold.
   *
   * @complexity{@math{\Theta(P)} where @math{P} is the number of permutations}
   *
   * @param sitePositions The site positions of each constuting atom
   *   stereopermutator, following the same ordering as compositeAlignment()
   * @param fittingReferences Pair of FittingReferences to constituting atom
   *   stereopermutators (order is irrelevant)
   * @param mode Mode altering the assignment of stereopermutations depending
   *   on geometric closeness to the idealized minimum
   */
  void fit(
    const SitePositionsPair& sitePositions,
    std::pair<FittingReferences, FittingReferences> fittingReferences,
    FittingMode mode = FittingMode::Thresholded
  );

  //! @overload
  void fit(
    const AngstromPositions& angstromWrapper,
    std::pair<FittingReferences, FittingReferences> fittingReferences,
    FittingMode mode = FittingMode::Thresholded
  );

  /*!
   * @brief Propagates the bond stereocenter's state through a possible ranking
   *   change on one of its constituting atom stereopermutators.
   *
   * @param oldPermutator The AtomStereopermutator that is changing at its
   *   state before a ranking change effects its current stereopermutation or
   *   substituent ranking
   * @param newPermutator The AtomStereopermutator that is changing after its
   *   internal state has been propagated through a ranking change
   * @param inner Inner graph for context
   * @param permutators Stereopermutator list for context
   *
   * @complexity{@math{O(S!)} where @math{S} is the size of the larger involved
   * shape}
   */
  void propagateGraphChange(
    const AtomStereopermutator::PropagatedState& oldPermutator,
    const AtomStereopermutator& newPermutator,
    const PrivateGraph& inner,
    const StereopermutatorList& permutators
  );

  //! Propagates invalidated atom indices through the removal of an index
  void propagateVertexRemoval(AtomIndex removedIndex);
//!@}

//!@name Information
//!@{
  /*! @brief Returns alignment parameter this was constructed with
   *
   * @complexity{@math{\Theta(1)}}
   */
  Alignment alignment() const;

  /*! @brief Returns the permutation index within the set of possible permutations, if set
   *
   * @complexity{@math{\Theta(1)}}
   * @returns whether the stereopermutator is assigned or not, and if so,
   * which assignment it is.
   */
  boost::optional<unsigned> assigned() const;

  /*! @brief Gives read-only access to the underlying Composite object
   *
   * @complexity{@math{\Theta(1)}}
   */
  const Stereopermutations::Composite& composite() const;

  //! The atom identifier alignment of the permutational composite
  std::pair<AtomIndex, AtomIndex> compositeAlignment() const;

  /*! @brief Angle between sites at stereopermutators in the current assignment
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param stereopermutatorA One constituting stereopermutator
   * @param siteIndexA A site index of @p stereopermutatorA
   * @param stereopermutatorB The other constituting stereopermutator
   * @param siteIndexB A site index of @p stereopermutatorB
   *
   * @throws std::logic_error If the stereopermutator is unassigned or if no
   *   dihedral is found for the passed sites
   */
  double dihedral(
    const AtomStereopermutator& stereopermutatorA,
    SiteIndex siteIndexA,
    const AtomStereopermutator& stereopermutatorB,
    SiteIndex siteIndexB
  ) const;

  /*! @brief Returns whether this stereopermutator has the same relative
   *   orientation as another stereopermutator
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool hasSameCompositeOrientation(const BondStereopermutator& other) const;

  /*! @brief Returns the index of permutation
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<unsigned> indexOfPermutation() const;

  /*! @brief Returns the number of possible assignments
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned numAssignments() const;

  /*! @brief Returns the number of possible stereopermutations
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned numStereopermutations() const;

  /*! @brief Returns an information string for diagnostic purposes
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::string info() const;

  /*! @brief Returns an information for ranking equality checking purposes
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::string rankInfo() const;

  /*! @brief Returns which bond this stereopermutator is placed on in the molecule
   *
   * @complexity{@math{\Theta(1)}}
   */
  BondIndex placement() const;
//!@}

//!@name Operators
//!@{
  bool operator < (const BondStereopermutator& other) const;
  /*!
   * @brief Compares whether the underlying composite and the assignment are
   *   equivalent to those of another bond stereopermutator
   */
  bool operator == (const BondStereopermutator& other) const;
  //! Inverts @p operator==
  bool operator != (const BondStereopermutator& other) const;
//!@}

private:
  struct Impl;
  std::unique_ptr<Impl> pImpl_;
};

} // namespace Molassembler
} // namespace Scine

#endif
