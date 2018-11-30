// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_BOND_stereopermutator_H
#define INCLUDE_MOLASSEMBLER_BOND_stereopermutator_H

#include "boost/optional/optional_fwd.hpp"

#include "molassembler/Types.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

#include <cmath>
#include <vector>
#include <string>
#include <memory>

/*! @file
 *
 * @brief Handle rotational arrangements of adjacent atom-centered symmetries
 *
 * Contains the BondStereopermutator class declaration, which models E/Z double bond
 * stereopermutators in molecules.
 */

namespace stereopermutation {

class Composite;

} // namespace stereopermutation

namespace molassembler {

class AngstromWrapper;
class OuterGraph;
class RankingInformation;

// Forward-declarations
class AtomStereopermutator;

namespace DistanceGeometry {
class SpatialModel;
struct ChiralityConstraint;
} // namespace DistanceGeometry

/**
 * @brief Handles specific relative arrangements of two atom stereopermutators
 *   joined by a bond
 *
 * This class exists to model rotational barriers in bond orders higher than
 * Single that join an arbitrary pair of idealized symmetries.
 */
class BondStereopermutator {
public:
  //! The volume tolerance of emitted chirality constraints
  static constexpr double chiralityConstraintTolerance = 0.1;
  //! An Assignment is accepted if the fit for each dihedral is below this value
  static constexpr double assignmentAcceptanceDihedralThreshold = M_PI / 36.0; // 5°

//!@name Special member functions
//!@{
  BondStereopermutator(BondStereopermutator&& other) noexcept;
  BondStereopermutator& operator = (BondStereopermutator&& other) noexcept;
  BondStereopermutator(const BondStereopermutator& other);
  BondStereopermutator& operator = (const BondStereopermutator& other);
  ~BondStereopermutator();

  BondStereopermutator() = delete;
  //! Constructs a bond stereopermutator on two atom stereopermutators
  BondStereopermutator(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    const BondIndex& edge
  );
//!@}

//!@name Modifiers
//!@{
  /*!
   * @brief Changes the assignment of the stereopermutator
   *
   * @param assignment The new assignment of the stereopermutator. May be
   *   @p boost::none, which sets the chiral state as indeterminate. Must be
   *   less than the number of assignments if not None.
   */
  void assign(boost::optional<unsigned> assignment);

  /*!
   * @brief Assign the Stereopermutator at random
   *
   * @note If the stereocenter is already assigned, it is reassigned.
   */
  void assignRandom();

  /*!
   * @brief Determines the assignment the permutator is in from positional
   *   information
   *
   * The assignment of this permutator is determined from three-dimensional
   * positions using penalties to modeled dihedral angles. An assignment is
   * considered matched if each dihedral that is part of the permutation is
   * within @p assignmentAcceptanceDihedralThreshold.
   *
   * @param angstromWrapper The positional information to extract the assignment
   *   from
   * @param stereopermutatorA One constituting atom stereopermutator
   * @param stereopermutatorB The other constituting atom stereopermutator
   */
  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB
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
   */
  void propagateGraphChange(
    const AtomStereopermutator& oldPermutator,
    const AtomStereopermutator& newPermutator
  );
//!@}

//!@name Information
//!@{
  /*!
   * @brief Returns the permutation index within the set of possible permutations, if set
   * @returns whether the stereopermutator is assigned or not, and if so,
   * which assignment it is.
   */
  boost::optional<unsigned> assigned() const;

  //! Gives read-only access to the underlying Composite object
  const stereopermutation::Composite& composite() const;

  /*!
   * @brief Returns whether this stereopermutator has the same relative
   *   orientation as another stereopermutator
   */
  bool hasSameCompositeOrientation(const BondStereopermutator& other) const;

  /*!
   * @brief Returns the index of permutation
   * @note For BondStereopermutators, the index of permutation and the
   *   assignment index are the same as there are no stereopermutation exclusion
   *   criteria yet.
   */
  boost::optional<unsigned> indexOfPermutation() const;

  /*!
   * @brief Returns the number of possible assignments
   * @note For BondStereopermutators, the number of stereopermutations and the
   *   number of assignments are the same as there are no stereopermutation
   *   exclusion criteria yet.
   */
  unsigned numAssignments() const;

  /*!
   * @brief Returns the number of possible stereopermutations
   * @note For BondStereopermutators, the number of stereopermutations and the
   *   number of assignments are the same as there are no stereopermutation
   *   exclusion criteria yet.
   */
  unsigned numStereopermutations() const;

  //! Returns an information string for diagnostic purposes
  std::string info() const;

  //! Returns an information for ranking equality checking purposes
  std::string rankInfo() const;

  //! Returns which bond this stereopermutator is placed on in the molecule
  BondIndex edge() const;
//!@}

//!@name Operators
//!@{
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
