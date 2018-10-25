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

namespace molassembler {

class AngstromWrapper;

// Forward-declarations
class AtomStereopermutator;

namespace DistanceGeometry {
class SpatialModel;
struct ChiralityConstraint;
} // namespace DistanceGeometry

class BondStereopermutator {
public:
  //! The volume tolerance of emitted chirality constraints
  static constexpr double chiralityConstraintTolerance = 0.1;
  //! An Assignment is accepted if the fit for each dihedral is below this value
  static constexpr double assignmentAcceptanceDihedralThreshold = M_PI / 60.0; // ~3°

//!@name Special member functions
//!@{
  BondStereopermutator(BondStereopermutator&& other) noexcept;
  BondStereopermutator& operator = (BondStereopermutator&& other) noexcept;
  BondStereopermutator(const BondStereopermutator& other);
  BondStereopermutator& operator = (const BondStereopermutator& other);
  ~BondStereopermutator();

  BondStereopermutator() = delete;
  /*!
   * @warning No haptic ligands allowed.
   * @note The left stereopermutator must have the smaller
   */
  BondStereopermutator(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    const BondIndex& edge
  );
//!@}

//!@name Modifiers
//!@{
  void assign(boost::optional<unsigned> assignment);

  void assignRandom();

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB
  );
//!@}

//!@name Information
//!@{
  boost::optional<unsigned> assigned() const;

  bool hasSameCompositeOrientation(const BondStereopermutator& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB
  ) const;

  std::string info() const;

  std::string rankInfo() const;

  BondIndex edge() const;

  void setModelInformation(
    DistanceGeometry::SpatialModel& model,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    double looseningMultiplier
  ) const;
//!@}

//!@name Operators
//!@{
  bool operator == (const BondStereopermutator& other) const;
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
