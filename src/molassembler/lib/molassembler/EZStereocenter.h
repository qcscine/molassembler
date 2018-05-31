#ifndef INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H
#define INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H

#include "Stereocenter.h"
#include "AngstromWrapper.h"

/*! @file
 *
 * Contains the EZStereocenter class declaration, which models E/Z double bond
 * stereocenters in molecules.
 */

/* TODO
 * - Could make dihedralAngleVariance constexpr, but unsure about effects on
 *   static initialization and if extra declaration required
 *   Also, consider if it violates single responsibility principle
 *
 * - Deprecate the entire class and polymorphic handling of both stereocenters
 *   throughout the library in favor of an overarching construct
 */

namespace molassembler {

namespace Stereocenters {

/*! Models stereogenic E/Z conformers in molecular graphs.
 *
 * Has either one or two assignments:
 *
 * - When on either side both substituents have equal rank, there is only one
 *   assignment.
 * - Otherwise, there are two assignments.
 * - In case there is only one substituent on one side, it is considered
 *   "unequal to nothingness". This is so that e.g. diazene::
 *
 *     H
 *       \
 *        N = N
 *              \
 *                H
 *
 *   is properly considered as having E and Z conformers
 *
 * Implementation-wise, _isZOption models the permutations as:
 * - Some true -> Z
 * - Some false -> E,
 * - None -> unassigned
 *
 * EZStereocenter is explicitly also instantiable on bonds that are not
 * stereogenic so that in DG, chirality constraints and dihedral limits that
 * result from the double bond are still collectible.
 */
class EZStereocenter : public Stereocenter {
public:
  using IndexPairsSet = std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >;

  static constexpr double chiralityConstraintTolerance = 1e-2;

private:
  struct DihedralLimits {
    const std::array<AtomIndexType, 4> indices;
    const double lower, upper;

    DihedralLimits(
      const std::array<AtomIndexType, 4>& indices,
      const std::pair<double, double>& limits
    ) : indices(indices),
        lower(limits.first),
        upper(limits.second)
    {
      assert(lower < upper);
    }
  };

/* State */
  //! Stores the indices of the minimal indices required to have a stereocenter
  AtomIndexType _leftCenter, _rightCenter;

  //! Store the ranking of indices
  RankingInformation _leftRanking, _rightRanking;

  //! Stores the current assignment state of the stereocenter
  boost::optional<bool> _isZOption;

/* Private members */
  static inline AtomIndexType highPriority(const RankingInformation& ranking) {
    /* NOTE: high priority is assigned from the back of the ranking since the list
     * is ordered ascending
     *
     * NOTE: back-back always works out to the right selection for all valid
     * number of ligands.
     */
    return ranking.ligands.at(
      // Highest ranked ligand index
      ranking.ligandsRanking.back().back()
    ).back();
  }

  static inline AtomIndexType lowPriority(const RankingInformation& ranking) {
    return ranking.ligands.at(
      // Lowest ranked ligand index
      ranking.ligandsRanking.front().front()
    ).front();
  }

  inline AtomIndexType _leftHighPriority() const {
    return highPriority(_leftRanking);
  }

  inline AtomIndexType _leftLowPriority() const {
    return lowPriority(_leftRanking);
  }

  inline AtomIndexType _rightHighPriority() const {
    return highPriority(_rightRanking);
  }

  inline AtomIndexType _rightLowPriority() const {
    return lowPriority(_rightRanking);
  }

  //! Generates the dihedral sequences of equal priority (i.e. high with high)
  std::vector<
    std::array<AtomIndexType, 4>
  > _equalPriorityDihedralSequences() const;

  //! Generates the dihedral sequences of unequal priortiy (i.e. high with low)
  std::vector<
    std::array<AtomIndexType, 4>
  > _differentPriorityDihedralSequences() const;

  /*!
   * Return the dihedral angle limits imposed by the underlying symmetries and
   * the current assignment.
   */
  std::vector<DihedralLimits> _dihedralLimits(
    const std::function<double(const AtomIndexType)> cycleMultiplierForIndex,
    const double looseningMultiplier
  ) const;

  static unsigned _numIndices(const RankingInformation& ranking);

public:
/* Constructors */
  EZStereocenter() = delete;
  EZStereocenter(
    const AtomIndexType& firstCenter,
    const RankingInformation& firstCenterRanking,
    const AtomIndexType& secondCenter,
    const RankingInformation& secondCenterRanking
  );

/* Modification */
  void addSubstituent(
    const AtomIndexType& center,
    RankingInformation centerRanking
  );

  void assign(boost::optional<unsigned> assignment) final;

  void assignRandom() final;

  void fit(const AngstromWrapper& angstromWrapper);

  void propagateGraphChange(
    RankingInformation firstCenterRanking,
    RankingInformation secondCenterRanking
  );

  void propagateVertexRemoval(const AtomIndexType removedIndex) final;

  void removeSubstituent(
    const AtomIndexType center,
    const AtomIndexType which
  );

/* Information */
  boost::optional<unsigned> assigned() const final;

  boost::optional<unsigned> indexOfPermutation() const final;

  unsigned numAssignments() const final;

  unsigned numStereopermutations() const final;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints() const final;

  const RankingInformation& getLeftRanking() const;
  const RankingInformation& getRightRanking() const;

  std::string info() const final;

  std::string rankInfo() const;

  std::vector<AtomIndexType> involvedAtoms() const final;

  void setModelInformation(
    DistanceGeometry::MoleculeSpatialModel& model,
    const std::function<double(const AtomIndexType)> cycleMultiplierForIndex,
    const double looseningMultiplier
  ) const final;

  Type type() const final;

/* Operators */
  bool operator == (const EZStereocenter& other) const;
  bool operator < (const EZStereocenter& other) const;
};

} // namespace Stereocenters

} // namespace molassembler

#endif
