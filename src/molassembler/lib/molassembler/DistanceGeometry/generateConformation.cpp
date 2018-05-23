#include <dlib/optimization.h>
#include <Eigen/Dense>
#include "temple/Containers.h"
#include "temple/Random.h"

#include "DistanceGeometry/dlibAdaptors.h"
#include "DistanceGeometry/dlibDebugAdaptors.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/RefinementProblem.h"
#include "DistanceGeometry/Error.h"
#include "DistanceGeometry/ExplicitGraph.h"

#include "GraphAlgorithms.h"

/* TODO
 * - Random assignment of unassigned stereocenters isn't ideal yet
 *   - Sequence randomization
 */

namespace molassembler {

namespace DistanceGeometry {

namespace predicates {

bool hasZeroPermutationsStereocenters(const Molecule& molecule) {
  return temple::any_of(
    molecule.getStereocenterList(),
    [](const auto& stereocenterPtr) -> bool {
      return stereocenterPtr -> numStereopermutations() == 0;
    }
  );
}

bool hasUnassignedStereocenters(const Molecule& mol) {
  return temple::any_of(
    mol.getStereocenterList(),
    [](const auto& stereocenterPtr) -> bool {
      return !stereocenterPtr -> assigned();
    }
  );
}

} // namespace predicates

namespace detail {

Delib::PositionCollection convertToPositionCollection(
  const Eigen::VectorXd& vectorizedPositions
) {
  const unsigned dimensionality = 4;
  assert(vectorizedPositions.size() % dimensionality == 0);

  Delib::PositionCollection positions;
  const unsigned N = vectorizedPositions.size() / dimensionality;

  for(unsigned i = 0; i < N; i++) {
    positions.push_back(
      Delib::Position {
        vectorizedPositions.template segment<3>(
          static_cast<Eigen::Index>(dimensionality) * i
        )
      }
    );
  }

  return positions;
}

Delib::PositionCollection convertToPositionCollection(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
) {
  const Eigen::Index dimensionality = 4;
  assert(vectorizedPositions.size() % dimensionality == 0);

  Delib::PositionCollection positions;
  const unsigned N = vectorizedPositions.size() / dimensionality;
  for(unsigned i = 0; i < N; i++) {
    positions.push_back(
      Delib::Position {
        vectorizedPositions(dimensionality * i),
        vectorizedPositions(dimensionality * i + 1),
        vectorizedPositions(dimensionality * i + 2)
      }
    );
  }

  return positions;
}


outcome::result<ChiralityConstraint> propagate(
  const DistanceBoundsMatrix& bounds,
  const Stereocenters::ChiralityConstraintPrototype& prototype
) {
  if(prototype.target == Stereocenters::ChiralityConstraintTarget::Flat) {
    // TODO introduce flatness tolerance
    return ChiralityConstraint {
      prototype.indices,
      0.0,
      0.0
    };
  }

  /* Calculate the upper and lower bounds on the volume. The target volume of
   * the chirality constraint created by the tetrahedron is calculated using
   * internal coordinates (the Cayley-Menger determinant), always leading to V
   * > 0, so depending on the current assignment, the sign of the result is
   * switched.  The formula used later in chirality constraint calculation for
   * explicit coordinates is adjusted by V' = 6 V to avoid an unnecessary
   * factor, so we do that here too:
   *
   *    288 V²  = |...|               | substitute V' = 6 V
   * -> 8 (V')² = |...|
   * ->      V' = sqrt(|...| / 8)
   *
   * where the Cayley-Menger determinant |...| is square symmetric:
   *
   *          |   0    1    1    1    1  |
   *          |        0  d12² d13² d14² |
   *  |...| = |             0  d23² d24² |
   *          |                  0  d34² |
   *          |  ...                  0  |
   *
   */

  Eigen::Matrix<double, 5, 5> lowerMatrix, upperMatrix;

  // Order of operations here is very important!

  lowerMatrix.row(0).setOnes();
  upperMatrix.row(0).setOnes();

  lowerMatrix.diagonal().setZero();
  upperMatrix.diagonal().setZero();

  for(unsigned i = 0; i < 4; i++) {
    for(unsigned j = i + 1; j < 4; j++) {
      const double& lowerBound = bounds.lowerBound(
        prototype.indices.at(i),
        prototype.indices.at(j)
      );

      const double& upperBound = bounds.upperBound(
        prototype.indices.at(i),
        prototype.indices.at(j)
      );

      lowerMatrix(
        static_cast<Eigen::Index>(i) + 1,
        static_cast<Eigen::Index>(j) + 1
      ) = lowerBound * lowerBound;

      upperMatrix(
        static_cast<Eigen::Index>(i) + 1,
        static_cast<Eigen::Index>(j) + 1
      ) = upperBound * upperBound;
    }
  }

  const double boundFromLower = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    lowerMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

  const double boundFromUpper = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    upperMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

#ifndef NDEBUG
  // Log in Debug builds, provided particular is set
  Log::log(Log::Particulars::PrototypePropagatorDebugInfo)
    << "Lower bound Cayley-menger matrix (determinant=" << boundFromLower <<  "): " << std::endl << lowerMatrix << std::endl
    << "Upper bound Cayley-menger matrix (determinant=" << boundFromUpper <<  "): " << std::endl << upperMatrix << std::endl;
#endif

  /* Check to ensure a volume is constructible from the passed bounds.
   * If the Cayley-Menger determinant is < 0, then it is impossible. The rather
   * lax criterion below follows from determinant calculation errors. Rather
   * be tolerant of approximately zero-volume chirality constraints (which MAY
   * happen) than have this fail.
   */
  if(boundFromLower < -1e-7 || boundFromUpper < -1e-7) {
    return DGError::GraphImpossible;
  }

  /* Although it is tempting to assume that the Cayley-Menger determinant using
   * the lower bounds is smaller than the one using upper bounds, this is NOT
   * true. We cannot a priori know which of both yields the lower or upper bounds
   * on the 3D volume, and hence must ensure only that the ordering is preserved
   * in the generation of the ChiralityConstraint, which checks that the lower
   * bound on the volume is certainly smaller than the upper one.
   *
   * You can check this assertion with a CAS. The relationship between both
   * determinants (where u_ij = l_ij + Δ) is wholly indeterminant, i.e. no
   * logical operator (<, >, <=, >=, ==) between both is true. It completely
   * depends on the individual values. Maybe in very specific cases one can
   * deduce some relationship, but not generally.
   */

  /* Abs is necessary to ensure that negative almost-zeros are interpreted as
   * positive near-zeros. Additionally, we do the adjustment by factor 8 here
   * (see above)
   */
  const double volumeFromLower = std::sqrt(std::fabs(boundFromLower) / 8);
  const double volumeFromUpper = std::sqrt(std::fabs(boundFromUpper) / 8);

  if(prototype.target == Stereocenters::ChiralityConstraintTarget::Negative) {
    /* For negative case, be aware inversion around 0 (*= -1) flips the inequality
     *      lower <  upper | *(-1)
     * ->  -lower > -upper
     *
     * So we construct it with the previous upper bound negated as the lower bound
     * and similarly for the previous lower bound.
     */
    return ChiralityConstraint {
      prototype.indices,
      - std::max(volumeFromLower, volumeFromUpper),
      - std::min(volumeFromLower, volumeFromUpper),
    };
  }

  // Regular case (Positive target)
  return ChiralityConstraint {
    prototype.indices,
    std::min(volumeFromLower, volumeFromUpper),
    std::max(volumeFromLower, volumeFromUpper),
  };
}

bool exceededFailureRatio(
  const unsigned& failures,
  const unsigned& numStructures,
  const unsigned& failureRatio
) noexcept {
  return static_cast<double>(failures) / numStructures >= failureRatio;
}

// Non-debug version of DG
outcome::result<
  std::vector<Delib::PositionCollection>
> runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const Partiality& metrizationOption,
  const bool& useYInversionTrick,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod
) {
  // TODO where to split control depending on loosening factor of spatial model?

  if(predicates::hasZeroPermutationsStereocenters(molecule)) {
    return DGError::ZeroPermutationStereocenters;
  }

  MoleculeDGInformation DGData;
  /* In case the molecule has unassigned stereocenters, we need to randomly
   * assign them in every step prior to generating the distance bounds matrix
   */
  bool regenerateEachStep = predicates::hasUnassignedStereocenters(molecule);
  if(!regenerateEachStep) {
    DGData = gatherDGInformation(
      molecule,
      distanceMethod
    );
  }

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function returns an error code. Remember that if an optimization is
   * considered a failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 0.1; // allow only 10% failures in release
  unsigned failures = 0;
  std::vector<Delib::PositionCollection> ensemble;
  ensemble.reserve(numStructures);

  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numStructures;
    currentStructureNumber += 1
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = molecule;

      do {
        for(const auto& stereocenterPtr : moleculeCopy.getStereocenterList()) {
          if(!stereocenterPtr->assigned()) {
            moleculeCopy.assignStereocenterRandomly(
              stereocenterPtr->involvedAtoms().front()
            );

            /* Jump out of for-loop and re-check if there are still unassigned
             * stereocenters since assigning a stereocenter can invalidate
             * the list iterators
             */
            break;
          }
        }
      } while(predicates::hasUnassignedStereocenters(moleculeCopy));

      // Fetch the DG data from the molecule with no unassigned stereocenters
      DGData = gatherDGInformation(
        moleculeCopy,
        distanceMethod
      );
    }

    ExplicitGraph explicitGraph {
      molecule,
      DGData.boundList
    };

    auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
    if(!distanceBoundsResult) {
      return distanceBoundsResult.as_failure();
    }

    DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

    /* No need to smooth the distance bounds, the graph type creates them
     * within the triangle inequality bounds
     */

    assert(distanceBounds.boundInconsistencies() == 0);

    auto distanceMatrixResult = explicitGraph.makeDistanceMatrix(metrizationOption);
    if(!distanceMatrixResult) {
      return distanceMatrixResult.as_failure();
    }

    // Make a metric matrix from a generated distances matrix
    MetricMatrix metric(
      std::move(distanceMatrixResult.value())
    );

    /* Get the chirality constraints by converting the prototypes found by the
     * collector into full ChiralityConstraints
     */
    std::vector<ChiralityConstraint> chiralityConstraints;

    for(const auto& prototype : DGData.chiralityConstraintPrototypes) {
      if(auto propagateResult = propagate(distanceBounds, prototype)) {
        chiralityConstraints.push_back(
          std::move(propagateResult.value())
        );
      } else {
        return propagateResult.as_failure();
      }
    }

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize and transfer to dlib
    errfValue<true>::Vector dlibPositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      )
    );

    if(useYInversionTrick) {
      /* If a count of chirality constraints reveals that more than half are
       * incorrect, we can invert the structure (by multiplying e.g. all y
       * coordinates with -1) and then have more than half of chirality
       * constraints correct! In the count, chirality constraints with a target
       * value of zero are not considered (this would skew the count as those
       * chirality constraints should not have to pass an energetic maximum to
       * converge properly as opposed to tetrahedra with volume).
       */
      if(
        errfDetail::proportionChiralityConstraintsCorrectSign(
          chiralityConstraints,
          dlibPositions
        ) < 0.5
      ) {
        // Invert y coordinates
        for(unsigned i = 0; i < distanceBounds.N(); ++i) {
          dlibPositions(
            static_cast<dlibIndexType>(i) * 4  + 1
          ) *= -1;
        }
      }
    }

    // Create the squared bounds matrix for use in refinement
    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      static_cast<Eigen::MatrixXd>(
        distanceBounds.access().cwiseProduct(distanceBounds.access())
      )
    );

    /* Refinement without penalty on fourth dimension only necessary if not all
     * chiral centers are correct. Of course, for molecules without chiral
     * centers at all, this stage is unnecessary
     */
    if(
      errfDetail::proportionChiralityConstraintsCorrectSign(
        chiralityConstraints,
        dlibPositions
      ) < 1
    ) {
      dlibAdaptors::iterationOrAllChiralitiesCorrectStrategy inversionStopStrategy {
        chiralityConstraints,
        10000
      };

      dlib::find_min(
        dlib::bfgs_search_strategy(),
        inversionStopStrategy,
        errfValue<false>(
          squaredBounds,
          chiralityConstraints
        ),
        errfGradient<false>(
          squaredBounds,
          chiralityConstraints
        ),
        dlibPositions,
        0
      );

      if(
        inversionStopStrategy.iterations > 10000
        || errfDetail::proportionChiralityConstraintsCorrectSign(
          chiralityConstraints,
          dlibPositions
        ) < 1.0
      ) { // Failure to invert
        failures += 1;
        if(exceededFailureRatio(failures, numStructures, failureRatio)) {
          return DGError::TooManyFailures;
        }
        continue; // this triggers a new structure to be generated
      }
    }

    dlibAdaptors::iterationOrGradientNormStopStrategy refinementStopStrategy {
      10000, // iteration limit
      1e-5 // gradient limit
    };

    // Refinement with penalty on fourth dimension is always necessary
    dlib::find_min(
      dlib::bfgs_search_strategy(),
      refinementStopStrategy,
      errfValue<true>(
        squaredBounds,
        chiralityConstraints
      ),
      errfGradient<true>(
        squaredBounds,
        chiralityConstraints
      ),
      dlibPositions,
      0
    );

    bool reachedMaxIterations = refinementStopStrategy.iterations >= 10000;
    bool notAllChiralitiesCorrect = errfDetail::proportionChiralityConstraintsCorrectSign(
      chiralityConstraints,
      dlibPositions
    ) < 1;
    bool structureAcceptable = errfDetail::finalStructureAcceptable(
      distanceBounds,
      chiralityConstraints,
      dlibPositions
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      failures += 1;

      if(exceededFailureRatio(failures, numStructures, failureRatio)) {
        return DGError::TooManyFailures;
      }
    } else {
      ensemble.emplace_back(
        detail::convertToPositionCollection(dlibPositions)
      );
    }
  }

  return ensemble;
}

// Debug version
std::list<RefinementData> debugDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const Partiality& metrizationOption,
  const bool& useYInversionTrick,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod
) {
  if(predicates::hasZeroPermutationsStereocenters(molecule)) {
    Log::log(Log::Level::Warning)
      << "This molecule has stereocenters with zero valid permutations!"
      << std::endl;
  }

  std::list<RefinementData> refinementList;

  /* In case the molecule has unassigned stereocenters that are not trivially
   * assignable (u/1 -> 0/1), random assignments have to be made prior to calling
   * gatherDGInformation (which creates the DistanceBoundsMatrix via the
   * MoleculeSpatialModel, which expects all stereocenters to be assigned).
   * Accordingly, gatherDGInformation has to be repeated in those cases, while
   * it is necessary only once in the other
   */

  /* There is also some degree of doubt about the relative frequencies of
   * assignments in stereopermutation. I should get on that too.
   */

  bool regenerateEachStep = predicates::hasUnassignedStereocenters(molecule);

  MoleculeDGInformation DGData;

  if(!regenerateEachStep) { // Collect once, keep all the time
    DGData = gatherDGInformation(
      molecule,
      distanceMethod
    );
  }

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a
   * failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 3; // more lenient than production, must be > 0
  unsigned failures = 0;

  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numStructures;
    currentStructureNumber += 1
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = molecule;

      do {
        for(const auto& stereocenterPtr : moleculeCopy.getStereocenterList()) {
          if(!stereocenterPtr->assigned()) {
            moleculeCopy.assignStereocenterRandomly(
              stereocenterPtr->involvedAtoms().front()
            );

            /* Jump out of for-loop and re-check if there are still unassigned
             * stereocenters since assigning a stereocenter can invalidate
             * the list iterators (because stereocenters may appear or disappear
             * on assignment)
             */
            break;
          }
        }
      } while(predicates::hasUnassignedStereocenters(moleculeCopy));

      // Fetch the DG data from the molecule with no unassigned stereocenters
      DGData = gatherDGInformation(
        moleculeCopy,
        distanceMethod
      );
    }

    std::list<RefinementStepData> refinementSteps;

    ExplicitGraph explicitGraph {
      molecule,
      DGData.boundList
    };

    auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
    if(!distanceBoundsResult) {
      Log::log(Log::Level::Warning) << "Failure in distance bounds matrix construction.\n";
      failures += 1;
      if(exceededFailureRatio(failures, numStructures, failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
        return refinementList;
      }

      continue;
    }

    DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

    /* No need to smooth the distance bounds, ExplicitGraph creates it
     * so that the triangle inequalities are fulfilled
     */

    assert(distanceBounds.boundInconsistencies() == 0);

    auto distanceMatrixResult = explicitGraph.makeDistanceMatrix(metrizationOption);
    if(!distanceMatrixResult) {
      Log::log(Log::Level::Warning) << "Failure in distance matrix construction.\n";
      failures += 1;
      if(exceededFailureRatio(failures, numStructures, failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
        return refinementList;
      }

      continue;
    }

    // Make a metric matrix from a generated distances matrix
    MetricMatrix metric(
      std::move(distanceMatrixResult.value())
    );

    /* Get the chirality constraints by converting the prototypes found by the
     * collector into full chiralityConstraints
     */
    std::vector<ChiralityConstraint> chiralityConstraints;

    for(const auto& prototype : DGData.chiralityConstraintPrototypes) {
      if(auto propagateResult = propagate(distanceBounds, prototype)) {
        chiralityConstraints.push_back(
          std::move(propagateResult.value())
        );
      } else {
        Log::log(Log::Level::Warning) << "Propagation of chirality constraint failed!\n";
        return refinementList;
      }
    }

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize and transfer to dlib
    errfValue<true>::Vector dlibPositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      )
    );

    if(useYInversionTrick) {
      /* If a count of chirality constraints reveals that more than half are
       * incorrect, we can invert the structure (by multiplying e.g. all y
       * coordinates with -1) and then have more than half of chirality
       * constraints correct! In the count, chirality constraints with a target
       * value of zero are not considered (this would skew the count as those
       * chirality constraints should not have to pass an energetic maximum to
       * converge properly as opposed to tetrahedra with volume).
       */
      if(
        errfDetail::proportionChiralityConstraintsCorrectSign(
          chiralityConstraints,
          dlibPositions
        ) < 0.5
      ) {
        // Invert y coordinates
        for(unsigned i = 0; i < distanceBounds.N(); i++) {
          dlibPositions(
            static_cast<dlibIndexType>(i) * 4 + 1
          ) *= -1;
        }
      }
    }

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      static_cast<Eigen::MatrixXd>(
        distanceBounds.access().cwiseProduct(distanceBounds.access())
      )
    );

    /* Refinement without penalty on fourth dimension only necessary if not all
     * chiral centers are correct. Of course, for molecules without chiral
     * centers at all, this stage is unnecessary
     */
    if(
      errfDetail::proportionChiralityConstraintsCorrectSign(
        chiralityConstraints,
        dlibPositions
      ) < 1
    ) {

      errfValue<false> valueFunctor {
        squaredBounds,
        chiralityConstraints
      };

      dlibAdaptors::debugIterationOrAllChiralitiesCorrectStrategy inversionStopStrategy {
        10000,
        refinementSteps,
        valueFunctor
      };

      dlib::find_min(
        dlib::bfgs_search_strategy(),
        inversionStopStrategy,
        valueFunctor,
        errfGradient<false>(
          squaredBounds,
          chiralityConstraints
        ),
        dlibPositions,
        0
      );

      if(
        inversionStopStrategy.iterations > 10000
        || errfDetail::proportionChiralityConstraintsCorrectSign(
          chiralityConstraints,
          dlibPositions
        ) < 1.0
      ) { // Failure to invert
        Log::log(Log::Level::Warning) << "First stage of refinement fails.\n";
        failures += 1;
        if(exceededFailureRatio(failures, numStructures, failureRatio)) {
          Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
          return refinementList;
        }
        continue; // this triggers a new structure to be generated
      }
    }

    errfValue<true> refinementValueFunctor {
      squaredBounds,
      chiralityConstraints
    };

    dlibAdaptors::debugIterationOrGradientNormStopStrategy refinementStopStrategy {
      10000, // iteration limit
      1e-5, // gradient limit
      refinementSteps,
      refinementValueFunctor
    };

    // Refinement with penalty on fourth dimension is always necessary
    dlib::find_min(
      dlib::bfgs_search_strategy(),
      refinementStopStrategy,
      refinementValueFunctor,
      errfGradient<true>(
        squaredBounds,
        chiralityConstraints
      ),
      dlibPositions,
      0
    );

    bool reachedMaxIterations = refinementStopStrategy.iterations >= 10000;
    bool notAllChiralitiesCorrect = errfDetail::proportionChiralityConstraintsCorrectSign(
      chiralityConstraints,
      dlibPositions
    ) < 1;
    bool structureAcceptable = errfDetail::finalStructureAcceptable(
      distanceBounds,
      chiralityConstraints,
      dlibPositions
    );

    RefinementData refinementData;
    refinementData.steps = std::move(refinementSteps);
    refinementData.constraints = chiralityConstraints;
    refinementData.isFailure = (reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable);

    refinementList.push_back(
      std::move(refinementData)
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      Log::log(Log::Level::Warning) << "Second stage of refinement fails.\n";
      if(reachedMaxIterations) {
        Log::log(Log::Level::Warning) << "- Reached max iterations.\n";
      }

      if(notAllChiralitiesCorrect) {
        Log::log(Log::Level::Warning) << "- Not all chirality constraints have the correct sign.\n";
      }

      if(!structureAcceptable) {
        Log::log(Log::Level::Warning) << "- The final structure is unacceptable.\n";
        if(Log::particulars.count(Log::Particulars::DGStructureAcceptanceFailures)) {
          errfDetail::explainAcceptanceFailure(
            distanceBounds,
            chiralityConstraints,
            dlibPositions
          );
        }
      }

      failures += 1;
      if(exceededFailureRatio(failures, numStructures, failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
        return refinementList;
      }
    }
  }

  return refinementList;
}

} // namespace detail

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const MoleculeSpatialModel::DistanceMethod& distanceMethod
) {
  // Initialize the return object
  MoleculeDGInformation data;

  // Generate a spatial model from the molecular graph and stereocenters
  MoleculeSpatialModel spatialModel {
    molecule,
    distanceMethod
  };

  // Generate the distance bounds from the spatial model
  spatialModel.addDefaultDihedrals();

  data.boundList = spatialModel.makeBoundList();

  // Extract the chirality constraint prototypes
  data.chiralityConstraintPrototypes = spatialModel.getChiralityPrototypes();

  return data;
}

outcome::result<
  std::vector<Delib::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures
) {
  return detail::runDistanceGeometry(molecule, numStructures);
}

outcome::result<Delib::PositionCollection> generateConformation(const Molecule& molecule) {
  if(auto result = detail::runDistanceGeometry(molecule, 1)) {
    const auto& conformationList = result.value();
    assert(conformationList.size() == 1);
    return conformationList.front();
  } else {
    return result.as_failure();
  }
}

} // namespace DistanceGeometry

} // namespace molassembler
