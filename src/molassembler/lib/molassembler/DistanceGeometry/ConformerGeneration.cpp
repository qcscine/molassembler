#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include <dlib/optimization.h>
#include <Eigen/Dense>
#include "Delib/Constants.h"

#include "molassembler/DistanceGeometry/dlibAdaptors.h"
#include "molassembler/DistanceGeometry/dlibDebugAdaptors.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/RefinementProblem.h"
#include "molassembler/DistanceGeometry/Error.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/Graph/GraphAlgorithms.h"

/* TODO
 * - Random assignment of unassigned stereocenters isn't ideal yet
 *   - Sequence randomization
 */

namespace molassembler {

namespace DistanceGeometry {

namespace detail {

AngstromWrapper convertToAngstromWrapper(
  const Eigen::VectorXd& vectorizedPositions
) {
  const unsigned dimensionality = 4;
  assert(vectorizedPositions.size() % dimensionality == 0);

  const unsigned N = vectorizedPositions.size() / dimensionality;
  AngstromWrapper angstromWrapper {N};

  for(unsigned i = 0; i < N; i++) {
    angstromWrapper.positions.at(i) = Delib::Position {
      vectorizedPositions.template segment<3>(
        static_cast<Eigen::Index>(dimensionality) * i
      )
    };
  }

  return angstromWrapper;
}

AngstromWrapper convertToAngstromWrapper(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
) {
  const Eigen::Index dimensionality = 4;
  assert(vectorizedPositions.size() % dimensionality == 0);

  const unsigned N = vectorizedPositions.size() / dimensionality;
  AngstromWrapper angstromWrapper {N};
  for(unsigned i = 0; i < N; i++) {
    angstromWrapper.positions.at(i) = Delib::Position {
      vectorizedPositions(dimensionality * i),
      vectorizedPositions(dimensionality * i + 1),
      vectorizedPositions(dimensionality * i + 2)
    };
  }

  return angstromWrapper;
}

Molecule narrow(Molecule moleculeCopy) {
  const auto& stereocenterList = moleculeCopy.stereocenters();

  do {
    /* If we change any stereocenter, we must jump out of for-loop and
     * re-check if there are still unassigned stereocenters since assigning
     * a stereocenter can invalidate stereocenter list iterators (because
     * stereocenters may appear or disappear on assignment due to ranking)
     */
    bool changedStereocenterFlag = false;

    for(const auto& atomStereocenter : stereocenterList.atomStereocenters()) {
      if(!atomStereocenter.assigned()) {
        moleculeCopy.assignStereocenterRandomly(
          atomStereocenter.centralIndex()
        );

        changedStereocenterFlag = true;
        break;
      }
    }

    if(changedStereocenterFlag) {
      // Re-check the do-while loop condition
      continue;
    }

    for(const auto& bondStereocenter : stereocenterList.bondStereocenters()) {
      if(!bondStereocenter.assigned()) {
        moleculeCopy.assignStereocenterRandomly(
          bondStereocenter.edge()
        );

        // Break this loop and re-check do-while loop condition
        break;
      }
    }
  } while(stereocenterList.hasUnassignedStereocenters());

  return moleculeCopy;
}

inline double looseningFactor(const unsigned failures, const unsigned numStructures) {
  // x ranges from 0 to failureRatio, which is currently 2
  double x = static_cast<double>(failures) / numStructures;

  /* Could consider an initially a stronger increase and then linear later, like
   * -1 / (x + 1)^2 + 3x / 8 + 2,
   *
   * but initially prefer linear from (0, 2) to (2, 2.5)
   */
  return 3 * x / 4 + 1;
}

inline bool exceededFailureRatio(
  const unsigned failures,
  const unsigned numStructures,
  const unsigned failureRatio
) noexcept {
  return static_cast<double>(failures) / numStructures >= failureRatio;
}

} // namespace detail

// Non-debug version of DG
outcome::result<
  std::vector<AngstromWrapper>
> run(
  const Molecule& molecule,
  const unsigned numStructures,
  const Partiality metrizationOption,
  const bool useYInversionTrick
) {
  if(molecule.stereocenters().hasZeroAssignmentStereocenters()) {
    return DGError::ZeroAssignmentStereocenters;
  }

  MoleculeDGInformation DGData;
  /* In case the molecule has unassigned stereocenters, we need to randomly
   * assign them in every step prior to generating the distance bounds matrix
   */
  bool regenerateEachStep = molecule.stereocenters().hasUnassignedStereocenters();
  if(!regenerateEachStep) {
    DGData = gatherDGInformation(molecule);
  }

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function returns an error code. Remember that if an optimization is
   * considered a failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 2; // allow for max 2x #conformations embedding failures
  unsigned failures = 0;
  bool loosenBounds = false;
  std::vector<AngstromWrapper> ensemble;
  ensemble.reserve(numStructures);

  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numStructures;
    currentStructureNumber += 1
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = detail::narrow(molecule);

      if(moleculeCopy.stereocenters().hasZeroAssignmentStereocenters()) {
        return DGError::ZeroAssignmentStereocenters;
      }

      // Fetch the DG data from the molecule with no unassigned stereocenters
      DGData = gatherDGInformation(
        moleculeCopy,
        detail::looseningFactor(failures, numStructures)
      );
    } else if(loosenBounds) {
      DGData = gatherDGInformation(
        molecule,
        detail::looseningFactor(failures, numStructures)
      );

      loosenBounds = false;
    }

    ExplicitGraph explicitGraph {
      molecule,
      DGData.bounds
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
          DGData.chiralityConstraints,
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
        DGData.chiralityConstraints,
        dlibPositions
      ) < 1
    ) {
      dlibAdaptors::iterationOrAllChiralitiesCorrectStrategy inversionStopStrategy {
        DGData.chiralityConstraints,
        10000
      };

      dlib::find_min(
        dlib::bfgs_search_strategy(),
        inversionStopStrategy,
        errfValue<false>(
          squaredBounds,
          DGData.chiralityConstraints
        ),
        errfGradient<false>(
          squaredBounds,
          DGData.chiralityConstraints
        ),
        dlibPositions,
        0
      );

      if(
        inversionStopStrategy.iterations > 10000
        || errfDetail::proportionChiralityConstraintsCorrectSign(
          DGData.chiralityConstraints,
          dlibPositions
        ) < 1.0
      ) { // Failure to invert
        failures += 1;
        loosenBounds = true;
        if(detail::exceededFailureRatio(failures, numStructures, failureRatio)) {
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
        DGData.chiralityConstraints
      ),
      errfGradient<true>(
        squaredBounds,
        DGData.chiralityConstraints
      ),
      dlibPositions,
      0
    );

    bool reachedMaxIterations = refinementStopStrategy.iterations >= 10000;
    bool notAllChiralitiesCorrect = errfDetail::proportionChiralityConstraintsCorrectSign(
      DGData.chiralityConstraints,
      dlibPositions
    ) < 1;
    bool structureAcceptable = errfDetail::finalStructureAcceptable(
      distanceBounds,
      DGData.chiralityConstraints,
      dlibPositions
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      failures += 1;
      loosenBounds = true;

      if(detail::exceededFailureRatio(failures, numStructures, failureRatio)) {
        return DGError::TooManyFailures;
      }
    } else {
      ensemble.emplace_back(
        detail::convertToAngstromWrapper(dlibPositions)
      );
    }
  }

  return ensemble;
}

// Debug version
std::list<RefinementData> debug(
  const Molecule& molecule,
  const unsigned numStructures,
  const Partiality metrizationOption,
  const bool useYInversionTrick
) {
  if(molecule.stereocenters().hasZeroAssignmentStereocenters()) {
    Log::log(Log::Level::Warning)
      << "This molecule has stereocenters with zero valid permutations!"
      << std::endl;
  }

  std::list<RefinementData> refinementList;

  /* In case the molecule has unassigned stereocenters that are not trivially
   * assignable (u/1 -> 0/1), random assignments have to be made prior to calling
   * gatherDGInformation (which creates the DistanceBoundsMatrix via the
   * SpatialModel, which expects all stereocenters to be assigned).
   * Accordingly, gatherDGInformation has to be repeated in those cases, while
   * it is necessary only once in the other
   */

  /* There is also some degree of doubt about the relative frequencies of
   * assignments in stereopermutation. I should get on that too.
   */

  bool regenerateEachStep = molecule.stereocenters().hasUnassignedStereocenters();

  MoleculeDGInformation DGData;
  std::string spatialModelGraphviz;

  if(!regenerateEachStep) { // Collect once, keep all the time
    DGData = gatherDGInformation(molecule, 1.0, spatialModelGraphviz);
  }

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a
   * failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 3; // more lenient than production, must be > 0
  unsigned failures = 0;
  bool loosenBounds = false;

  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numStructures;
    currentStructureNumber += 1
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = detail::narrow(molecule);

      if(moleculeCopy.stereocenters().hasZeroAssignmentStereocenters()) {
        Log::log(Log::Level::Warning)
          << "After setting stereocenters at random, this molecule has "
          << "stereocenters with zero valid permutations!"
          << std::endl;
      }

      // Fetch the DG data from the molecule with no unassigned stereocenters
      DGData = gatherDGInformation(
        moleculeCopy,
        detail::looseningFactor(failures, numStructures),
        spatialModelGraphviz
      );
    } else if(loosenBounds) {
      DGData = gatherDGInformation(
        molecule,
        detail::looseningFactor(failures, numStructures),
        spatialModelGraphviz
      );

      loosenBounds = false;
    }

    std::list<RefinementStepData> refinementSteps;

    ExplicitGraph explicitGraph {
      molecule,
      DGData.bounds
    };

    auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
    if(!distanceBoundsResult) {
      Log::log(Log::Level::Warning) << "Failure in distance bounds matrix construction: "
        << distanceBoundsResult.error().message() << "\n";
      failures += 1;
      loosenBounds = true;
      if(detail::exceededFailureRatio(failures, numStructures, failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG. Sample spatial model written to 'DG-failure-spatial-model.dot'.\n";
        if(regenerateEachStep) {
          auto moleculeCopy = detail::narrow(molecule);

          SpatialModel model {moleculeCopy};
          model.writeGraphviz("DG-failure-spatial-model.dot");
        } else {
          SpatialModel model {molecule};
          model.writeGraphviz("DG-failure-spatial-model.dot");
        }
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
      loosenBounds = true;
      if(detail::exceededFailureRatio(failures, numStructures, failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
        return refinementList;
      }

      continue;
    }

    // Make a metric matrix from a generated distances matrix
    MetricMatrix metric(
      std::move(distanceMatrixResult.value())
    );

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
          DGData.chiralityConstraints,
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
        DGData.chiralityConstraints,
        dlibPositions
      ) < 1
    ) {

      errfValue<false> valueFunctor {
        squaredBounds,
        DGData.chiralityConstraints
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
          DGData.chiralityConstraints
        ),
        dlibPositions,
        0
      );

      if(
        inversionStopStrategy.iterations > 10000
        || errfDetail::proportionChiralityConstraintsCorrectSign(
          DGData.chiralityConstraints,
          dlibPositions
        ) < 1.0
      ) { // Failure to invert
        Log::log(Log::Level::Warning)
          << "[" << currentStructureNumber << "]: "
          << "First stage of refinement fails. Loosening factor was "
          << detail::looseningFactor(failures, numStructures)
          <<  "\n";
        failures += 1;
        loosenBounds = true;
        if(detail::exceededFailureRatio(failures, numStructures, failureRatio)) {
          Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
          return refinementList;
        }
        continue; // this triggers a new structure to be generated
      }
    }

    errfValue<true> refinementValueFunctor {
      squaredBounds,
      DGData.chiralityConstraints
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
        DGData.chiralityConstraints
      ),
      dlibPositions,
      0
    );

    bool reachedMaxIterations = refinementStopStrategy.iterations >= 10000;
    bool notAllChiralitiesCorrect = errfDetail::proportionChiralityConstraintsCorrectSign(
      DGData.chiralityConstraints,
      dlibPositions
    ) < 1;
    bool structureAcceptable = errfDetail::finalStructureAcceptable(
      distanceBounds,
      DGData.chiralityConstraints,
      dlibPositions
    );

    if(Log::particulars.count(Log::Particulars::DGFinalErrorContributions)) {
      errfDetail::explainFinalContributions(
        distanceBounds,
        DGData.chiralityConstraints,
        dlibPositions
      );
    }

    RefinementData refinementData;
    refinementData.steps = std::move(refinementSteps);
    refinementData.constraints = DGData.chiralityConstraints;
    refinementData.looseningFactor = detail::looseningFactor(failures, numStructures);
    refinementData.isFailure = (reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable);
    refinementData.spatialModelGraphviz = spatialModelGraphviz;

    refinementList.push_back(
      std::move(refinementData)
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      Log::log(Log::Level::Warning)
        << "[" << currentStructureNumber << "]: "
        << "Second stage of refinement fails. Loosening factor was "
        << detail::looseningFactor(failures, numStructures)
        <<  "\n";
      if(reachedMaxIterations) {
        Log::log(Log::Level::Warning) << "- Reached max iterations.\n";
      }

      if(notAllChiralitiesCorrect) {
        Log::log(Log::Level::Warning) << "- Not all chirality constraints have the correct sign.\n";
      }

      if(!structureAcceptable) {
        Log::log(Log::Level::Warning) << "- The final structure is unacceptable.\n";
        if(Log::isSet(Log::Particulars::DGStructureAcceptanceFailures)) {
          errfDetail::explainAcceptanceFailure(
            distanceBounds,
            DGData.chiralityConstraints,
            dlibPositions
          );
        }
      }

      failures += 1;
      loosenBounds = true;
      if(detail::exceededFailureRatio(failures, numStructures, failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
        return refinementList;
      }
    }
  }

  return refinementList;
}

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const double looseningFactor
) {
  // Generate a spatial model from the molecular graph and stereocenters
  SpatialModel spatialModel {molecule, looseningFactor};

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makeBoundsList();
  data.chiralityConstraints = spatialModel.getChiralityConstraints();

  return data;
}

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const double looseningFactor,
  std::string& spatialModelGraphvizString
) {
  // Generate a spatial model from the molecular graph and stereocenters
  SpatialModel spatialModel {molecule, looseningFactor};
  spatialModelGraphvizString = spatialModel.dumpGraphviz();

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makeBoundsList();
  data.chiralityConstraints = spatialModel.getChiralityConstraints();

  return data;
}

} // namespace DistanceGeometry

} // namespace molassembler
