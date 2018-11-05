// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

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
#include "molassembler/Utils/QuaternionFit.h"

namespace molassembler {

namespace DistanceGeometry {

namespace detail {

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

AngstromWrapper convertToAngstromWrapper(const Eigen::MatrixXd& positionsMatrix) {
  assert(positionsMatrix.cols() == 3);
  const unsigned N = positionsMatrix.rows();
  AngstromWrapper angstromWrapper {N};
  for(unsigned i = 0; i < N; ++i) {
    angstromWrapper.positions.at(i) = Delib::Position {
      positionsMatrix.row(i).transpose()
    };
  }
  return angstromWrapper;
}

Eigen::MatrixXd fitAndSetFixedPositions(
  const dlib::matrix<double, 0, 1>& vectorizedPositions,
  const Configuration& configuration
) {
  /* Fixed positions postprocessing:
   * - Rotate and translate the generated coordinates towards the fixed
   *   positions indicated for each.
   * - Assuming the fit isn't absolutely exact, overwrite the existing
   *   positions with the fixed ones
   */
  // Convert the vectorized dlib positions into a positions matrix
  const unsigned dimensionality = 4;
  const unsigned N = vectorizedPositions.size() / dimensionality;
  Eigen::MatrixXd positionMatrix(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    positionMatrix.row(i) = Eigen::RowVector3d(
      vectorizedPositions(dimensionality * i),
      vectorizedPositions(dimensionality * i + 1),
      vectorizedPositions(dimensionality * i + 2)
    );
  }

  /* Construct a reference matrix from the fixed positions (these are still in
   * bohr!) and a weights vector
   */
  Eigen::MatrixXd referenceMatrix = Eigen::MatrixXd::Zero(N, 3);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
  for(const auto& indexPositionPair : configuration.fixedPositions) {
    referenceMatrix.row(indexPositionPair.first) = indexPositionPair.second.asEigenVector().transpose();
    weights(indexPositionPair.first) = 1;
  }

  referenceMatrix *= Delib::angstrom_per_bohr;

  // Perform the QuaternionFit
  Scine::Utils::QuaternionFit fit(referenceMatrix, positionMatrix, weights);

  return fit.getFittedData();
}

Molecule narrow(Molecule molecule) {
  const auto& stereopermutatorList = molecule.stereopermutators();

  do {
    /* If we change any stereopermutator, we must re-check if there are still
     * unassigned stereopermutators since assigning a stereopermutator can
     * invalidate the entire stereopermutator list (because stereopermutators
     * may appear or disappear on assignment due to ranking)
     */
    std::vector<AtomIndex> unassignedAtomStereopermutators;

    for(const auto& atomStereopermutator : stereopermutatorList.atomStereopermutators()) {
      if(!atomStereopermutator.assigned()) {
        unassignedAtomStereopermutators.push_back(
          atomStereopermutator.centralIndex()
        );
      }
    }

    if(!unassignedAtomStereopermutators.empty()) {
      unsigned choice = temple::random::getSingle<unsigned>(
        0,
        unassignedAtomStereopermutators.size() - 1,
        randomnessEngine()
      );

      molecule.assignStereopermutatorRandomly(
        unassignedAtomStereopermutators.at(choice)
      );

      // Re-check the loop condition
      continue;
    }

    std::vector<BondIndex> unassignedBondStereopermutators;

    for(const auto& bondStereopermutator : stereopermutatorList.bondStereopermutators()) {
      if(!bondStereopermutator.assigned()) {
        unassignedBondStereopermutators.push_back(
          bondStereopermutator.edge()
        );
      }
    }

    if(!unassignedBondStereopermutators.empty()) {
      unsigned choice = temple::random::getSingle<unsigned>(
        0,
        unassignedBondStereopermutators.size() - 1,
        randomnessEngine()
      );

      molecule.assignStereopermutatorRandomly(
        unassignedBondStereopermutators.at(choice)
      );
    }

  } while(stereopermutatorList.hasUnassignedStereopermutators());

  return molecule;
}

inline double looseningFactor(
  const unsigned failures,
  const unsigned numConformers,
  const double failureRatio
) {
  // x ranges from 0 to failureRatio
  double x = static_cast<double>(failures) / numConformers;

  return (1 + x / failureRatio);
}

inline bool exceededFailureRatio(
  const unsigned failures,
  const unsigned numConformers,
  const unsigned failureRatio
) noexcept {
  return static_cast<double>(failures) / numConformers >= failureRatio;
}

} // namespace detail

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration,
  const double looseningFactor
) {
  // Generate a spatial model from the molecular graph and stereopermutators
  SpatialModel spatialModel {molecule, configuration, looseningFactor};

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makeBoundsList();
  data.chiralityConstraints = spatialModel.getChiralityConstraints();

  return data;
}

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration,
  const double looseningFactor,
  std::string& spatialModelGraphvizString
) {
  // Generate a spatial model from the molecular graph and stereopermutators
  SpatialModel spatialModel {molecule, configuration, looseningFactor};
  spatialModelGraphvizString = spatialModel.dumpGraphviz();

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makeBoundsList();
  data.chiralityConstraints = spatialModel.getChiralityConstraints();

  return data;
}

// Debug version
std::list<RefinementData> debug(
  const Molecule& molecule,
  const unsigned numConformers,
  const Configuration& configuration
) {
  if(molecule.stereopermutators().hasZeroAssignmentStereopermutators()) {
    Log::log(Log::Level::Warning)
      << "This molecule has stereopermutators with zero valid permutations!"
      << std::endl;
  }

  SpatialModel::checkFixedPositionsPreconditions(molecule, configuration);

  std::list<RefinementData> refinementList;

  /* In case the molecule has unassigned stereopermutators that are not trivially
   * assignable (u/1 -> 0/1), random assignments have to be made prior to calling
   * gatherDGInformation (which creates the DistanceBoundsMatrix via the
   * SpatialModel, which expects all stereopermutators to be assigned).
   * Accordingly, gatherDGInformation has to be repeated in those cases, while
   * it is necessary only once in the other
   */

  /* There is also some degree of doubt about the relative frequencies of
   * assignments in stereopermutation. I should get on that too.
   */

  bool regenerateEachStep = molecule.stereopermutators().hasUnassignedStereopermutators();

  MoleculeDGInformation DGData;
  std::string spatialModelGraphviz;

  if(!regenerateEachStep) { // Collect once, keep all the time
    DGData = gatherDGInformation(molecule, configuration, 1.0, spatialModelGraphviz);
  }

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a
   * failure is dependent only on the stopping criteria!
   */
  unsigned failures = 0;
  bool loosenBounds = false;

  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numConformers;
    currentStructureNumber += 1
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = detail::narrow(molecule);

      if(moleculeCopy.stereopermutators().hasZeroAssignmentStereopermutators()) {
        Log::log(Log::Level::Warning)
          << "After setting stereopermutators at random, this molecule has "
          << "stereopermutators with zero valid permutations!"
          << std::endl;
      }

      // Fetch the DG data from the molecule with no unassigned stereopermutators
      DGData = gatherDGInformation(
        moleculeCopy,
        configuration,
        detail::looseningFactor(failures, numConformers, configuration.failureRatio),
        spatialModelGraphviz
      );
    } else if(loosenBounds) {
      DGData = gatherDGInformation(
        molecule,
        configuration,
        detail::looseningFactor(failures, numConformers, configuration.failureRatio),
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
      if(detail::exceededFailureRatio(failures, numConformers, configuration.failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG. Sample spatial model written to 'DG-failure-spatial-model.dot'.\n";

        double looseningFactor = detail::looseningFactor(failures, numConformers, configuration.failureRatio);

        if(regenerateEachStep) {
          auto moleculeCopy = detail::narrow(molecule);

          SpatialModel model {moleculeCopy, configuration, looseningFactor};
          model.writeGraphviz("DG-failure-spatial-model.dot");
        } else {
          SpatialModel model {molecule, configuration, looseningFactor};
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

    auto distanceMatrixResult = explicitGraph.makeDistanceMatrix(configuration.partiality);
    if(!distanceMatrixResult) {
      Log::log(Log::Level::Warning) << "Failure in distance matrix construction.\n";
      failures += 1;
      loosenBounds = true;
      if(detail::exceededFailureRatio(failures, numConformers, configuration.failureRatio)) {
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

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      static_cast<Eigen::MatrixXd>(
        distanceBounds.access().cwiseProduct(distanceBounds.access())
      )
    );

    /* Our embedded coordinates are four dimensional. Now we want to make sure
     * that all chiral constraints are correct, allowing the structure to expand
     * into the fourth spatial dimension if necessary to allow inversion.
     *
     * This stage of refinement is only needed if not all chirality constraints
     * are already correct (or there are none).
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
        configuration.refinementStepLimit,
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

      // Handle inversion failure (hit step limit)
      if(
        inversionStopStrategy.iterations > configuration.refinementStepLimit
        || errfDetail::proportionChiralityConstraintsCorrectSign(
          DGData.chiralityConstraints,
          dlibPositions
        ) < 1.0
      ) {
        Log::log(Log::Level::Warning)
          << "[" << currentStructureNumber << "]: "
          << "First stage of refinement fails. Loosening factor was "
          << detail::looseningFactor(failures, numConformers, configuration.failureRatio)
          <<  "\n";
        failures += 1;
        loosenBounds = true;
        if(detail::exceededFailureRatio(failures, numConformers, configuration.failureRatio)) {
          Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
          return refinementList;
        }
        continue; // this triggers a new structure to be generated
      }
    }

    /* Set up the second stage of refinement where we compress out the fourth
     * dimension that we allowed expansion into to invert the chiralities.
     */

    errfValue<true> refinementValueFunctor {
      squaredBounds,
      DGData.chiralityConstraints
    };

    dlibAdaptors::debugIterationOrGradientNormStopStrategy refinementStopStrategy {
      configuration.refinementStepLimit,
      configuration.refinementGradientTarget,
      refinementSteps,
      refinementValueFunctor
    };

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

    bool reachedMaxIterations = refinementStopStrategy.iterations >= configuration.refinementStepLimit;
    bool notAllChiralitiesCorrect = errfDetail::proportionChiralityConstraintsCorrectSign(
      DGData.chiralityConstraints,
      dlibPositions
    ) < 1;
    bool structureAcceptable = errfDetail::finalStructureAcceptable(
      distanceBounds,
      DGData.chiralityConstraints,
      dlibPositions
    );

    if(Log::particulars.count(Log::Particulars::DGFinalErrorContributions) > 0) {
      errfDetail::explainFinalContributions(
        distanceBounds,
        DGData.chiralityConstraints,
        dlibPositions
      );
    }

    RefinementData refinementData;
    refinementData.steps = std::move(refinementSteps);
    refinementData.constraints = DGData.chiralityConstraints;
    refinementData.looseningFactor = detail::looseningFactor(failures, numConformers, configuration.failureRatio);
    refinementData.isFailure = (reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable);
    refinementData.spatialModelGraphviz = spatialModelGraphviz;

    refinementList.push_back(
      std::move(refinementData)
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      Log::log(Log::Level::Warning)
        << "[" << currentStructureNumber << "]: "
        << "Second stage of refinement fails. Loosening factor was "
        << detail::looseningFactor(failures, numConformers, configuration.failureRatio)
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
      if(detail::exceededFailureRatio(failures, numConformers, configuration.failureRatio)) {
        Log::log(Log::Level::Warning) << "Exceeded failure ratio in debug DG.\n";
        return refinementList;
      }
    }
  }

  return refinementList;
}


// Non-debug version of DG
outcome::result<
  std::vector<AngstromWrapper>
> run(
  const Molecule& molecule,
  const unsigned numConformers,
  const Configuration& configuration
) {
  if(molecule.stereopermutators().hasZeroAssignmentStereopermutators()) {
    return DGError::ZeroAssignmentStereopermutators;
  }

  MoleculeDGInformation DGData;

  /* In case the molecule has unassigned stereopermutators, we need to randomly
   * assign them for each conformer generated prior to generating the distance
   * bounds matrix
   */
  bool regenerateEachStep = molecule.stereopermutators().hasUnassignedStereopermutators();

  /* If there are no unassigned stereocenters, DGData does not have to be
   * regenerated. If there are no failures, DGData stays the same for the rest
   * of the procedure.
   */
  if(!regenerateEachStep) {
    DGData = gatherDGInformation(molecule, configuration);
  }

  /* Allow up to failures / numConformers = 2 of errors.
   * function returns an error code.
   */
  unsigned failures = 0;
  // Track whether a failure should lead to a bounds loosening
  bool loosenBounds = false;
  std::vector<AngstromWrapper> ensemble;
  ensemble.reserve(numConformers);

  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numConformers;
    currentStructureNumber += 1
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = detail::narrow(molecule);

      if(moleculeCopy.stereopermutators().hasZeroAssignmentStereopermutators()) {
        return DGError::ZeroAssignmentStereopermutators;
      }

      // Fetch the DG data from the molecule with no unassigned stereopermutators
      DGData = gatherDGInformation(
        moleculeCopy,
        configuration,
        detail::looseningFactor(failures, numConformers, configuration.failureRatio)
      );
    } else if(loosenBounds) {
      DGData = gatherDGInformation(
        molecule,
        configuration,
        detail::looseningFactor(failures, numConformers, configuration.failureRatio)
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

    auto distanceMatrixResult = explicitGraph.makeDistanceMatrix(configuration.partiality);
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
        configuration.refinementStepLimit
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
        inversionStopStrategy.iterations >= configuration.refinementStepLimit
        || errfDetail::proportionChiralityConstraintsCorrectSign(
          DGData.chiralityConstraints,
          dlibPositions
        ) < 1.0
      ) { // Failure to invert
        failures += 1;
        loosenBounds = true;
        if(detail::exceededFailureRatio(failures, numConformers, configuration.failureRatio)) {
          return DGError::TooManyFailures;
        }
        continue; // this triggers a new structure to be generated
      }
    }

    dlibAdaptors::iterationOrGradientNormStopStrategy refinementStopStrategy {
      configuration.refinementStepLimit,
      configuration.refinementGradientTarget
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

    bool reachedMaxIterations = refinementStopStrategy.iterations >= configuration.refinementStepLimit;
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

      if(detail::exceededFailureRatio(failures, numConformers, configuration.failureRatio)) {
        return DGError::TooManyFailures;
      }
    } else {
      if(!configuration.fixedPositions.empty()) {
        auto positionMatrix = detail::fitAndSetFixedPositions(dlibPositions, configuration);

        ensemble.emplace_back(
          detail::convertToAngstromWrapper(positionMatrix)
        );
      } else {
        ensemble.emplace_back(
          detail::convertToAngstromWrapper(dlibPositions)
        );
      }
    }
  }

  return ensemble;
}

} // namespace DistanceGeometry

} // namespace molassembler
