/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include <Eigen/Dense>
#include "Utils/Constants.h"
#include "Utils/Typenames.h"

#include "molassembler/DistanceGeometry/EigenRefinement.h"
#include "molassembler/DistanceGeometry/Error.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/RefinementMeta.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "Utils/Math/QuaternionFit.h"

#include "temple/Optimization/LBFGS.h"
#include "temple/Optionals.h"
#include "temple/Random.h"

#include <iostream>

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {
namespace detail {

Eigen::MatrixXd gather(const Eigen::VectorXd& vectorizedPositions) {
  constexpr unsigned dimensionality = 4;
  const unsigned N = vectorizedPositions.size() / dimensionality;
  Eigen::MatrixXd positionMatrix(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    positionMatrix.row(i) = vectorizedPositions.template segment<3>(dimensionality * i);
  }
  return positionMatrix;
}

AngstromWrapper convertToAngstromWrapper(const Eigen::MatrixXd& positions) {
  assert(positions.cols() == 3);
  const unsigned N = positions.rows();
  AngstromWrapper angstromWrapper {N};
  for(unsigned i = 0; i < N; ++i) {
    angstromWrapper.positions.row(i) = Scine::Utils::Position {
      positions.row(i).transpose()
    };
  }
  return angstromWrapper;
}

Eigen::MatrixXd fitAndSetFixedPositions(
  const Eigen::MatrixXd& positions,
  const Configuration& configuration
) {
  /* Fixed positions postprocessing:
   * - Rotate and translate the generated coordinates towards the fixed
   *   positions indicated for each.
   *
   * Maybe?
   * - Assuming the fit isn't absolutely exact, overwrite the existing
   *   positions with the fixed ones
   */
  /* Construct a reference matrix from the fixed positions (these are still in
   * bohr!) and a weights vector
   */
  assert(positions.cols() == 3);
  const unsigned N = positions.rows();
  Eigen::MatrixXd referenceMatrix = Eigen::MatrixXd::Zero(N, 3);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
  for(const auto& indexPositionPair : configuration.fixedPositions) {
    referenceMatrix.row(indexPositionPair.first) = indexPositionPair.second.transpose();
    weights(indexPositionPair.first) = 1;
  }

  referenceMatrix *= Scine::Utils::Constants::angstrom_per_bohr;

  // Perform the QuaternionFit
  Scine::Utils::QuaternionFit fit(referenceMatrix, positions, weights);

  return fit.getFittedData();
}

Molecule narrow(Molecule molecule, random::Engine& engine) {
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
        engine
      );

      molecule.assignStereopermutatorRandomly(
        unassignedAtomStereopermutators.at(choice),
        engine
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
        engine
      );

      molecule.assignStereopermutatorRandomly(
        unassignedBondStereopermutators.at(choice),
        engine
      );
    }

  } while(stereopermutatorList.hasUnassignedStereopermutators());

  return molecule;
}

template<typename EigenRefinementType>
struct InversionOrIterLimitStop {
  const unsigned iterLimit;
  const EigenRefinementType& refinementFunctorReference;


  using VectorType = typename EigenRefinementType::VectorType;
  using FloatType = typename EigenRefinementType::FloatingPointType;

  InversionOrIterLimitStop(
    const unsigned passIter,
    const EigenRefinementType& functor
  ) : iterLimit(passIter),
      refinementFunctorReference(functor)
  {}

  template<typename StepValues>
  bool shouldContinue(unsigned iteration, const StepValues& /* step */) {
    return (
      iteration < iterLimit
      && refinementFunctorReference.proportionChiralConstraintsCorrectSign < 1.0
    );
  }
};

template<typename FloatType>
struct GradientOrIterLimitStop {
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

  template<typename StepValues>
  bool shouldContinue(unsigned iteration, const StepValues& step) {
    return (
      iteration < iterLimit
      && step.gradients.current.template cast<double>().norm() > gradNorm
    );
  }

  unsigned iterLimit = 10000;
  double gradNorm = 1e-5;
};

} // namespace detail

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration
) {
  // Generate a spatial model from the molecular graph and stereopermutators
  SpatialModel spatialModel {molecule, configuration};

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makePairwiseBounds();
  data.chiralConstraints = spatialModel.getChiralConstraints();
  data.dihedralConstraints = spatialModel.getDihedralConstraints();

  return data;
}

outcome::result<AngstromWrapper> refine(
  Eigen::MatrixXd embeddedPositions,
  const DistanceBoundsMatrix& distanceBounds,
  const Configuration& configuration,
  const std::shared_ptr<MoleculeDGInformation>& DGDataPtr
) {
  /* Refinement problem compile-time settings
   * - Dimensionality four is needed to ensure chiral constraints invert
   *   nicely
   * - FloatType double is helpful for refinement stability, and float
   *   doesn't affect speed
   * - Using the alternative SIMD implementations of the refinement problems
   *   hardly affects speed at all, in fact, it commonly worsens it.
   */
  constexpr unsigned dimensionality = 4;
  using FloatType = double;
  constexpr bool SIMD = false;

  using FullRefinementType = EigenRefinementProblem<dimensionality, FloatType, SIMD>;
  using VectorType = typename FullRefinementType::VectorType;

  // Vectorize positions
  VectorType transformedPositions = Eigen::Map<Eigen::VectorXd>(
    embeddedPositions.data(),
    embeddedPositions.cols() * embeddedPositions.rows()
  ).template cast<FloatType>().eval();

  const unsigned N = transformedPositions.size() / dimensionality;

  const auto squaredBounds = static_cast<Eigen::MatrixXd>(
    distanceBounds.access().cwiseProduct(distanceBounds.access())
  );

  FullRefinementType refinementFunctor {
    squaredBounds,
    DGDataPtr->chiralConstraints,
    DGDataPtr->dihedralConstraints
  };

  /* If a count of chiral constraints reveals that more than half are
   * incorrect, we can invert the structure (by multiplying e.g. all y
   * coordinates with -1) and then have more than half of chirality
   * constraints correct! In the count, chiral constraints with a target
   * value of zero are not considered (this would skew the count as those
   * chiral constraints should not have to pass an energetic maximum to
   * converge properly as opposed to tetrahedra with volume).
   */
  double initiallyCorrectChiralConstraints = refinementFunctor.calculateProportionChiralConstraintsCorrectSign(transformedPositions);
  if(initiallyCorrectChiralConstraints < 0.5) {
    // Invert y coordinates
    for(unsigned i = 0; i < N; ++i) {
      transformedPositions(dimensionality * i + 1) *= -1;
    }

    initiallyCorrectChiralConstraints = 1 - initiallyCorrectChiralConstraints;
  }

  /* Refinement without penalty on fourth dimension only necessary if not all
   * chiral centers are correct. Of course, for molecules without chiral
   * centers at all, this stage is unnecessary
   */
  unsigned firstStageIterations = 0;
  if(initiallyCorrectChiralConstraints < 1) {
    detail::InversionOrIterLimitStop<FullRefinementType> inversionChecker {
      configuration.refinementStepLimit,
      refinementFunctor
    };

    temple::LBFGS<FloatType, 32> optimizer;

    try {
      auto result = optimizer.minimize(
        transformedPositions,
        refinementFunctor,
        inversionChecker
      );
      firstStageIterations = result.iterations;
    } catch(std::runtime_error& e) {
      return DGError::RefinementException;
    }

    if(firstStageIterations >= configuration.refinementStepLimit) {
      return DGError::RefinementMaxIterationsReached;
    }

    if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1.0) {
      return DGError::RefinedChiralsWrong;
    }
  }

  /* Set up the second stage of refinement where we compress out the fourth
   * dimension that we allowed expansion into to invert the chiralities.
   */
  refinementFunctor.compressFourthDimension = true;

  unsigned secondStageIterations = 0;
  detail::GradientOrIterLimitStop<FloatType> gradientChecker;
  gradientChecker.gradNorm = 1e-3;
  gradientChecker.iterLimit = configuration.refinementStepLimit - firstStageIterations;

  try {
    temple::LBFGS<FloatType, 32> optimizer;

    auto result = optimizer.minimize(
      transformedPositions,
      refinementFunctor,
      gradientChecker
    );
    secondStageIterations = result.iterations;
  } catch(std::out_of_range& e) {
    return DGError::RefinementException;
  }

  // Max iterations reached
  if(secondStageIterations >= gradientChecker.iterLimit) {
    return DGError::RefinementMaxIterationsReached;
  }

  // Not all chiral constraints have the right sign
  if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1) {
    return DGError::RefinedChiralsWrong;
  }

  /* Add dihedral terms and refine again */
  unsigned thirdStageIterations = 0;
  gradientChecker = detail::GradientOrIterLimitStop<FloatType> {};
  gradientChecker.gradNorm = 1e-3;
  gradientChecker.iterLimit = (
    configuration.refinementStepLimit
    - firstStageIterations
    - secondStageIterations
  );

  refinementFunctor.dihedralTerms = true;

  try {
    temple::LBFGS<FloatType, 32> optimizer;

    auto result = optimizer.minimize(
      transformedPositions,
      refinementFunctor,
      gradientChecker
    );
    thirdStageIterations = result.iterations;
  } catch(std::out_of_range& e) {
    return DGError::RefinementException;
  }

  if(thirdStageIterations >= gradientChecker.iterLimit) {
    return DGError::RefinementMaxIterationsReached;
  }

  // Structure inacceptable
  if(
    !finalStructureAcceptable(
      refinementFunctor,
      distanceBounds,
      transformedPositions
    )
  ) {
    return DGError::RefinedStructureInacceptable;
  }

  auto gatheredPositions = detail::gather(transformedPositions);

  if(!configuration.fixedPositions.empty()) {
    return detail::convertToAngstromWrapper(
      detail::fitAndSetFixedPositions(gatheredPositions, configuration)
    );
  }

  return detail::convertToAngstromWrapper(gatheredPositions);
}

outcome::result<AngstromWrapper> generateConformer(
  const Molecule& molecule,
  const Configuration& configuration,
  std::shared_ptr<MoleculeDGInformation>& DGDataPtr,
  bool regenerateDGDataEachStep,
  random::Engine& engine
) {
  if(regenerateDGDataEachStep) {
    auto moleculeCopy = detail::narrow(molecule, engine);

    if(moleculeCopy.stereopermutators().hasZeroAssignmentStereopermutators()) {
      return DGError::ZeroAssignmentStereopermutators;
    }

    DGDataPtr = std::make_shared<MoleculeDGInformation>(
      gatherDGInformation(moleculeCopy, configuration)
    );
  }

  ExplicitGraph explicitGraph {
    molecule.graph().inner(),
    DGDataPtr->bounds
  };

  // Get distance bounds matrix from the graph
  auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
  if(!distanceBoundsResult) {
    return distanceBoundsResult.as_failure();
  }

  DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

  /* There should be no need to smooth the distance bounds, because the graph
   * type ought to create them within the triangle inequality bounds:
   */
  assert(distanceBounds.boundInconsistencies() == 0);

  // Generate a distances matrix from the graph
  auto distanceMatrixResult = explicitGraph.makeDistanceMatrix(
    engine,
    configuration.partiality
  );
  if(!distanceMatrixResult) {
    return distanceMatrixResult.as_failure();
  }

  // Make a metric matrix from the distances matrix
  MetricMatrix metric(
    std::move(distanceMatrixResult.value())
  );

  // Get a position matrix by embedding the metric matrix
  auto embeddedPositions = metric.embed();

  /* Refinement */
  return refine(
    std::move(embeddedPositions),
    distanceBounds,
    configuration,
    DGDataPtr
  );
}

std::vector<
  outcome::result<AngstromWrapper>
> run(
  const Molecule& molecule,
  const unsigned numConformers,
  const Configuration& configuration,
  const boost::optional<unsigned> seedOption
) {
  using ReturnType = std::vector<
    outcome::result<AngstromWrapper>
  >;

  // In case there are zero assignment stereopermutators, we give up immediately
  if(molecule.stereopermutators().hasZeroAssignmentStereopermutators()) {
    return ReturnType(numConformers, DGError::ZeroAssignmentStereopermutators);
  }

#ifdef _OPENMP
  /* Ensure the molecule's mutable properties are already generated so none are
   * generated on threaded const-access.
   */
  molecule.graph().inner().populateProperties();
#endif

  /* In case the molecule has unassigned stereopermutators, we need to randomly
   * assign them for each conformer generated prior to generating the distance
   * bounds matrix. If not, then modelling data can be kept across all
   * conformer generation runs since no randomness has entered the equation.
   */
  auto DGDataPtr = std::make_shared<MoleculeDGInformation>();
  bool regenerateEachStep = molecule.stereopermutators().hasUnassignedStereopermutators();
  if(!regenerateEachStep) {
    *DGDataPtr = gatherDGInformation(molecule, configuration);
  }

  ReturnType results(numConformers, static_cast<DGError>(0));

  auto engineOption = temple::optionals::map(
    seedOption,
    [](unsigned seed) -> random::Engine {
      random::Engine engine;
      engine.seed(seed);
      return engine;
    }
  );

  std::reference_wrapper<random::Engine> backgroundEngineWrapper = randomnessEngine();
  if(engineOption) {
    backgroundEngineWrapper = engineOption.value();
  }
  random::Engine& backgroundEngine = backgroundEngineWrapper.get();

#pragma omp parallel
  {
#ifdef _OPENMP
    /* We have to distribute pseudo-randomness into each thread reproducibly
     * and want to avoid having to guard the global PRNG against access from
     * multiple threads, so we provide each thread its own Engine.
     */
    const unsigned nThreads = omp_get_num_threads();
    std::vector<random::Engine> randomnessEngines(nThreads);
#endif

    /* Each thread has its own DGDataPtr, for the following reason: If we do
     * not need to regenerate the SpatialModel data, then having all threads
     * share access the underlying data to generate conformers is fine. If,
     * otherwise, we do need to regenerate the SpatialModel data for each
     * conformer, then each thread will reset its pointer to its self-generated
     * SpatialModel data, creating thread-private state.
     */
#pragma omp for firstprivate(DGDataPtr)
    for(unsigned i = 0; i < numConformers; ++i) {
      // Get thread-specific randomness engine reference
#ifdef _OPENMP
      random::Engine& engine = randomnessEngines.at(
        omp_get_thread_num()
      );

      // Re-seed the thread-local PRNG engine for each conformer
#pragma omp critical
      engine.seed(backgroundEngine());

#else
      random::Engine engine;
      engine.seed(backgroundEngine());
#endif

      /* We have to handle any and all exceptions here bceause this is a
       * parallel environment
       */
      try {
        // Generate the conformer
        auto conformerResult = generateConformer(
          molecule,
          configuration,
          DGDataPtr,
          regenerateEachStep,
          engine
        );

        results.at(i) = std::move(conformerResult);
      } catch(std::exception& e) {
#pragma omp critical(collectConformer)
        {
          std::cerr << "WARNING: Uncaught exception in conformer generation: " << e.what() << "\n";
        }
      } // end catch
    } // end pragma omp for private(DGDataPtr)
  } // end pragma omp parallel

  return results;
}

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine
