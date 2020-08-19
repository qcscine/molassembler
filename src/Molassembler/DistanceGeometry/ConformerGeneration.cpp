/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/DistanceGeometry/ConformerGeneration.h"

#include <Eigen/Dense>
#include "Utils/Constants.h"
#include "Utils/Typenames.h"

#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/DistanceGeometry/EigenRefinement.h"
#include "Molassembler/DistanceGeometry/Error.h"
#include "Molassembler/DistanceGeometry/ExplicitBoundsGraph.h"
#include "Molassembler/DistanceGeometry/MetricMatrix.h"
#include "Molassembler/DistanceGeometry/RefinementMeta.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Utils/Math/QuaternionFit.h"

#include "Molassembler/Temple/Optimization/Lbfgs.h"
#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/Random.h"

#include <iostream>

namespace Scine {
namespace Molassembler {
namespace DistanceGeometry {
namespace Detail {

Eigen::MatrixXd gather(const Eigen::VectorXd& vectorizedPositions) {
  constexpr unsigned dimensionality = 4;
  const unsigned N = vectorizedPositions.size() / dimensionality;
  Eigen::MatrixXd positionMatrix(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    positionMatrix.row(i) = vectorizedPositions.template segment<3>(dimensionality * i);
  }
  return positionMatrix;
}

AngstromPositions convertToAngstromPositions(const Eigen::MatrixXd& positions) {
  assert(positions.cols() == 3);
  const unsigned N = positions.rows();
  AngstromPositions angstromWrapper {N};
  for(unsigned i = 0; i < N; ++i) {
    angstromWrapper.positions.row(i) = Utils::Position {
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

  referenceMatrix *= Utils::Constants::angstrom_per_bohr;

  // Perform the QuaternionFit
  Utils::QuaternionFit fit(referenceMatrix, positions, weights);

  return fit.getFittedData();
}

Molecule narrow(Molecule molecule, Random::Engine& engine) {
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
          atomStereopermutator.placement()
        );
      }
    }

    if(!unassignedAtomStereopermutators.empty()) {
      molecule.assignStereopermutatorRandomly(
        Temple::Random::pick(unassignedAtomStereopermutators, engine),
        engine
      );

      // Re-check the loop condition
      continue;
    }

    std::vector<BondIndex> unassignedBondStereopermutators;

    for(const auto& bondStereopermutator : stereopermutatorList.bondStereopermutators()) {
      if(!bondStereopermutator.assigned()) {
        unassignedBondStereopermutators.push_back(
          bondStereopermutator.placement()
        );
      }
    }

    if(!unassignedBondStereopermutators.empty()) {
      molecule.assignStereopermutatorRandomly(
        Temple::Random::pick(unassignedBondStereopermutators, engine),
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

} // namespace Detail

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

outcome::result<AngstromPositions> refine(
  Eigen::MatrixXd embeddedPositions,
  const DistanceBoundsMatrix& distanceBounds,
  const Configuration& configuration,
  const std::shared_ptr<MoleculeDGInformation>& DgDataPtr
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
    DgDataPtr->chiralConstraints,
    DgDataPtr->dihedralConstraints
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
    Detail::InversionOrIterLimitStop<FullRefinementType> inversionChecker {
      configuration.refinementStepLimit,
      refinementFunctor
    };

    Temple::Lbfgs<FloatType, 32> optimizer;

    try {
      auto result = optimizer.minimize(
        transformedPositions,
        refinementFunctor,
        inversionChecker
      );
      firstStageIterations = result.iterations;
    } catch(std::runtime_error& e) {
      return DgError::RefinementException;
    }

    if(firstStageIterations >= configuration.refinementStepLimit) {
      return DgError::RefinementMaxIterationsReached;
    }

    if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1.0) {
      return DgError::RefinedChiralsWrong;
    }
  }

  /* Set up the second stage of refinement where we compress out the fourth
   * dimension that we allowed expansion into to invert the chiralities.
   */
  refinementFunctor.compressFourthDimension = true;

  unsigned secondStageIterations = 0;
  Detail::GradientOrIterLimitStop<FloatType> gradientChecker;
  gradientChecker.gradNorm = 1e-3;
  gradientChecker.iterLimit = configuration.refinementStepLimit - firstStageIterations;

  try {
    Temple::Lbfgs<FloatType, 32> optimizer;

    auto result = optimizer.minimize(
      transformedPositions,
      refinementFunctor,
      gradientChecker
    );
    secondStageIterations = result.iterations;
  } catch(std::out_of_range& e) {
    return DgError::RefinementException;
  }

  // Max iterations reached
  if(secondStageIterations >= gradientChecker.iterLimit) {
    return DgError::RefinementMaxIterationsReached;
  }

  // Not all chiral constraints have the right sign
  if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1) {
    return DgError::RefinedChiralsWrong;
  }

  /* Add dihedral terms and refine again */
  unsigned thirdStageIterations = 0;
  gradientChecker = Detail::GradientOrIterLimitStop<FloatType> {};
  gradientChecker.gradNorm = 1e-3;
  gradientChecker.iterLimit = (
    configuration.refinementStepLimit
    - firstStageIterations
    - secondStageIterations
  );

  refinementFunctor.dihedralTerms = true;

  try {
    Temple::Lbfgs<FloatType, 32> optimizer;

    auto result = optimizer.minimize(
      transformedPositions,
      refinementFunctor,
      gradientChecker
    );
    thirdStageIterations = result.iterations;
  } catch(std::out_of_range& e) {
    return DgError::RefinementException;
  }

  if(thirdStageIterations >= gradientChecker.iterLimit) {
    return DgError::RefinementMaxIterationsReached;
  }

  // Structure inacceptable
  if(
    !finalStructureAcceptable(
      refinementFunctor,
      distanceBounds,
      transformedPositions
    )
  ) {
    return DgError::RefinedStructureInacceptable;
  }

  auto gatheredPositions = Detail::gather(transformedPositions);

  if(!configuration.fixedPositions.empty()) {
    return Detail::convertToAngstromPositions(
      Detail::fitAndSetFixedPositions(gatheredPositions, configuration)
    );
  }

  return Detail::convertToAngstromPositions(gatheredPositions);
}

outcome::result<AngstromPositions> generateConformer(
  const Molecule& molecule,
  const Configuration& configuration,
  std::shared_ptr<MoleculeDGInformation>& DgDataPtr,
  bool regenerateDGDataEachStep,
  Random::Engine& engine
) {
  if(regenerateDGDataEachStep) {
    auto moleculeCopy = Detail::narrow(molecule, engine);

    if(moleculeCopy.stereopermutators().hasZeroAssignmentStereopermutators()) {
      return DgError::ZeroAssignmentStereopermutators;
    }

    DgDataPtr = std::make_shared<MoleculeDGInformation>(
      gatherDGInformation(moleculeCopy, configuration)
    );
  }

  ExplicitBoundsGraph explicitGraph {
    molecule.graph().inner(),
    DgDataPtr->bounds
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
    DgDataPtr
  );
}

std::vector<
  outcome::result<AngstromPositions>
> run(
  const Molecule& molecule,
  const unsigned numConformers,
  const Configuration& configuration,
  const boost::optional<unsigned> seedOption
) {
  using ReturnType = std::vector<
    outcome::result<AngstromPositions>
  >;

  // In case there are zero assignment stereopermutators, we give up immediately
  if(molecule.stereopermutators().hasZeroAssignmentStereopermutators()) {
    return ReturnType(numConformers, DgError::ZeroAssignmentStereopermutators);
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
  auto DgDataPtr = std::make_shared<MoleculeDGInformation>();
  bool regenerateEachStep = molecule.stereopermutators().hasUnassignedStereopermutators();
  if(!regenerateEachStep) {
    *DgDataPtr = gatherDGInformation(molecule, configuration);
  }

  ReturnType results(numConformers, static_cast<DgError>(0));

  /* If a seed is supplied, the global prng state is not to be advanced.
   * We create a random engine from the seed here if a seed is supplied.
   */
  auto engineOption = Temple::Optionals::map(
    seedOption,
    [](unsigned seed) { return Random::Engine(seed); }
  );

  /* Now we need something referencing either our new engine or the global prng
   * engine.
   */
  std::reference_wrapper<Random::Engine> backgroundEngineWrapper = randomnessEngine();
  if(engineOption) {
    backgroundEngineWrapper = engineOption.value();
  }
  Random::Engine& backgroundEngine = backgroundEngineWrapper.get();

  /* We have to distribute pseudo-randomness into each thread reproducibly
   * and want to avoid having to guard the global PRNG against access from
   * multiple threads, so we provide each thread its own Engine and pre-generate
   * each conformer's individual seed in the sequential section.
   */
#ifdef _OPENMP
  const unsigned nThreads = omp_get_max_threads();
#else
  const unsigned nThreads = 1;
#endif

  std::vector<Random::Engine> randomnessEngines(nThreads);
  const auto seeds = Temple::Random::getN<int>(
    0,
    std::numeric_limits<int>::max(),
    numConformers,
    backgroundEngine
  );

  /* Each thread has its own DgDataPtr, for the following reason: If we do
   * not need to regenerate the SpatialModel data, then having all threads
   * share access the underlying data to generate conformers is fine. If,
   * otherwise, we do need to regenerate the SpatialModel data for each
   * conformer, then each thread will reset its pointer to its self-generated
   * SpatialModel data, creating thread-private state.
   */
#pragma omp parallel for firstprivate(DgDataPtr) schedule(dynamic)
  for(unsigned i = 0; i < numConformers; ++i) {
    // Get thread-specific randomness engine reference
#ifdef _OPENMP
    Random::Engine& engine = randomnessEngines.at(
      omp_get_thread_num()
    );
#else
    Random::Engine& engine = randomnessEngines.front();
#endif

    // Re-seed the thread-local PRNG engine for each conformer
    engine.seed(seeds.at(i));

    /* We have to handle any and all exceptions here bceause this is a parallel
     * environment and exceptions are not propagated anywhere
     */
    try {
      // Generate the conformer
      auto conformerResult = generateConformer(
        molecule,
        configuration,
        DgDataPtr,
        regenerateEachStep,
        engine
      );

      results.at(i) = std::move(conformerResult);
    } catch(std::exception& e) {
#pragma omp critical(outputWarning)
      {
        std::cerr << "WARNING: Uncaught exception in conformer generation: " << e.what() << "\n";
      }
      results.at(i) = DgError::UnknownException;
    } // end catch
  } // end pragma omp for private(DgDataPtr)

  return results;
}

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine
