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
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "Utils/Math/QuaternionFit.h"

#include "temple/Optimization/LBFGS.h"
#include "temple/Random.h"

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
  const Configuration& configuration,
  random::Engine& engine
) {
  // Generate a spatial model from the molecular graph and stereopermutators
  SpatialModel spatialModel {molecule, configuration, engine};

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makePairwiseBounds();
  data.chiralConstraints = spatialModel.getChiralConstraints();
  data.dihedralConstraints = spatialModel.getDihedralConstraints();

  return data;
}

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration,
  std::string& spatialModelGraphvizString
) {
  // Generate a spatial model from the molecular graph and stereopermutators
  SpatialModel spatialModel {molecule, configuration, randomnessEngine()};
  spatialModelGraphvizString = spatialModel.dumpGraphviz();

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makePairwiseBounds();
  data.chiralConstraints = spatialModel.getChiralConstraints();
  data.dihedralConstraints = spatialModel.getDihedralConstraints();

  return data;
}

// Debug version
std::list<RefinementData> debugRefinement(
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
    DGData = gatherDGInformation(molecule, configuration, spatialModelGraphviz);
  }

  unsigned failures = 0;
  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber < numConformers;
    ++currentStructureNumber
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = detail::narrow(molecule, randomnessEngine());

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
        spatialModelGraphviz
      );
    }

    std::list<RefinementStepData> refinementSteps;

    ExplicitGraph explicitGraph {
      molecule.graph().inner(),
      DGData.bounds
    };

    auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
    if(!distanceBoundsResult) {
      Log::log(Log::Level::Warning) << "Failure in distance bounds matrix construction: "
        << distanceBoundsResult.error().message() << "\n";
      failures += 1;

      if(regenerateEachStep) {
        auto moleculeCopy = detail::narrow(molecule, randomnessEngine());

        SpatialModel model {moleculeCopy, configuration, randomnessEngine()};
        model.writeGraphviz("DG-failure-spatial-model-" + std::to_string(currentStructureNumber) + ".dot");
      } else {
        SpatialModel model {molecule, configuration, randomnessEngine()};
        model.writeGraphviz("DG-failure-spatial-model-" + std::to_string(currentStructureNumber) + ".dot");
      }

      continue;
    }

    DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

    /* No need to smooth the distance bounds, ExplicitGraph creates it
     * so that the triangle inequalities are fulfilled
     */

    assert(distanceBounds.boundInconsistencies() == 0);

    auto distanceMatrixResult = explicitGraph.makeDistanceMatrix(
      randomnessEngine(),
      configuration.partiality
    );
    if(!distanceMatrixResult) {
      Log::log(Log::Level::Warning) << "Failure in distance matrix construction.\n";
      failures += 1;
    }

    // Make a metric matrix from a generated distances matrix
    MetricMatrix metric(
      std::move(distanceMatrixResult.value())
    );

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

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
      DGData.chiralConstraints,
      DGData.dihedralConstraints
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

    struct Observer {
      const Molecule& mol;
      FullRefinementType& functor;
      std::list<RefinementStepData>& steps;

      Observer(const Molecule& passMolecule, FullRefinementType& passFunctor, std::list<RefinementStepData>& passSteps)
        : mol(passMolecule), functor(passFunctor), steps(passSteps) {}

      using VectorType = typename FullRefinementType::VectorType;

      void operator() (const VectorType& positions) {
        FloatType distanceError = 0, chiralError = 0,
                  dihedralError = 0, fourthDimensionError = 0;

        VectorType gradient;
        gradient.resize(positions.size());
        gradient.setZero();

        // Call functions individually to get separate error values
        functor.distanceContributions(positions, distanceError, gradient);
        functor.chiralContributions(positions, chiralError, gradient);
        functor.dihedralContributions(positions, dihedralError, gradient);
        functor.fourthDimensionContributions(positions, fourthDimensionError, gradient);

        steps.emplace_back(
          positions.template cast<double>().eval(),
          distanceError,
          chiralError,
          dihedralError,
          fourthDimensionError,
          gradient.template cast<double>().eval(),
          functor.proportionChiralConstraintsCorrectSign,
          functor.compressFourthDimension
        );
      }
    };

    Observer observer(molecule, refinementFunctor, refinementSteps);

    /* Our embedded coordinates are (dimensionality) dimensional. Now we want
     * to make sure that all chiral constraints are correct, allowing the
     * structure to expand into the fourth spatial dimension if necessary to
     * allow inversion.
     *
     * This stage of refinement is only needed if not all chiral constraints
     * are already correct (or there are none).
     */
    unsigned firstStageIterations = 0;
    if(initiallyCorrectChiralConstraints < 1.0) {
      detail::InversionOrIterLimitStop<FullRefinementType> inversionChecker {
        configuration.refinementStepLimit,
        refinementFunctor
      };

      temple::LBFGS<FloatType, 32> optimizer;

      try {
        auto result = optimizer.minimize(
          transformedPositions,
          refinementFunctor,
          inversionChecker,
          observer
        );
        firstStageIterations = result.iterations;
      } catch(std::runtime_error& e) {
        Log::log(Log::Level::Warning)
          << "Non-finite contributions to dihedral error function gradient.\n";
        failures += 1;
        continue;
      }

      // Handle inversion failure (hit step limit)
      if(
        firstStageIterations >= configuration.refinementStepLimit
        || refinementFunctor.proportionChiralConstraintsCorrectSign < 1.0
      ) {
        Log::log(Log::Level::Warning)
          << "[" << currentStructureNumber << "]: "
          << "First stage of refinement fails. Loosening factor was "
          << configuration.spatialModelLoosening
          <<  "\n";
        failures += 1;
        continue; // this triggers a new structure to be generated
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
        gradientChecker,
        observer
      );
      secondStageIterations = result.iterations;
    } catch(std::out_of_range& e) {
      Log::log(Log::Level::Warning)
        << "Non-finite contributions to dihedral error function gradient.\n";
      failures += 1;
      continue;
    }

    if(secondStageIterations >= gradientChecker.iterLimit) {
        Log::log(Log::Level::Warning)
          << "[" << currentStructureNumber << "]: "
          << "Second stage of refinement fails!\n";
        failures += 1;

        // Collect refinement data
        RefinementData refinementData;
        refinementData.steps = std::move(refinementSteps);
        refinementData.constraints = DGData.chiralConstraints;
        refinementData.looseningFactor = configuration.spatialModelLoosening;
        refinementData.isFailure = true;
        refinementData.spatialModelGraphviz = spatialModelGraphviz;

        refinementList.push_back(
          std::move(refinementData)
        );

        if(Log::particulars.count(Log::Particulars::DGFinalErrorContributions) > 0) {
          explainFinalContributions(
            refinementFunctor,
            distanceBounds,
            transformedPositions
          );
        }

        continue; // this triggers a new structure to be generated
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
        gradientChecker,
        observer
      );
      thirdStageIterations = result.iterations;
    } catch(std::out_of_range& e) {
      Log::log(Log::Level::Warning)
        << "Non-finite contributions to dihedral error function gradient.\n";
      failures += 1;
      continue;
    }

    bool reachedMaxIterations = thirdStageIterations >= gradientChecker.iterLimit;
    bool notAllChiralitiesCorrect = refinementFunctor.proportionChiralConstraintsCorrectSign < 1;
    bool structureAcceptable = finalStructureAcceptable(
      refinementFunctor,
      distanceBounds,
      transformedPositions
    );

    if(Log::particulars.count(Log::Particulars::DGFinalErrorContributions) > 0) {
      explainFinalContributions(
        refinementFunctor,
        distanceBounds,
        transformedPositions
      );
    }

    RefinementData refinementData;
    refinementData.steps = std::move(refinementSteps);
    refinementData.constraints = DGData.chiralConstraints;
    refinementData.looseningFactor = configuration.spatialModelLoosening;
    refinementData.isFailure = (reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable);
    refinementData.spatialModelGraphviz = spatialModelGraphviz;

    refinementList.push_back(
      std::move(refinementData)
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      Log::log(Log::Level::Warning)
        << "[" << currentStructureNumber << "]: "
        << "Third stage of refinement fails. Loosening factor was "
        << configuration.spatialModelLoosening
        <<  "\n";
      if(reachedMaxIterations) {
        Log::log(Log::Level::Warning) << "- Reached max iterations.\n";
      }

      if(notAllChiralitiesCorrect) {
        Log::log(Log::Level::Warning) << "- Not all chiral constraints have the correct sign.\n";
      }

      if(!structureAcceptable) {
        Log::log(Log::Level::Warning) << "- The final structure is unacceptable.\n";
        if(Log::isSet(Log::Particulars::DGStructureAcceptanceFailures)) {
          explainAcceptanceFailure(
            refinementFunctor,
            distanceBounds,
            transformedPositions
          );
        }
      }

      failures += 1;
    }
  }

  return refinementList;
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
      gatherDGInformation(moleculeCopy, configuration, engine)
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
  const Configuration& configuration
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
    *DGDataPtr = gatherDGInformation(molecule, configuration, randomnessEngine());
  }

  ReturnType results;
  results.reserve(numConformers);

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
      engine.seed(randomnessEngine()());

#else
      random::Engine& engine = randomnessEngine();
#endif

      /* We have to handle any and all exceptions here if this is a parallel
       * environment
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

#pragma omp critical(collectConformer)
        {
          results.push_back(std::move(conformerResult));
        }
      } catch(...) {
#pragma omp critical(collectConformer)
        {
          // Add an unknown error
          results.push_back(static_cast<DGError>(0));
        }
      } // end catch
    } // end pragma omp for private(DGDataPtr)
  } // end pragma omp parallel

  return results;
}

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine
