#ifndef INCLUDE_MOLASSEMBLER_ANALYSIS_DEBUG_EIGEN_REFINEMENT_H
#define INCLUDE_MOLASSEMBLER_ANALYSIS_DEBUG_EIGEN_REFINEMENT_H

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/RefinementMeta.h"
#include "molassembler/DistanceGeometry/EigenRefinement.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "molassembler/IO.h"

#include "temple/LBFGS.h"
// #include "Utils/Optimizer/GradientBased/LBFGS.h"

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {

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

  bool checkConvergence(
    const VectorType& /* parameters */,
    FloatType /* value */,
    const VectorType& /* gradients */
  ) {
    return refinementFunctorReference.proportionChiralConstraintsCorrectSign >= 1.0;
  }

  bool checkMaxIterations(unsigned currentIteration) {
    return currentIteration >= iterLimit;
  }
};

//template<typename EigenRefinementType>
//struct InversionOrIterLimitStop final : public Utils::GradientBasedCheck {
//  const unsigned iterLimit;
//  const EigenRefinementType& refinementFunctorReference;
//
//
//  using VectorType = typename EigenRefinementType::VectorType;
//  using FloatType = typename EigenRefinementType::FloatingPointType;
//
//  InversionOrIterLimitStop(
//    const unsigned passIter,
//    const EigenRefinementType& functor
//  ) : iterLimit(passIter),
//      refinementFunctorReference(functor)
//  {}
//
//  virtual bool checkConvergence(
//    const VectorType& /* parameters */,
//    FloatType /* value */,
//    const VectorType& /* gradients */
//  ) final {
//    return refinementFunctorReference.proportionChiralConstraintsCorrectSign >= 1.0;
//  }
//
//  virtual bool checkMaxIterations(unsigned currentIteration) final {
//    return currentIteration >= iterLimit;
//  }
//
//  void addSettingsDescriptors(Utils::UniversalSettings::DescriptorCollection& /* collection */) final {}
//  void applySettings(const Utils::Settings& /* s */) final {}
//};

template<typename FloatType>
struct GradientOrIterLimitStop {
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

  bool checkConvergence(
    const VectorType& /* parameters */,
    FloatType /* value */,
    const VectorType& gradients
  ) {
    return gradients.template cast<double>().norm() <= gradNorm;
  }

  bool checkMaxIterations(unsigned currentIteration) {
    return currentIteration >= iterLimit;
  }

  unsigned iterLimit = 10000;
  double gradNorm = 1e-5;
};

struct LBFGSOptimizerParameters {
  // These are the default LBFGS parameters
  double c1 = 1e-4;
  double c2 = 0.9;
  double stepLength = 0.1;
};

template<unsigned dimensionality, typename FloatType, bool SIMD>
std::list<RefinementData> debugEigenRefinement(
  const Molecule& molecule,
  const unsigned numConformers,
  const Configuration& configuration,
  LBFGSOptimizerParameters optimizerParameters = {}
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
        auto moleculeCopy = detail::narrow(molecule);

        SpatialModel model {moleculeCopy, configuration};
        model.writeGraphviz("DG-failure-spatial-model-" + std::to_string(currentStructureNumber) + ".dot");
      } else {
        SpatialModel model {molecule, configuration};
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

/* Here is where Eigen specific stuff starts ----------------------- */
    using FullRefinementType = EigenRefinementProblem<dimensionality, FloatType, SIMD>;
    using VectorType = typename FullRefinementType::VectorType;

    // Vectorize positions
    VectorType transformedPositions = Eigen::Map<Eigen::VectorXd>(
      embeddedPositions.data(),
      embeddedPositions.cols() * embeddedPositions.rows()
    ).template cast<FloatType>().eval();

    const unsigned N = transformedPositions.size() / dimensionality;

    auto squaredBounds = static_cast<Eigen::MatrixXd>(
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
          dlib::mat(positions.template cast<double>().eval()),
          distanceError,
          chiralError,
          dihedralError,
          fourthDimensionError,
          dlib::mat(gradient.template cast<double>().eval()),
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
      InversionOrIterLimitStop<FullRefinementType> inversionChecker {
        configuration.refinementStepLimit,
        refinementFunctor
      };

      temple::LBFGS<FloatType, 64> optimizer;
      optimizer.c1 = optimizerParameters.c1;
      optimizer.c2 = optimizerParameters.c2;
      optimizer.stepLength = optimizerParameters.stepLength;

      try {
        firstStageIterations = optimizer.optimize(
          transformedPositions,
          refinementFunctor,
          inversionChecker,
          observer
        );
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
    GradientOrIterLimitStop<FloatType> gradientChecker;
    gradientChecker.gradNorm = 1e-3;
    gradientChecker.iterLimit = configuration.refinementStepLimit - firstStageIterations;

    try {
      temple::LBFGS<FloatType, 64> optimizer;
      optimizer.c1 = optimizerParameters.c1;
      optimizer.c2 = optimizerParameters.c2;
      optimizer.stepLength = optimizerParameters.stepLength;

      secondStageIterations = optimizer.optimize(
        transformedPositions,
        refinementFunctor,
        gradientChecker,
        observer
      );
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
    gradientChecker = GradientOrIterLimitStop<FloatType> {};
    gradientChecker.gradNorm = 1e-3;
    gradientChecker.iterLimit = (
      configuration.refinementStepLimit
      - firstStageIterations
      - secondStageIterations
    );

    refinementFunctor.dihedralTerms = true;

    try {
      temple::LBFGS<FloatType, 64> optimizer;
      optimizer.c1 = optimizerParameters.c1;
      optimizer.c2 = optimizerParameters.c2;
      optimizer.stepLength = optimizerParameters.stepLength;

      thirdStageIterations = optimizer.optimize(
        transformedPositions,
        refinementFunctor,
        gradientChecker,
        observer
      );
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

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine

#endif
