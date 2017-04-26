#include "DistanceGeometry/generateConformation.h"

#include "DistanceGeometry/MetricMatrix.h"
#include "DistanceGeometry/DGRefinementProblem.h"
#include "DistanceGeometry/BFSConstraintCollector.h"
#include "AdjacencyListAlgorithms.h"
#include "TreeAlgorithms.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

namespace MoleculeManip {

namespace DistanceGeometry {

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
        vectorizedPositions.template segment<3>(dimensionality * i)
      }
    );
  }

  return positions;
}

// Non-debug version of DG
std::list<Delib::PositionCollection> runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const bool& useYInversionTrick,
  const BFSConstraintCollector::DistanceMethod& distanceMethod
) {
  std::list<Delib::PositionCollection> ensemble;

  const auto DGData = gatherDGInformation(
    molecule,
    distanceMethod
  );

  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGConjugatedGradientDescentSolver;

  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  stopCriteria.fDelta = 1e-5;

  DGConjugatedGradientDescentSolver.setStopCriteria(stopCriteria);

  unsigned failures = 0;

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a 
   * failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 0.1; // allow only 10% failures in release

  // Begin
  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - failures < numStructures;
    currentStructureNumber += 1
  ) {
    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix(
      metrization
    );
    
    /* Get the chirality constraints by converting the prototypes found by the
     * collector into full chiralityConstraints
     */
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      /* The propagator needs a function that gives the distances between
       * pairs of atoms, in this case we use the distances matrix we generated
       * from the bounds (must satisfy triangle inequalities). The Propagator
       * then has a unary function call operator that takes prototypes and spits
       * out full chiralityConstraints.
       */
      detail::makePropagator(
        [&distancesMatrix](
          const AtomIndexType& i,
          const AtomIndexType& j
        ) -> double {
          return distancesMatrix(
            std::min(i, j),
            std::max(i, j)
          );
        }
      )
    );

    /* Instantiantiate the refinement problem and its solver, set the stop 
     * criteria
     */
    DGRefinementProblem<double> problem(
      chiralityConstraints,
      DGData.distanceBounds
    );

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize the positions for use with cppoptlib
    Eigen::VectorXd vectorizedPositions {
      Eigen::Map<Eigen::VectorXd>(
        embeddedPositions.data(),
        embeddedPositions.cols() * embeddedPositions.rows()
      )
    };

    if(useYInversionTrick) {
      /* If a count of chirality constraints reveals that more than half are
       * incorrect, we can invert the structure (by multiplying e.g. all y
       * coordinates with -1) and then have more than half of chirality
       * constraints correct! In the count, chirality constraints with a target
       * value of zero are not considered (this would skew the count as those
       * chirality constraints should not have to pass an energetic maximum to
       * converge properly as opposed to tetrahedra with volume).
       */
      if(!problem.moreThanHalfChiralityConstraintsCorrect(vectorizedPositions)) {
        problem.invertY(vectorizedPositions);
      }
    }

    // Run the actual minimization
    DGConjugatedGradientDescentSolver.minimize(problem, vectorizedPositions);

    // What to do if the optimization fails
    if(DGConjugatedGradientDescentSolver.status() == cppoptlib::Status::Continue) {
      failures += 1;

      if( // fail-if
        (
          static_cast<double>(failures)
          / numStructures 
        ) >= failureRatio
      ) {
        throw std::runtime_error("Refinement failures exceeded threshold!");
      }
    } else {
      // Make a count of correct chirality constraints
      auto count = problem.countCorrectChiralityConstraints(vectorizedPositions);

      if(count.incorrectNonZeroChiralityConstraints > 0) {
        failures += 1;
        
        if( // fail-if
          (
            static_cast<double>(failures)
            / numStructures 
          ) >= failureRatio
        ) {
          throw std::runtime_error("Refinement failures exceeded threshold!");
        }
      } else {
        // Add the result to positions
        ensemble.emplace_back(
          detail::convertToPositionCollection(vectorizedPositions)
        );
      }
    }
  }

  return ensemble;
}

// Debug version
DGDebugData debugDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const bool& useYInversionTrick,
  const BFSConstraintCollector::DistanceMethod& distanceMethod
) {
  DGDebugData resultObject;

  const auto DGData = gatherDGInformation(
    molecule,
    distanceMethod
  );

  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGConjugatedGradientDescentSolver;

  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  stopCriteria.fDelta = 1e-5;

  DGConjugatedGradientDescentSolver.setStopCriteria(stopCriteria);

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a 
   * failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 3; // must be > 0

  // Begin
  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - resultObject.failures < numStructures;
    currentStructureNumber += 1
  ) {
    std::list<RefinementStepData> refinementSteps;

    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix(
      metrization
    );
    
    /* Get the chirality constraints by converting the prototypes found by the
     * collector into full chiralityConstraints
     */
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      /* The propagator needs a function that gives the distances between
       * pairs of atoms, in this case we use the distances matrix we generated
       * from the bounds (must satisfy triangle inequalities). The Propagator
       * then has a unary function call operator that takes prototypes and spits
       * out full chiralityConstraints.
       */
      detail::makePropagator(
        [&distancesMatrix](
          const AtomIndexType& i,
          const AtomIndexType& j
        ) -> double {
          return distancesMatrix(
            std::min(i, j),
            std::max(i, j)
          );
        }
      )
    );

    /* Instantiantiate the refinement problem and its solver, set the stop 
     * criteria
     */
    DGRefinementProblem<double> problem(
      chiralityConstraints,
      DGData.distanceBounds
    );

    auto getGradient = [&](const Eigen::VectorXd& positions){
      Eigen::VectorXd gradientVector( positions.size() );
      problem.gradient(positions, gradientVector);
      return gradientVector;
    };

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize the positions for use with cppoptlib
    Eigen::VectorXd vectorizedPositions {
      Eigen::Map<Eigen::VectorXd>(
        embeddedPositions.data(),
        embeddedPositions.cols() * embeddedPositions.rows()
      )
    };

    if(useYInversionTrick) {
      // Add the structure before inversion
      refinementSteps.emplace_back(
        vectorizedPositions,
        problem.value(vectorizedPositions),
        getGradient(vectorizedPositions),
        false
      );

      /* If a count of chirality constraints reveals that more than half are
       * incorrect, we can invert the structure (by multiplying e.g. all y
       * coordinates with -1) and then have more than half of chirality
       * constraints correct! In the count, chirality constraints with a target
       * value of zero are not considered (this would skew the count as those
       * chirality constraints should not have to pass an energetic maximum to
       * converge properly as opposed to tetrahedra with volume).
       */
      if(!problem.moreThanHalfChiralityConstraintsCorrect(vectorizedPositions)) {
        problem.invertY(vectorizedPositions);
      }
    }

    // add a refinement step
    refinementSteps.emplace_back(
      vectorizedPositions,
      problem.value(vectorizedPositions),
      getGradient(vectorizedPositions),
      false
    );

    // initial step
    auto stepResult = DGConjugatedGradientDescentSolver.step(problem, vectorizedPositions);
    unsigned iterations = 1;

    // add a refinement step
    refinementSteps.emplace_back(
      vectorizedPositions,
      problem.value(vectorizedPositions),
      getGradient(vectorizedPositions),
      problem.compress
    );

    while(stepResult.status == cppoptlib::Status::Continue && iterations < 1e5) {
      stepResult = DGConjugatedGradientDescentSolver.step(problem, vectorizedPositions);

      iterations += 1;

      refinementSteps.emplace_back(
        vectorizedPositions,
        problem.value(vectorizedPositions),
        getGradient(vectorizedPositions),
        problem.compress
      );
    }

    // What to do if the optimization fails
    // Make a count of correct chirality constraints
    auto count = problem.countCorrectChiralityConstraints(vectorizedPositions);
    if(count.incorrectNonZeroChiralityConstraints > 0
      || iterations == 1e5
    ) {
      resultObject.failures += 1;
      
      if( // fail-if
        (
          static_cast<double>(resultObject.failures)
          / numStructures 
        ) >= failureRatio
      ) {
        throw std::runtime_error("Refinement failures exceeded threshold!");
      }
    } else {
      // Add the result to the debug data
      RefinementData refinementData;
      refinementData.steps = std::move(refinementSteps);
      refinementData.constraints = problem.constraints;

      resultObject.refinements.emplace_back(
        std::move(refinementData)
      );
    }
  }

  return resultObject;
}

} // namespace detail

MoleculeDGInformation::MoleculeDGInformation(const unsigned& N) : distanceBounds(N) {}

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const BFSConstraintCollector::DistanceMethod& distanceMethod
) {
  const auto& adjacencies = molecule.getAdjacencyList();

  MoleculeDGInformation data {adjacencies.numAtoms()};

  BFSConstraintCollector collector {
    adjacencies,
    molecule.stereocenters,
    data.distanceBounds,
    distanceMethod
  };

  /* For every atom in the molecule, generate a tree of height 3 to traverse
   * with our collector. This leads to some repeated work that could be 
   * optimized out at some point (TODO).
   */
  for(AtomIndexType i = 0; i < molecule.getNumAtoms(); i++) {
    // Make the tree
    auto rootPtr = AdjacencyListAlgorithms::makeTree(
      adjacencies,
      i,
      3 // max Depth of 3 to limit to up to dihedral length chains
    );

#ifndef NDEBUG
    Log::log(Log::Particulars::gatherDGInformationTrees) << "On atom " << i
      << rootPtr << std::endl;
#endif

    // Gather the distance constraints via visitor side effect
    TreeAlgorithms::BFSVisit(
      rootPtr,
      collector,
      3
    );
  }

  // Fix 0-length lower bounds and smooth the bounds once
  collector.finalizeBoundsMatrix();

  // Extract the chirality constraint prototypes
  data.chiralityConstraintPrototypes = collector.getChiralityPrototypes();

  return data;
}

std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization
) {
  return detail::runDistanceGeometry(
    molecule,
    numStructures,
    metrization
  );
}

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization
) {
  auto list = detail::runDistanceGeometry(
    molecule,
    1,
    metrization
  );
  
  assert(list.size() == 1);

  return *list.begin();
}

} // namespace DistanceGeometry

} // namespace MoleculeManip
