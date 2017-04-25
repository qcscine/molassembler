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
  const Eigen::VectorXd& vectorizedPositions,
  const EmbeddingOption& embedding
) {
  const unsigned dimensionality = static_cast<unsigned>(embedding);
  assert(vectorizedPositions.size() % dimensionality == 0);

  Delib::PositionCollection positions;
  const unsigned N = vectorizedPositions.size() / dimensionality;

  for(unsigned i = 0; i < N; i++) {
    positions.push_back(
      Delib::Position {
        getEigen<3>(
          vectorizedPositions,
          i,
          dimensionality
        )
      }
    );
  }

  return positions;
}

DGResult runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding,
  const bool& useYInversionTrick,
  const BFSConstraintCollector::DistanceMethod& distanceMethod
) {
  DGResult resultObject;

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
  const double failureRatio = 0.1; // must be 0 < x <= 1

  // Begin
  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - resultObject.statistics.failures < numStructures;
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
    auto embeddedPositions = metric.embed(embedding);

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
      resultObject.statistics.failures += 1;

      if( // fail-if
        (
          static_cast<double>(resultObject.statistics.failures)
          / currentStructureNumber 
        ) >= failureRatio
      ) {
        throw std::runtime_error("Refinement failures exceeded threshold!");
      }
    } else {
      // Make a count of correct chirality constraints
      auto count = problem.countCorrectChiralityConstraints(vectorizedPositions);

      if(count.incorrectNonZeroChiralityConstraints > 0) {
        resultObject.statistics.failures += 1;
        
        if( // fail-if
          (
            static_cast<double>(resultObject.statistics.failures)
            / currentStructureNumber 
          ) >= failureRatio
        ) {
          throw std::runtime_error("Refinement failures exceeded threshold!");
        }
      } else {
        // Add the result to positions
        resultObject.ensemble.emplace_back(
          detail::convertToPositionCollection(
            vectorizedPositions,
            embedding
          )
        );
      }
    }
  }

  // for every column in embedded?
  //   positions.push_back(Delib::Position([Eigen::Vector3d] column))
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

  BFSConstraintCollector collector(
    adjacencies,
    molecule.stereocenters,
    data.distanceBounds,
    distanceMethod
  );

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

    // Gather the distance constraints via visitor side effect
    TreeAlgorithms::BFSVisit(
      rootPtr,
      collector,
      3
    );
  }

  data.chiralityConstraintPrototypes = collector.getChiralityPrototypes();

  return data;
}

std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  return detail::runDistanceGeometry(
    molecule,
    numStructures,
    metrization,
    embedding
  ).ensemble;
}

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  auto list = detail::runDistanceGeometry(
    molecule,
    1,
    metrization,
    embedding
  ).ensemble;

  return *list.begin();
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
