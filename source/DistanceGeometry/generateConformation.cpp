#include "DistanceGeometry/generateConformation.h"

#include "DistanceGeometry/MetricMatrix.h"
#include "DistanceGeometry/DGRefinementProblem.h"
#include "DistanceGeometry/BFSConstraintCollector.h"
#include "AdjacencyListAlgorithms.h"
#include "TreeAlgorithms.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

#include "Log.h"

namespace MoleculeManip {

namespace DistanceGeometry {

namespace detail {

Eigen::Vector3d getPos(
  const Eigen::MatrixXd& positions,
  const AtomIndexType& index
) {
  Eigen::Vector3d retv;
  retv = positions.col(index);
  return retv;
}

double evaluateChiralityConstraint(
  const ChiralityConstraint& chiralityConstraint,
  const Eigen::MatrixXd& positions
) {
  AtomIndexType i, j, k, l;
  std::tie(i, j, k, l, std::ignore) = chiralityConstraint;

  // V = (1 - 4) * [ (2 - 4) x (3 - 4) ]
  return (
    (
      getPos(positions, i)
      - getPos(positions, l)
    ).dot(
      (
       getPos(positions, j)
       - getPos(positions, l)
      ).cross(
        getPos(positions, k)
        - getPos(positions, l)
      )
    )
  );
}

bool moreThanHalfChiralityConstraintsIncorrect(
  const Eigen::MatrixXd& positions,
  const std::vector<ChiralityConstraint>& chiralityConstraints
) {
  unsigned totalNonZeroConstraints = 0, incorrectNonZeroConstraints = 0;
  for(const auto& chiralityConstraint : chiralityConstraints) {
    auto& target = std::get<4>(chiralityConstraint);

    if(target != 0.0) {
      totalNonZeroConstraints += 1;
    }

    auto eval = evaluateChiralityConstraint(
      chiralityConstraint,
      positions
    );

    if( // can this be simplified? -> sign bit XOR?
      ( eval < 0 && target > 0)
      || (eval > 0 && target < 0)
    ) {
      incorrectNonZeroConstraints += 1;
    }
  }

  return (
    // if there are no non-zero constraints, return immediately
    totalNonZeroConstraints == 0 
    || ( // otherwise, do a proper check
      static_cast<double>(incorrectNonZeroConstraints) / totalNonZeroConstraints 
      > 0.5
    )
  );

}

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

PrototypePropagator::PrototypePropagator(const Eigen::MatrixXd& distancesMatrix) : distancesMatrix(distancesMatrix) {}

ChiralityConstraint PrototypePropagator::operator () (
  const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
) {
  using namespace Stereocenters;

  if(prototype.second == Stereocenter::ChiralityConstraintTarget::Flat) {
    return ChiralityConstraint {
      prototype.first[0],
      prototype.first[1],
      prototype.first[2],
      prototype.first[3],
      0.0
    };
  } else {
    /*
     * Calculate the volume. The target volume of the chirality constraint
     * created by the tetrahedron is calculated using internal coordinates (the
     * Cayley-Menger determinant), always leading to V > 0, so depending on the
     * current assignment, the sign of the result is switched.  The formula
     * used later in chirality constraint calculation for explicit coordinates
     * is adjusted by V' = 6 V to avoid an unnecessary factor, so we do that
     * here too:
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

    Eigen::Matrix<double, 5, 5> cayleyMenger;
    cayleyMenger.setZero();

    for(unsigned i = 0; i < 4; i++) {
      for(unsigned j = i + 1; j < 4; j++) {

        cayleyMenger(i + 1, j + 1) = pow(
          distancesMatrix(
            std::min(
              prototype.first[i],
              prototype.first[j]
            ),
            std::max(
              prototype.first[i],
              prototype.first[j]
            )
          ),
          2
        );
      }
    }

    // top row of cayleyMenger matrix
    for(unsigned i = 1; i < 5; i++) {
      cayleyMenger(0, i) = 1;
    }

#ifndef NDEBUG
      // Log in Debug builds, provided particular is set
      Log::log(Log::Particulars::PrototypePropagatorDebugInfo) 
        << "distances Matrix: " << distancesMatrix << std::endl
        << "Cayley-menger matrix: " << cayleyMenger << std::endl;
#endif

    auto determinant = static_cast<
      Eigen::Matrix<double, 5, 5>
    >(
      cayleyMenger.selfadjointView<Eigen::Upper>()
    ).determinant();

    /* If this fails, then we have not done distance selection properly. There
     * must be triangle inequality violations in the distances matrix.
     */
    assert(determinant >= 0);

    auto chiralityTarget = sqrt(
      determinant / 8.0
    );

    if(
      prototype.second 
      == Stereocenters::Stereocenter::ChiralityConstraintTarget::Negative
    ) {
      chiralityTarget *= -1;
    }

    return ChiralityConstraint {
      prototype.first[0],
      prototype.first[1],
      prototype.first[2],
      prototype.first[3],
      chiralityTarget
    };
  }
}

std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  const auto DGData = gatherDGInformation(molecule);

  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGConjugatedGradientDescentSolver;

  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  // TODO this will need adjustment when some experience exists
  stopCriteria.iterations = 1000; 
  stopCriteria.fDelta = 1e-5;

  DGConjugatedGradientDescentSolver.setStopCriteria(stopCriteria);

  // Final initializations
  std::list<Delib::PositionCollection> positions;
  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a 
   * failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 0.1; // must be 0 < x <= 1
  unsigned optimizationFailures = 0;

  // Begin
  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber - optimizationFailures < numStructures;
    currentStructureNumber += 1
  ) {
    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix(
      metrization
    );
    
    /* Get the chirality constraints by converting the prototypes found by the
     * collector into full chiralityConstraints using the distances matrix
     */
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      /* Partial application with distances matrix so we have a unary function
       * to perform the mapping with
       */
      detail::PrototypePropagator {distancesMatrix}
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

    /* If a count of chirality constraints reveals that more than half are
     * incorrect, we can invert the structure (by multiplying e.g. all y 
     * coordinates with -1) and then have more than half of chirality 
     * constraints correct! In the count, chirality constraints with a target
     * value of zero are not considered (this would skew the count as those
     * chirality constraints should not have to pass an energetic maximum to 
     * converge properly as opposed to tetrahedra with volume).
     */
    if(detail::moreThanHalfChiralityConstraintsIncorrect(
      embeddedPositions,
      chiralityConstraints
    )) {
      embeddedPositions.row(2) *= -1;
    }

    // Vectorize the positions for use with cppoptlib
    Eigen::VectorXd vectorizedPositions {
      Eigen::Map<Eigen::VectorXd>(
        embeddedPositions.data(),
        embeddedPositions.cols() * embeddedPositions.rows()
      )
    };

    // Run the actual minimization
    DGConjugatedGradientDescentSolver.minimize(problem, vectorizedPositions);

    // What to do if the optimization fails
    if(DGConjugatedGradientDescentSolver.status() == cppoptlib::Status::Continue) {
      optimizationFailures += 1;

      if(
        static_cast<double>(optimizationFailures) / currentStructureNumber 
        >= failureRatio
      ) {
        throw std::runtime_error("Refinement failures exceeded threshold!");
      }
    } else {
      // Add the result to positions
      positions.emplace_back(
        detail::convertToPositionCollection(
          vectorizedPositions,
          embedding
        )
      );
    }
  }

  // for every column in embedded?
  //   positions.push_back(Delib::Position([Eigen::Vector3d] column))
  return positions;
}


} // eo namespace detail

MoleculeDGInformation::MoleculeDGInformation(const unsigned& N) : distanceBounds(N) {}

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule
) {
  const auto& adjacencies = molecule.getAdjacencyList();

  MoleculeDGInformation data {adjacencies.numAtoms()};

  BFSConstraintCollector collector(
    adjacencies,
    molecule.stereocenters,
    data.distanceBounds
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
  return detail::generateEnsemble(
    molecule,
    numStructures,
    metrization,
    embedding
  );
}

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding
) {
  auto list = detail::generateEnsemble(
    molecule,
    1,
    metrization,
    embedding
  );

  return *list.begin();
}

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip
