#include "DistanceGeometry/generateConformation.h"

#include "DistanceGeometry/MetricMatrix.h"
#include "DistanceGeometry/RefinementProblem.h"
#include "DistanceGeometry/dlibAdaptors.h"
#include "DistanceGeometry/dlibDebugAdaptors.h"
#include "DistanceGeometry/BFSConstraintCollector.h"
#include "AdjacencyListAlgorithms.h"
#include "TreeAlgorithms.h"

#include <dlib/optimization.h>
#include <Eigen/Dense>

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

Delib::PositionCollection convertToPositionCollection(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
) {
  const unsigned dimensionality = 4;
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

ChiralityConstraint propagate(
  const DistanceBoundsMatrix& bounds,
  const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
) {
  using namespace Stereocenters;

  if(prototype.second == Stereocenter::ChiralityConstraintTarget::Flat) {
    return ChiralityConstraint {
      prototype.first,
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
        prototype.first.at(i),
        prototype.first.at(j)
      );

      const double& upperBound = bounds.upperBound(
        prototype.first.at(i),
        prototype.first.at(j)
      );

      upperMatrix(i + 1, j + 1) = upperBound * upperBound;
      lowerMatrix(i + 1, j + 1) = lowerBound * lowerBound;
    }
  }

  const double lowerBound = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    lowerMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

  const double upperBound = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    upperMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

#ifndef NDEBUG
  // Log in Debug builds, provided particular is set
  Log::log(Log::Particulars::PrototypePropagatorDebugInfo) 
    << "Lower bound Cayley-menger matrix (determinant=" << lowerBound <<  "): " << std::endl << lowerMatrix << std::endl
    << "Upper bound Cayley-menger matrix (determinant=" << upperBound <<  "): " << std::endl << upperMatrix << std::endl;
#endif

  assert(lowerBound > -1e-7 && upperBound > -1e-7);
  assert(lowerBound < upperBound);

  /* For negative case, be aware inversion around 0 (*= -1) flips the inequality
   *      lower <  upper | *(-1)
   * ->  -lower > -upper
   *
   * So we construct it with the previous upper bound negated as the lower bound
   * and similarly for the previous lower bound.
   */
  if(prototype.second == Stereocenter::ChiralityConstraintTarget::Negative) {
    /* Abs is necessary to ensure that negative almost-zeros are interpreted as
     * positive near-zeros. Additionally, we do the adjustment by factor 8 here
     * (see above)
     */
    return ChiralityConstraint {
      prototype.first,
      - sqrt(std::fabs(upperBound) / 8.0),
      - sqrt(std::fabs(lowerBound) / 8.0)
    };
  }

  // Regular case (Positive target)
  return ChiralityConstraint {
    prototype.first,
    sqrt(std::fabs(lowerBound) / 8.0),
    sqrt(std::fabs(upperBound) / 8.0)
  };
}

void checkFailureRatio(
  const unsigned& failures,
  const unsigned& numStructures,
  const unsigned& failureRatio
) {
  if(static_cast<double>(failures) / numStructures >= failureRatio) {
    throw std::runtime_error("Refinement failures exceeded threshold!");
  }
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

  // Always:
//  stopCriteria.iterations = 1e5;
//  stopCriteria.gradNorm = 1e-5;

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
      [&DGData](
        const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
      ) -> ChiralityConstraint {
        return propagate(
          DGData.distanceBounds,
          prototype
        );
      }
    );

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

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
        for(unsigned i = 0; i < DGData.distanceBounds.N; i++) {
          dlibPositions(4 * i + 1) *= -1;
        }
      }
    }

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
          DGData.distanceBounds,
          chiralityConstraints
        ),
        errfGradient<false>(
          DGData.distanceBounds,
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
        checkFailureRatio(failures, numStructures, failureRatio);
        continue; // this triggers a new structure to be generated
      }
    }

    dlibAdaptors::iterationOrGradientNormStopStrategy refinementStopStrategy {
      10000, // iteration limit
      1e-5 // gradient limit
    };

    errfValue<true> refinementValueFunctor {
      DGData.distanceBounds,
      chiralityConstraints
    };

    // Refinement with penalty on fourth dimension is always necessary
    dlib::find_min(
      dlib::bfgs_search_strategy(),
      refinementStopStrategy,
      refinementValueFunctor,
      errfGradient<true>(
        DGData.distanceBounds,
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
      DGData.distanceBounds,
      chiralityConstraints,
      dlibPositions
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      failures += 1;

      if(static_cast<double>(failures) / numStructures >= failureRatio) {
        throw std::runtime_error("Refinement failures exceeded threshold!");
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

  /* If the ratio of failures/total optimizations exceeds this value,
   * the function throws. Remember that if an optimization is considered a 
   * failure is dependent only on the stopping criteria!
   */
  const double failureRatio = 3; // more lenient than production, must be > 0
  unsigned failures = 0;

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
      [&DGData](
        const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
      ) -> ChiralityConstraint {
        return propagate(
          DGData.distanceBounds,
          prototype
        );
      }
    );

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

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
        for(unsigned i = 0; i < DGData.distanceBounds.N; i++) {
          dlibPositions(4 * i + 1) *= -1;
        }
      }
    }

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
        DGData.distanceBounds,
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
          DGData.distanceBounds,
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
        // TODO Debug Log!
        checkFailureRatio(failures, numStructures, failureRatio);
        continue; // this triggers a new structure to be generated
      }
    }

    errfValue<true> refinementValueFunctor {
      DGData.distanceBounds,
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
        DGData.distanceBounds,
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
      DGData.distanceBounds,
      chiralityConstraints,
      dlibPositions
    );

    if(reachedMaxIterations || notAllChiralitiesCorrect || !structureAcceptable) {
      failures += 1;
      checkFailureRatio(failures, numStructures, failureRatio);
    } else {
      // Add the result to the debug data
      RefinementData refinementData;
      refinementData.steps = std::move(refinementSteps);
      refinementData.constraints = chiralityConstraints;

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
