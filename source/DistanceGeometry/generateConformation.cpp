#include "DistanceGeometry/generateConformation.h"

#include "DistanceGeometry/MetricMatrix.h"
#include "DistanceGeometry/RefinementProblem.h"
#include "DistanceGeometry/dlibAdaptors.h"
#include "DistanceGeometry/dlibDebugAdaptors.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "GraphAlgorithms.h"
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
        vectorizedPositions.template segment<3>(
          static_cast<Eigen::Index>(dimensionality) * i
        )
      }
    );
  }

  return positions;
}

Delib::PositionCollection convertToPositionCollection(
  const dlib::matrix<double, 0, 1>& vectorizedPositions
) {
  const Eigen::Index dimensionality = 4;
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
  const Stereocenters::ChiralityConstraintPrototype& prototype
) {
  if(prototype.target == Stereocenters::ChiralityConstraintTarget::Flat) {
    // TODO introduce flatness tolerance
    return ChiralityConstraint {
      prototype.indices,
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
        prototype.indices.at(i),
        prototype.indices.at(j)
      );

      const double& upperBound = bounds.upperBound(
        prototype.indices.at(i),
        prototype.indices.at(j)
      );

      lowerMatrix(
        static_cast<Eigen::Index>(i) + 1,
        static_cast<Eigen::Index>(j) + 1
      ) = lowerBound * lowerBound;

      upperMatrix(
        static_cast<Eigen::Index>(i) + 1,
        static_cast<Eigen::Index>(j) + 1
      ) = upperBound * upperBound;
    }
  }

  const double boundFromLower = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    lowerMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

  const double boundFromUpper = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    upperMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

#ifndef NDEBUG
  // Log in Debug builds, provided particular is set
  Log::log(Log::Particulars::PrototypePropagatorDebugInfo) 
    << "Lower bound Cayley-menger matrix (determinant=" << boundFromLower <<  "): " << std::endl << lowerMatrix << std::endl
    << "Upper bound Cayley-menger matrix (determinant=" << boundFromUpper <<  "): " << std::endl << upperMatrix << std::endl;
#endif

  /* Check to ensure a volume is constructible from the passed bounds.
   * If the Cayley-Menger determinant is < 0, then it is impossible. The rather
   * lax criterion below follows from determinant calculation errors. Rather
   * be tolerant of approximately zero-volume chirality constraints (which MAY
   * happen) than have this fail.
   */
  assert(boundFromLower > -1e-7 && boundFromUpper > -1e-7);

  /* Although it is tempting to assume that the Cayley-Menger determinant using
   * the lower bounds is smaller than the one using upper bounds, this is NOT
   * true. We cannot a priori know which of both yields the lower or upper bounds
   * on the 3D volume, and hence must ensure only that the ordering is preserved
   * in the generation of the ChiralityConstraint, which checks that the lower
   * bound on the volume is certainly smaller than the upper one.
   *
   * You can check this assertion with a CAS. The relationship between both
   * determinants (where u_ij = l_ij + Δ) is wholly indeterminant, i.e. no
   * logical operator (<, >, <=, >=, ==) between both is true. It completely
   * depends on the individual values. Maybe in very specific cases one can
   * deduce some relationship, but not generally.
   */

  /* Abs is necessary to ensure that negative almost-zeros are interpreted as
   * positive near-zeros. Additionally, we do the adjustment by factor 8 here
   * (see above)
   */
  const double volumeFromLower = std::sqrt(std::fabs(boundFromLower) / 8);
  const double volumeFromUpper = std::sqrt(std::fabs(boundFromUpper) / 8);

  if(prototype.target == Stereocenters::ChiralityConstraintTarget::Negative) {
    /* For negative case, be aware inversion around 0 (*= -1) flips the inequality
     *      lower <  upper | *(-1)
     * ->  -lower > -upper
     *
     * So we construct it with the previous upper bound negated as the lower bound
     * and similarly for the previous lower bound.
     */
    return ChiralityConstraint {
      prototype.indices,
      - std::max(volumeFromLower, volumeFromUpper),
      - std::min(volumeFromLower, volumeFromUpper),
    };
  }

  // Regular case (Positive target)
  return ChiralityConstraint {
    prototype.indices,
    std::min(volumeFromLower, volumeFromUpper),
    std::max(volumeFromLower, volumeFromUpper),
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
  const MoleculeSpatialModel::DistanceMethod& distanceMethod
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
        const Stereocenters::ChiralityConstraintPrototype& prototype
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
          dlibPositions(
            static_cast<dlibIndexType>(i) * 4  + 1
          ) *= -1;
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
  const MoleculeSpatialModel::DistanceMethod& distanceMethod
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
        const Stereocenters::ChiralityConstraintPrototype& prototype
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
          dlibPositions(
            static_cast<dlibIndexType>(i) * 4 + 1
          ) *= -1;
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
  const MoleculeSpatialModel::DistanceMethod& distanceMethod
) {
  const auto& adjacencies = molecule.getAdjacencyList();

  // Initialize the return object
  MoleculeDGInformation data {adjacencies.numAtoms()};

  // Generate a spatial model from the molecular graph and stereocenters
  MoleculeSpatialModel spatialModel {
    adjacencies,
    molecule.stereocenters,
    distanceMethod
  };

  // Generate the distance bounds from the spatial model
  spatialModel.addDefaultDihedrals();
  data.distanceBounds = spatialModel.makeDistanceBounds();

  // Extract the chirality constraint prototypes
  data.chiralityConstraintPrototypes = spatialModel.getChiralityPrototypes();

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
