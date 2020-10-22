/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 *
 * This file reimplements Distance Geometry with more invasive debug structures
 * and data collection.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "boost/functional/hash.hpp"

#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/DistanceGeometry/ConformerGeneration.h"
#include "Molassembler/DistanceGeometry/EigenRefinement.h"
#include "Molassembler/DistanceGeometry/ExplicitBoundsGraph.h"
#include "Molassembler/DistanceGeometry/MetricMatrix.h"
#include "Molassembler/DistanceGeometry/RefinementMeta.h"
#include "Molassembler/DistanceGeometry/TetrangleSmoothing.h"
#include "Molassembler/IO.h"
#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/Log.h"
#include "Molassembler/GraphAlgorithms.h"

#include "Molassembler/Temple/Adaptors/Enumerate.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Optimization/Lbfgs.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include <fstream>
#include <iomanip>

namespace Scine {
namespace Molassembler {
namespace {

template<typename T, typename U>
std::pair<double, unsigned> covarianceImpl(const T& t, const U& u) {
  // Single-pass covariance calculation from Wikipedia
  double mu_x = 0;
  double mu_y = 0;
  double C = 0;
  unsigned N = 0;

  Temple::forEach(
    Temple::Adaptors::zip(t, u),
    [&](const double x, const double y) {
      N += 1;
      const double dx = x - mu_x;
      mu_x += dx / N;
      mu_y += (y - mu_y) / N;
      C += dx * (y - mu_y);
    }
  );

  return std::make_pair(C, N);
}

template<typename T, typename U>
double sampleCovariance(const T& t, const U& u) {
  auto parts = covarianceImpl(t, u);
  if(parts.second < 2) {
    throw std::logic_error("Need at least two values to calculate covariance");
  }
  return parts.first / (parts.second - 1);
}

} // namespace

namespace DistanceGeometry {

struct TermTraceVisitor {
  using IndexPair = std::pair<AtomIndex, AtomIndex>;

  template<typename Map, typename Key>
  static void addToMap(Map& map, const Key& key, double value) {
    const auto findIter = map.find(key);
    if(findIter == map.end()) {
      std::vector<double> values;
      values.push_back(value);
      map.emplace(key, values);
    } else {
      findIter->second.push_back(value);
    }
  }

  void distanceTerm(unsigned i, unsigned j, double value) {
    addToMap(distanceTerms, IndexPair {i, j}, value);
  }

  void chiralTerm(const ChiralConstraint::SiteSequence& sites, double value) {
    addToMap(chiralTerms, sites, value);
  }

  void dihedralTerm(const DihedralConstraint::SiteSequence& sites, double value) {
    addToMap(dihedralTerms, sites, value);
  }

  void fourthDimensionTerm(unsigned i, double value) {
    addToMap(fourthDimensionTerms, i, value);
  }

  using NameCovariancePair = std::pair<std::string, double>;
  using CovarianceContainer = std::vector<NameCovariancePair>;

  template<typename Key, typename DataContainer, typename Stringifier>
  static CovarianceContainer accumulateHighestCovariance(
    CovarianceContainer covariances,
    const DataContainer& container,
    const Key& deduplicationKey,
    const std::vector<double>& referenceData,
    unsigned count,
    Stringifier&& stringifier
  ) {
    return Temple::accumulate(
      container,
      std::move(covariances),
      [&](auto best, const auto& iterPair) {
        // Avoid calculating autocovariances
        if(iterPair.first == deduplicationKey) {
          return best;
        }

        const double covariance = sampleCovariance(iterPair.second, referenceData);

        if(best.size() < count) {
          best.emplace_back(stringifier(iterPair.first), covariance);
          return best;
        }

        const auto findIter = std::find_if(
          std::begin(best),
          std::end(best),
          [&](const auto& currentPair) {
            return currentPair.second < covariance;
          }
        );

        if(findIter != std::end(best)) {
          *findIter = NameCovariancePair {
            stringifier(iterPair.first),
            covariance
          };
        }

        return best;
      }
    );
  }

  void highestCovariance(std::ostream& os, unsigned i, unsigned j, unsigned count) const {
    const IndexPair distanceKey {i, j};
    const auto& referenceValues = distanceTerms.at(distanceKey);

    auto best = accumulateHighestCovariance(
      CovarianceContainer {},
      distanceTerms,
      IndexPair {i, j},
      referenceValues,
      count,
      [](auto x) { return "Distance " + Temple::stringify(x); }
    );

    best = accumulateHighestCovariance(
      std::move(best),
      chiralTerms,
      ChiralConstraint::SiteSequence {},
      referenceValues,
      count,
      [](auto x) { return "Chiral " + Temple::stringify(x); }
    );

    best = accumulateHighestCovariance(
      std::move(best),
      dihedralTerms,
      DihedralConstraint::SiteSequence {},
      referenceValues,
      count,
      [](auto x) { return "Dihedral " + Temple::stringify(x); }
    );

    best = accumulateHighestCovariance(
      std::move(best),
      fourthDimensionTerms,
      std::numeric_limits<unsigned>::max(),
      referenceValues,
      count,
      [](auto x) { return "4D " + Temple::stringify(x); }
    );

    Temple::sort(best, [](const auto& lhs, const auto& rhs) {
      return lhs.second > rhs.second;
    });

    for(const auto& nameCovariancePair : best) {
      os << nameCovariancePair.second << ": " << nameCovariancePair.first << "\n";
    }
  }

  std::unordered_map<IndexPair, std::vector<double>, boost::hash<IndexPair>> distanceTerms;
  std::unordered_map<ChiralConstraint::SiteSequence, std::vector<double>, boost::hash<ChiralConstraint::SiteSequence>> chiralTerms;
  std::unordered_map<DihedralConstraint::SiteSequence, std::vector<double>, boost::hash<DihedralConstraint::SiteSequence>> dihedralTerms;
  std::unordered_map<AtomIndex, std::vector<double>> fourthDimensionTerms;
};

/**
 * @brief Writes messages to the log explaining why a structure was rejected
 *
 * @see finalStructureAcceptable for criteria on acceptability
 *
 * A final structure is acceptable if
 * - All distance bounds are within 0.5 of either the lower or upper boundary
 * - All chiral constraints are within 0.5 of either the lower or upper boundary
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param refinement The refinement functor
 * @param bounds The distance bounds
 * @param positions The final positions from a refinement
 */
template<class RefinementType, typename PositionType>
void explainAcceptanceFailure(
  const RefinementType& refinement,
  const DistanceBoundsMatrix& bounds,
  const PositionType& positions,
  const TermTraceVisitor& termTrace
) {
  struct AcceptanceFailureExplainer {
    const double deviationThreshold = 0.5;
    bool earlyExit = false;
    bool value = true;
    std::ostream& log = Log::log(Log::Particulars::DgStructureAcceptanceFailures);
    const DistanceBoundsMatrix& distanceBounds;
    const TermTraceVisitor& trace;

    AcceptanceFailureExplainer(const DistanceBoundsMatrix& passBounds, const TermTraceVisitor& termTrace)
      : distanceBounds(passBounds), trace(termTrace) {}

    void distanceOverThreshold(AtomIndex i, AtomIndex j, double distance) {
      log << "Distance constraint " << i << " - " << j << " : ["
        << distanceBounds.lowerBound(i, j) << ", " << distanceBounds.upperBound(i, j)
        << "] deviation over threshold, distance is : " << distance << "\n";
      trace.highestCovariance(log, i, j, 5);
    }

    void chiralOverThreshold(const ChiralConstraint& chiral, double volume) {
      if(!chiral.targetVolumeIsZero()) {
        log << "Chiral constraint " << Temple::stringifyContainer(chiral.sites) << " : ["
          << chiral.lower << ", " << chiral.upper
          << "] deviation over threshold, volume is : " << volume << "\n";
      }
    }

    void dihedralOverThreshold(const DihedralConstraint& dihedral, double angle, double term) {
      log << "Dihedral constraint " << Temple::stringifyContainer(dihedral.sites)
        << " : [" << dihedral.lower << ", " << dihedral.upper
        << "], angle is : " << angle << ", term is: " << term << "\n";
    }
  };

  refinement.visitUnfulfilledConstraints(
    bounds,
    positions,
    AcceptanceFailureExplainer {bounds, termTrace}
  );
}

/**
 * @brief Writes messages to the log explaining the final contributions to
 *   the error function
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param refinement The refinement functor
 * @param bounds The distance bounds
 * @param positions The final positions from a refinement
 */
template<class RefinementType, typename PositionType>
void explainFinalContributions(
  const RefinementType& refinement,
  const DistanceBoundsMatrix& bounds,
  const PositionType& positions
) {
  struct FinalContributionsExplainer {
    const double deviationThreshold = 0;
    bool earlyExit = false;
    bool value = true;
    std::ostream& log = Log::log(Log::Particulars::DgFinalErrorContributions);
    const DistanceBoundsMatrix& distanceBounds;

    FinalContributionsExplainer(const DistanceBoundsMatrix& passBounds)
      : distanceBounds(passBounds) {}

    void distanceOverThreshold(AtomIndex i, AtomIndex j, double distance) {
      log << "Distance constraint " << i << " - " << j << " : ["
        << distanceBounds.lowerBound(i, j) << ", " << distanceBounds.upperBound(i, j)
        << "] is unsatisfied, distance is: " << distance << "\n";
    }

    void chiralOverThreshold(const ChiralConstraint& chiral, double volume) {
      log << "Chiral constraint " << Temple::stringify(chiral.sites) << " : ["
        << chiral.lower << ", " << chiral.upper
        << "] is unsatisfied, volume is: " << volume << "\n";
    }

    void dihedralOverThreshold(const DihedralConstraint& dihedral, double angle, double term) {
      log << "Dihedral constraint " << Temple::stringify(dihedral.sites)
        << " : [" << dihedral.lower << ", " << dihedral.upper
        << "] is unsatisfied, angle is: " << angle << ", term is " << term << "\n";
    }
  };

  refinement.visitUnfulfilledConstraints(
    bounds,
    positions,
    FinalContributionsExplainer {bounds}
  );
}



/**
 * @brief Debug class containing a step from refinement
 */
struct RefinementStepData {
  Eigen::VectorXd positions;
  double distanceError;
  double chiralError;
  double dihedralError;
  double fourthDimError;
  Eigen::VectorXd gradient;
  double proportionCorrectChiralConstraints;
  bool compress;

  RefinementStepData(
    const Eigen::VectorXd& passPositions,
    const double passDistanceError,
    const double passChiralError,
    const double passDihedralError,
    const double passFourthDimError,
    const Eigen::VectorXd& passGradient,
    const double passProportionCorrectChiralConstraints,
    const bool passCompress
  ) : positions(passPositions),
      distanceError(passDistanceError),
      chiralError(passChiralError),
      dihedralError(passDihedralError),
      fourthDimError(passFourthDimError),
      gradient(passGradient),
      proportionCorrectChiralConstraints(passProportionCorrectChiralConstraints),
      compress(passCompress)
  {}
};

/**
 * @brief Debug data class containing data on a refinement
 */
struct RefinementData {
  std::list<RefinementStepData> steps;
  std::vector<ChiralConstraint> constraints;
  double looseningFactor;
  bool isFailure;
  std::string spatialModelGraphviz;
};

namespace Detail {

template<typename EigenRefinementType>
struct InversionOrIterLimitStop {
  using VectorType = typename EigenRefinementType::VectorType;
  using FloatType = typename EigenRefinementType::FloatingPointType;

  InversionOrIterLimitStop(
    const unsigned passIter,
    const EigenRefinementType& functor
  ) : iterLimit(passIter),
      refinementFunctorReference(functor)
  {}

  template<typename StepValues>
  bool shouldContinue(unsigned iteration, const StepValues& step) {
    return (
      iteration < iterLimit
      && refinementFunctorReference.proportionChiralConstraintsCorrectSign < 1.0
      && (step.parameters.proposed - step.parameters.current).norm() > minParameterDiffNorm
    );
  }

  const unsigned iterLimit;
  const EigenRefinementType& refinementFunctorReference;
  double minParameterDiffNorm = 1e-3;
};

template<typename FloatType>
struct GradientOrIterLimitStop {
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

  template<typename StepValues>
  bool shouldContinue(unsigned iteration, const StepValues& step) {
    return (
      iteration < iterLimit
      && step.gradients.current.template cast<double>().norm() > gradNorm
      && (step.parameters.proposed - step.parameters.current).norm() > minParameterDiffNorm
    );
  }

  unsigned iterLimit = 10000;
  double gradNorm = 1e-5;
  double minParameterDiffNorm = 1e-3;
};

template<unsigned dimensionality>
Eigen::Vector3d averagePosition(
  const Eigen::VectorXd& linearPositions,
  const std::vector<AtomIndex>& indices
) {
  if(indices.size() == 1) {
    return linearPositions.template segment<3>(dimensionality * indices.front());
  }

  Eigen::Vector3d average = Eigen::Vector3d::Zero();
  for(AtomIndex i : indices) {
    average += linearPositions.template segment<3>(dimensionality * i);
  }
  average /= indices.size();
  return average;
}

template<unsigned dimensionality>
void twistRotatableDihedrals(
  Eigen::Ref<Eigen::VectorXd> positions,
  const std::vector<DihedralConstraint>& constraints,
  const MoleculeDGInformation::GroupMapType& rotatableGroupMap
) {
  for(const DihedralConstraint& constraint : constraints) {
    const AtomIndex j = constraint.sites.at(1).front();
    const AtomIndex k = constraint.sites.at(2).front();
    BondIndex bond {j, k};

    const auto findIter = rotatableGroupMap.find(bond);
    if(findIter == std::end(rotatableGroupMap)) {
      continue;
    }

    const MoleculeDGInformation::RotatableGroup& group = findIter->second;

    const Eigen::Vector3d jVector = positions.template segment<3>(dimensionality * j);
    const Eigen::Vector3d kVector = positions.template segment<3>(dimensionality * k);

    const double measuredDihedral = Cartesian::dihedral(
      averagePosition<dimensionality>(positions, constraint.sites.at(0)),
      jVector,
      kVector,
      averagePosition<dimensionality>(positions, constraint.sites.at(3))
    );

    const double targetDihedral = Cartesian::dihedralAverage(
      constraint.lower,
      constraint.upper
    );

    Eigen::Vector3d axisVector = (kVector - jVector).normalized();
    if(group.side == j) {
      axisVector *= -1;
    }

    const Eigen::Matrix3d rotation = Eigen::AngleAxisd {
      targetDihedral - measuredDihedral,
      axisVector
    }.toRotationMatrix();

    const Eigen::Vector3d sideVector = positions.template segment<3>(dimensionality * group.side);
    for(AtomIndex i : group.vertices) {
      positions.template segment<3>(dimensionality * i) = sideVector + rotation * (
        positions.template segment<3>(dimensionality * i) - sideVector
      );
    }
  }
}

} // namespace Detail

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const Configuration& configuration,
  std::string& spatialModelGraphvizString,
  const bool applyTetrangleSmoothing,
  const bool printBounds
) {
  // Generate a spatial model from the molecular graph and stereopermutators
  SpatialModel spatialModel {molecule, configuration};
  spatialModelGraphvizString = spatialModel.dumpGraphviz();

  // Extract gathered data
  MoleculeDGInformation data;
  data.bounds = spatialModel.makePairwiseBounds();
  data.chiralConstraints = spatialModel.getChiralConstraints();
  data.dihedralConstraints = spatialModel.getDihedralConstraints();
  data.rotatableGroups = MoleculeDGInformation::make(data.dihedralConstraints, molecule);

  if(applyTetrangleSmoothing) {
    /* Add implicit lower and upper bounds */
    const AtomIndex N = molecule.graph().N();
    for(AtomIndex i = 0; i < N; ++i) {
      for(AtomIndex j = i + 1; j < N; ++j) {
        double& lower = data.bounds(j, i);
        double& upper = data.bounds(i, j);

        if(lower == 0.0 && upper == 0.0) {
          double vdwLowerBound = (
            AtomInfo::vdwRadius(molecule.graph().elementType(i))
            + AtomInfo::vdwRadius(molecule.graph().elementType(j))
          );

          lower = vdwLowerBound;
          upper = 100.0;
        }
      }
    }

    // Triangle smooth
    DistanceBoundsMatrix::smooth(data.bounds);
    // Tetrangle smooth
    unsigned iterations = tetrangleSmooth(data.bounds);
    std::cout << "Applied " << iterations << " iterations of tetrangle smoothing\n";
  }

  if(printBounds) {
    const AtomIndex N = molecule.graph().N();
    for(AtomIndex i = 0; i < N; ++i) {
      auto iGraphDistances = distance(i, molecule.graph());

      for(AtomIndex j = i + 1; j < N; ++j) {
        const double lower = data.bounds(j, i);
        const double upper = data.bounds(i, j);

        if(lower != 0.0 || upper != 0.0) {
          std::cout << i << " - " << j
            << " = [" << lower << ", " << upper << "], width = " << (upper - lower)
            << ", graph distance " << iGraphDistances.at(j) << "\n";
        }
      }
    }
  }

  return data;
}

/*! @brief Logging, not throwing mostly identical implementation to run()
 *
 * A logging, not throwing, mostly identical implementation of
 * runDistanceGeometry that returns detailed intermediate data from
 * refinements, while run returns only the final conformers, which may
 * also be translated and rotated to satisfy fixed position constraints.
 */
std::list<RefinementData> debugRefinement(
  const Molecule& molecule,
  unsigned numConformers,
  const Configuration& configuration,
  bool applyTetrangleSmoothing,
  bool printBounds
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

  MoleculeDGInformation DgData;
  std::string spatialModelGraphviz;

  if(!regenerateEachStep) { // Collect once, keep all the time
    DgData = gatherDGInformation(
      molecule,
      configuration,
      spatialModelGraphviz,
      applyTetrangleSmoothing,
      printBounds
    );
  }

  unsigned failures = 0;
  for(
    unsigned currentStructureNumber = 0;
    // Failed optimizations do not count towards successful completion
    currentStructureNumber < numConformers;
    ++currentStructureNumber
  ) {
    if(regenerateEachStep) {
      auto moleculeCopy = Detail::narrow(molecule, randomnessEngine());

      if(moleculeCopy.stereopermutators().hasZeroAssignmentStereopermutators()) {
        Log::log(Log::Level::Warning)
          << "After setting stereopermutators at random, this molecule has "
          << "stereopermutators with zero valid permutations!"
          << std::endl;
      }

      // Fetch the DG data from the molecule with no unassigned stereopermutators
      DgData = gatherDGInformation(
        moleculeCopy,
        configuration,
        spatialModelGraphviz,
        applyTetrangleSmoothing,
        printBounds
      );
    }

    std::list<RefinementStepData> refinementSteps;

    ExplicitBoundsGraph explicitGraph {
      molecule.graph().inner(),
      DgData.bounds
    };

    auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
    if(!distanceBoundsResult) {
      Log::log(Log::Level::Warning) << "Failure in distance bounds matrix construction: "
        << distanceBoundsResult.error().message() << "\n";
      failures += 1;

      if(regenerateEachStep) {
        auto moleculeCopy = Detail::narrow(molecule, randomnessEngine());

        SpatialModel model {moleculeCopy, configuration};
        model.writeGraphviz("DG-failure-spatial-model-" + std::to_string(currentStructureNumber) + ".dot");
      } else {
        SpatialModel model {molecule, configuration};
        model.writeGraphviz("DG-failure-spatial-model-" + std::to_string(currentStructureNumber) + ".dot");
      }

      continue;
    }

    DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

    /* No need to smooth the distance bounds, ExplicitBoundsGraph creates it
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
      DgData.chiralConstraints,
      DgData.dihedralConstraints
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
      TermTraceVisitor& termTracer;

      Observer(const Molecule& passMolecule, FullRefinementType& passFunctor, std::list<RefinementStepData>& passSteps, TermTraceVisitor& passTracer)
        : mol(passMolecule), functor(passFunctor), steps(passSteps), termTracer(passTracer) {}

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

        functor.visitTerms(positions, termTracer);
      }
    };

    TermTraceVisitor termTracer;
    Observer observer(molecule, refinementFunctor, refinementSteps, termTracer);

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
      Detail::InversionOrIterLimitStop<FullRefinementType> inversionChecker {
        configuration.refinementStepLimit,
        refinementFunctor
      };

      Temple::Lbfgs<FloatType, 32> optimizer;

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
    Detail::GradientOrIterLimitStop<FloatType> gradientChecker;
    gradientChecker.gradNorm = 1e-3;
    gradientChecker.iterLimit = configuration.refinementStepLimit - firstStageIterations;

    try {
      Temple::Lbfgs<FloatType, 32> optimizer;

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
        refinementData.constraints = DgData.chiralConstraints;
        refinementData.looseningFactor = configuration.spatialModelLoosening;
        refinementData.isFailure = true;
        refinementData.spatialModelGraphviz = spatialModelGraphviz;

        refinementList.push_back(
          std::move(refinementData)
        );

        if(Log::particulars.count(Log::Particulars::DgFinalErrorContributions) > 0) {
          explainFinalContributions(
            refinementFunctor,
            distanceBounds,
            transformedPositions
          );
        }

        continue; // this triggers a new structure to be generated
    }

    /* Twist all freely rotatable dihedrals to their target values to avoid
     * conflicts between distance and dihedral errors to prevent rotations to
     * target values.
     */
    Detail::twistRotatableDihedrals<dimensionality>(
      transformedPositions,
      DgData.dihedralConstraints,
      DgData.rotatableGroups
    );

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

    if(Log::particulars.count(Log::Particulars::DgFinalErrorContributions) > 0) {
      explainFinalContributions(
        refinementFunctor,
        distanceBounds,
        transformedPositions
      );
    }

    RefinementData refinementData;
    refinementData.steps = std::move(refinementSteps);
    refinementData.constraints = DgData.chiralConstraints;
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
        if(Log::isSet(Log::Particulars::DgStructureAcceptanceFailures)) {
          explainAcceptanceFailure(
            refinementFunctor,
            distanceBounds,
            transformedPositions,
            termTracer
          );
        }
      }

      failures += 1;
    }
  }

  return refinementList;
}


} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine

using namespace std::string_literals;
using namespace Scine;
using namespace Molassembler;

void writeProgressFile(
  const Molecule& mol,
  const std::string& baseFilename,
  const unsigned index,
  const Eigen::VectorXd& positions
) {
  const std::string filename = baseFilename + "-" + std::to_string(index) + ".mol";
  AngstromPositions angstromWrapper = DistanceGeometry::Detail::convertToAngstromPositions(
    DistanceGeometry::Detail::gather(positions)
  );
  IO::write(filename, mol, angstromWrapper);
}

void writeProgressFiles(
  const Molecule& mol,
  const std::string& baseFilename,
  const DistanceGeometry::RefinementData& refinementData
) {
  /* Write the progress file */
  std::string progressFilename = baseFilename + "-progress.csv"s;
  std::ofstream progressFile (progressFilename);

  progressFile << std::scientific;

  for(const auto& refinementStep : refinementData.steps) {
    progressFile
      << refinementStep.distanceError << ","
      << refinementStep.chiralError << ","
      << refinementStep.dihedralError << ","
      << refinementStep.fourthDimError << ","
      << refinementStep.gradient.norm() << ","
      << static_cast<unsigned>(refinementStep.compress) << ","
      << refinementStep.proportionCorrectChiralConstraints << "\n";
  }

  progressFile.close();

  const unsigned maxProgressFiles = 100;

  if(refinementData.steps.size() > maxProgressFiles) {
    // Determine 100 roughly equispaced conformations to write to POV files
    double stepLength = static_cast<double>(refinementData.steps.size()) / maxProgressFiles;
    auto listIter = refinementData.steps.begin();
    unsigned currentIndex = 0;
    for(unsigned i = 0; i < maxProgressFiles; ++i) {
      unsigned targetIndex = std::floor(i * stepLength);
      assert(targetIndex >= currentIndex && targetIndex < refinementData.steps.size());
      std::advance(listIter, targetIndex - currentIndex);
      currentIndex = targetIndex;

      writeProgressFile(
        mol,
        baseFilename,
        i,
        listIter->positions
      );
    }
  } else {
    for(const auto enumPair : Temple::Adaptors::enumerate(refinementData.steps)) {
      writeProgressFile(
        mol,
        baseFilename,
        enumPair.index,
        enumPair.value.positions
      );
    }
  }

  // Write the graphviz representation of that structure number's spatial model
  std::string graphvizFilename = baseFilename + "-spatial-model.dot"s;
  std::ofstream graphvizfile (graphvizFilename);
  graphvizfile << refinementData.spatialModelGraphviz;
  graphvizfile.close();
}

const std::string partialityChoices =
  "  0 - Four-Atom Metrization\n"
  "  1 - 10% Metrization\n"
  "  2 - All (default)\n";

int main(int argc, char* argv[]) {
/* Set program options from command-line arguments */
  // Defaults
  unsigned nStructures = 1;

  bool showFinalContributions = false;
  bool showChiralConstraints = false;
  bool applyTetrangleSmoothing = false;
  bool printBounds = false;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help,h", "Produce help message")
    (
      "num_conformers,n",
      boost::program_options::value<unsigned>(),
      "Set number of structures to generate"
    )
    (
      "input,i",
      boost::program_options::value<std::string>(),
      "Filename to read molecule from or SMILES string"
    )
    (
      "partiality,p",
      boost::program_options::value<unsigned>(),
      "Set metrization partiality option (Default: full)"
    )
    (
      "steps,s",
      boost::program_options::value<unsigned>(),
      "Alter the maximum number of refinement steps (Default: 10'000)"
    )
    (
      "final_contributions,f",
      boost::program_options::bool_switch(&showFinalContributions),
      "Show the final contributions to the refinement error functions"
    )
    (
      "chirals,c",
      boost::program_options::bool_switch(&showChiralConstraints),
      "Show chiral constraint representations"
    )
    (
      "tetrangle,t",
      boost::program_options::bool_switch(&applyTetrangleSmoothing),
      "Apply tetrangle smoothing once, prior to distance matrix generation"
    )
    (
      "bounds,b",
      boost::program_options::bool_switch(&printBounds),
      "Print the distance bounds matrix"
    )
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::positional_options_description positional_description;
  positional_description.add("input", 1);
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv).
    options(options_description).
    positional(positional_description).
    style(
      boost::program_options::command_line_style::unix_style
      | boost::program_options::command_line_style::allow_long_disguise
    ).run(),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  // Manage the results
  if(options_variables_map.count("help") > 0 || options_variables_map.count("input") != 1) {
    std::cout << options_description << std::endl;
    return 0;
  }

  if(options_variables_map.count("num_conformers") > 0) {
    unsigned argN = options_variables_map["num_conformers"].as<unsigned>();
    if(argN == 0) {
      std::cout << "Specified to generate zero structures. Exiting."
        << std::endl;
      return 0;
    }

    nStructures = argN;
  }

  DistanceGeometry::Partiality metrizationOption = DistanceGeometry::Partiality::All;
  if(options_variables_map.count("partiality") > 0) {
    unsigned index =  options_variables_map["partiality"].as<unsigned>();

    if(index > 2) {
      std::cout << "Specified metrization option is out of bounds. Valid choices are:\n"
        << partialityChoices;
      return 0;
    }

    metrizationOption = static_cast<DistanceGeometry::Partiality>(index);
  }

  Log::particulars.insert(Log::Particulars::DgStructureAcceptanceFailures);

  if(showFinalContributions) {
    Log::particulars.insert(Log::Particulars::DgFinalErrorContributions);
  }

  unsigned nSteps = 10000;
  if(options_variables_map.count("steps") > 0) {
    nSteps = options_variables_map["steps"].as<unsigned>();
  }

/* Generating work */
  std::string baseName;
  Molecule mol;

  const std::string input = options_variables_map["input"].as<std::string>();
  if(boost::filesystem::exists(input)) {
    try {
      mol = IO::read(input);
    } catch(...) {
      std::cout << "Input exists as filename, but could not be read!\n";
      return 1;
    }

    boost::filesystem::path filepath {input};
    baseName = filepath.stem().string();
  } else {
    std::cout << "Input is not found as file. Interpreting as SMILES string\n";

    // Input is possibly a SMILES string
    try {
      mol = IO::Experimental::parseSmilesSingleMolecule(input);
    } catch(...) {
      std::cout << "Input could not be interpreted as a SMILES string.\n";
      return 1;
    }

    baseName = "smiles";
  }

  std::ofstream graphFile(baseName +  "-graph.dot");
  graphFile << mol.dumpGraphviz();
  graphFile.close();

  DistanceGeometry::Configuration DgConfiguration;
  DgConfiguration.partiality = metrizationOption;
  DgConfiguration.refinementStepLimit = nSteps;

  auto debugData = DistanceGeometry::debugRefinement(
    mol,
    nStructures,
    DgConfiguration,
    applyTetrangleSmoothing,
    printBounds
  );

  for(const auto& enumPair : Temple::Adaptors::enumerate(debugData)) {
    const auto& structNum = enumPair.index;
    const auto& refinementData = enumPair.value;

    std::string structBaseName = baseName + "-"s + std::to_string(structNum);

    writeProgressFiles(mol, structBaseName, refinementData);

    IO::write(
      structBaseName + "-last.mol"s,
      mol,
      DistanceGeometry::Detail::convertToAngstromPositions(
        DistanceGeometry::Detail::gather(refinementData.steps.back().positions)
      )
    );

    if(showChiralConstraints) {
      std::cout << "Chiral constraints (four atom index sets and bounds on signed volume in conformer) of refinement " << structNum << ":\n";
      for(const auto& constraint : refinementData.constraints) {
        std::cout << Temple::stringify(constraint.sites) << " -> [" << constraint.lower << ", " << constraint.upper << "]\n";
      }
    }
  }

  auto failures = Temple::sum(
    Temple::map(
      debugData,
      [](const auto& refinementData) -> unsigned {
        return static_cast<unsigned>(refinementData.isFailure);
      }
    )
  );

  if(failures > 0) {
    std::cout << "WARNING: " << failures << " refinements failed.\n";
  }
}
