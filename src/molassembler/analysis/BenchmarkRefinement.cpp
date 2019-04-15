/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/EigenRefinement.h"
#include "molassembler/DistanceGeometry/DlibRefinement.h"
#include "molassembler/DistanceGeometry/DlibAdaptors.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "temple/LBFGS.h"

#include "molassembler/IO.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Dense>

#include <chrono>
#include <iomanip>

using namespace Scine;
using namespace molassembler;

constexpr std::size_t nExperiments = 100;
constexpr std::size_t refinementStepLimit = 10000;
constexpr double refinementGradientTarget = 1e-5;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

struct TimingFunctor {
  virtual boost::optional<unsigned> value(
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    const Eigen::MatrixXd& positions
  ) = 0;

  virtual std::string name() = 0;
};

struct FunctorResults {
  unsigned count = 0;
  double timingAverage = 0;
  double timingStddev = 0;
  double iterationsAverage = 0;
  double iterationsStddev = 0;
};

template<size_t N>
std::vector<FunctorResults> timeFunctors(
  const Molecule& molecule,
  const std::vector<
    std::unique_ptr<TimingFunctor>
  >& functors
) {
  using namespace std::chrono;

  time_point<steady_clock> start, end;

  struct Counter {
    std::vector<double> timings;
    std::vector<double> iterations;
  };

  std::vector<Counter> counters(functors.size());

  for(unsigned n = 0; n < N; ++n) {
    // Prepare new data for the functors in every experiment
    DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

    const auto boundsList = spatialModel.makeBoundsList();

    const auto chiralConstraints = spatialModel.getChiralConstraints();
    const auto dihedralConstraints = spatialModel.getDihedralConstraints();

    DistanceGeometry::ExplicitGraph explicitGraph {
      molecule,
      boundsList
    };

    auto distancesMatrixResult = explicitGraph.makeDistanceMatrix(randomnessEngine());
    if(!distancesMatrixResult) {
      throw std::runtime_error(distancesMatrixResult.error().message());
    }

    auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
    if(!distanceBoundsResult) {
      throw std::runtime_error("Failure in distance bounds matrix construction: " + distanceBoundsResult.error().message());
    }

    DistanceGeometry::DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

    Eigen::MatrixXd squaredBounds = distanceBounds.access().cwiseProduct(
      distanceBounds.access()
    );

    auto metricMatrix = DistanceGeometry::MetricMatrix(
      std::move(distancesMatrixResult.value())
    );

    auto embeddedPositions = metricMatrix.embed();

    // Run all functors on the same
    for(unsigned functorIndex = 0; functorIndex < functors.size(); ++functorIndex) {
      start = steady_clock::now();
      auto resultOption = functors.at(functorIndex)->value(
        squaredBounds,
        chiralConstraints,
        dihedralConstraints,
        embeddedPositions
      );
      end = steady_clock::now();

      if(resultOption) {
        counters.at(functorIndex).timings.push_back(duration_cast<microseconds>(end - start).count());
        counters.at(functorIndex).iterations.push_back(resultOption.value());
      }
    }
  }

  return temple::map(
    counters,
    [](const Counter& counter) -> FunctorResults {
      FunctorResults result;

      result.count = counter.timings.size();
      if(!counter.timings.empty()) {
        result.timingAverage = temple::average(counter.timings);
        result.timingStddev = temple::stddev(counter.timings, result.timingAverage);
        result.iterationsAverage = temple::average(counter.iterations);
        result.iterationsStddev = temple::stddev(counter.iterations, result.iterationsAverage);
      }

      return result;
    }
  );
}

struct DlibFunctor final : public TimingFunctor {
  boost::optional<unsigned> value (
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    const Eigen::MatrixXd& positions
  ) final {
    using namespace DistanceGeometry;

    Eigen::MatrixXd positionCopy = positions;

    const unsigned N = positions.size() / 4;
    unsigned iterationCount = 0;

    /* Transform positions to dlib space */
    ErrorFunctionValue::Vector dlibPositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
        Eigen::Map<Eigen::VectorXd>(
          positionCopy.data(),
          positionCopy.cols() * positionCopy.rows()
        )
      )
    );

    // Transform squared bounds to dlib type
    const dlib::matrix<double, 0, 0> dlibSquaredBounds = dlib::mat(squaredBounds);

    /* Instantiate functors */
    ErrorFunctionValue valueFunctor {
      dlibSquaredBounds,
      chiralConstraints,
      dihedralConstraints
    };

    ErrorFunctionGradient gradientFunctor {
      dlibSquaredBounds,
      chiralConstraints,
      dihedralConstraints
    };

    /* Perform inversion trick */
    double initiallyCorrectChiralConstraints = valueFunctor
      .calculateProportionChiralConstraintsCorrectSign(dlibPositions);

    if(initiallyCorrectChiralConstraints < 0.5) {
      // Invert y coordinates
      for(unsigned i = 0; i < N; i++) {
        dlibPositions(
          static_cast<dlibIndexType>(i) * 4 + 1
        ) *= -1;
      }

      initiallyCorrectChiralConstraints = 1 - initiallyCorrectChiralConstraints;
    }

    /* Stage one: Invert chirals */
    if(initiallyCorrectChiralConstraints < 1) {
      unsigned firstStageIterations;

      dlibAdaptors::IterationOrAllChiralitiesCorrectStrategy inversionStopStrategy {
        firstStageIterations,
        valueFunctor,
        refinementStepLimit
      };

      try {
        dlib::find_min(
          dlib::bfgs_search_strategy(),
          inversionStopStrategy,
          valueFunctor,
          gradientFunctor,
          dlibPositions,
          0
        );

      } catch(...) {
        return boost::none;
      }

      /* Handle errors */
      if(firstStageIterations >= refinementStepLimit) {
        return boost::none;
      }

      if(valueFunctor.proportionChiralConstraintsCorrectSign < 1.0) {
        return boost::none;
      }

      iterationCount += firstStageIterations;
    }

    unsigned secondStageIterations;
    dlibAdaptors::IterationOrGradientNormStopStrategy refinementStopStrategy {
      secondStageIterations,
      refinementStepLimit,
      refinementGradientTarget
    };

    /* Stage two: refine, compressing the fourth dimension */
    valueFunctor.compressFourthDimension = true;
    gradientFunctor.compressFourthDimension = true;
    try {
      dlib::find_min(
        dlib::bfgs_search_strategy(),
        refinementStopStrategy,
        valueFunctor,
        gradientFunctor,
        dlibPositions,
        0
      );
    } catch(std::out_of_range& e) {
      return boost::none;
    }

    /* Error conditions */
    // Max iterations reached
    if(secondStageIterations >= refinementStepLimit) {
      return boost::none;
    }

    // Not all chiral constraints have the right sign
    if(valueFunctor.proportionChiralConstraintsCorrectSign < 1) {
      return boost::none;
    }

    iterationCount += secondStageIterations;

    return iterationCount;
  }

  std::string name() final {
    return "Dlib";
  }
};

// template<typename EigenRefinementType>
// struct InversionOrIterLimitStop final : public Utils::GradientBasedCheck {
//   const EigenRefinementType& refinementFunctorReference;
//   using VectorType = typename EigenRefinementType::VectorType;
//   using FloatType = typename EigenRefinementType::FloatingPointType;
//
//   InversionOrIterLimitStop(const EigenRefinementType& functor) : refinementFunctorReference(functor) {}
//
//   virtual bool checkConvergence(
//     const VectorType& /* parameters */,
//     FloatType /* value */,
//     const VectorType& /* gradients */
//   ) final {
//     return refinementFunctorReference.proportionChiralConstraintsCorrectSign >= 1.0;
//   }
//
//   virtual bool checkMaxIterations(unsigned currentIteration) final {
//     return currentIteration >= refinementStepLimit;
//   }
//
//   void addSettingsDescriptors(Utils::UniversalSettings::DescriptorCollection& /* collection */) final {}
//   void applySettings(const Utils::Settings& /* s */) final {}
// };

template<typename EigenRefinementType>
struct InversionOrIterLimitStop {
  const EigenRefinementType& refinementFunctorReference;

  InversionOrIterLimitStop(const EigenRefinementType& functor) : refinementFunctorReference(functor) {}

  template<typename VectorType, typename FloatType>
  bool checkConvergence(
    const Eigen::Ref<VectorType>& /* parameters */,
    FloatType /* value */,
    const VectorType& /* gradients */
  ) {
    return refinementFunctorReference.proportionChiralConstraintsCorrectSign >= 1.0;
  }

  bool checkMaxIterations(unsigned currentIteration) {
    return currentIteration >= refinementStepLimit;
  }
};

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
    return currentIteration >= refinementStepLimit;
  }

  double gradNorm = 1e-5;
};

struct OptimizerParameters {
  // These are the default LBFGS parameters
  unsigned maxm = 50;
  double c1 = 1e-4;
  double c2 = 0.9;
  double stepLength = 1.0;
};

template<unsigned dimensionality, typename FloatType, bool SIMD>
boost::optional<unsigned> eigenRefine(
  const Eigen::MatrixXd& squaredBounds,
  const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
  const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
  const Eigen::MatrixXd& positions,
  OptimizerParameters optimizerParameters = {}
) {
  unsigned iterationCount = 0;

  using FullRefinementType = DistanceGeometry::EigenRefinementProblem<dimensionality, FloatType, SIMD>;

  const unsigned N = positions.size() / dimensionality;
  /* Transfer positions into vector form */

  Eigen::MatrixXd copiedPositions = positions;
  typename FullRefinementType::VectorType transformedPositions = Eigen::Map<Eigen::VectorXd>(
    copiedPositions.data(),
    copiedPositions.cols() * copiedPositions.rows()
  ).template cast<FloatType>().eval();

  FullRefinementType refinementFunctor {
    squaredBounds,
    chiralConstraints,
    dihedralConstraints
  };

  double initiallyCorrectChiralConstraints = refinementFunctor.calculateProportionChiralConstraintsCorrectSign(transformedPositions);

  if(initiallyCorrectChiralConstraints < 0.5) {
    for(unsigned i = 0; i < N; ++i) {
      transformedPositions(dimensionality * i + 1) *= -1;
    }
    initiallyCorrectChiralConstraints = 1 - initiallyCorrectChiralConstraints;
  }

  temple::LBFGS<FloatType, 32> optimizer;
  optimizer.c1 = optimizerParameters.c1;
  optimizer.c2 = optimizerParameters.c2;
  optimizer.stepLength = optimizerParameters.stepLength;

  /* First stage: Invert all chiral constraints */
  if(initiallyCorrectChiralConstraints < 1) {
    InversionOrIterLimitStop<FullRefinementType> inversionChecker {refinementFunctor};

    unsigned iterations;
    try {
      iterations = optimizer.optimize(
        transformedPositions,
        refinementFunctor,
        inversionChecker
      );
    } catch (...) {
      return boost::none;
    }

    if(iterations >= refinementStepLimit) {
      return boost::none;
    }

    if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1.0) {
      return boost::none;
    }

    iterationCount += iterations;
  }

  /* Second stage: Refine */
  refinementFunctor.compressFourthDimension = true;

  GradientOrIterLimitStop<FloatType> gradientChecker;
  /*Utils::GradientBasedCheck gradientChecker;
  gradientChecker.maxIter = refinementStepLimit;
  gradientChecker.gradNorm = 1e-5;*/

  unsigned iterations;
  try {
    iterations = optimizer.optimize(
      transformedPositions,
      refinementFunctor,
      gradientChecker
    );
  } catch (...) {
    return boost::none;
  }

  if(iterations >= refinementStepLimit) {
    return boost::none;
  }

  if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1) {
    return boost::none;
  }

  iterationCount += iterations;

  return iterationCount;

  /*Utils::PositionCollection finalPositions(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    finalPositions.row(i) = transformedPositions.segment<3>(dimensionality * i);
  }

  return AngstromWrapper {std::move(finalPositions), LengthUnit::Angstrom};*/
}


template<unsigned dimensionality, typename FloatType, bool SIMD>
struct EigenFunctor final : public TimingFunctor {
  boost::optional<unsigned> value(
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    const Eigen::MatrixXd& positions
  ) final {
    return eigenRefine<dimensionality, FloatType, SIMD>(
      squaredBounds,
      chiralConstraints,
      dihedralConstraints,
      positions
    );
  }

  std::string name() final {
    std::string composite = "Eigen<" + std::to_string(dimensionality) +", ";

    if(std::is_same<FloatType, double>::value) {
      composite += "dbl, ";
    } else {
      composite += "flt, ";
    }

    if(SIMD) {
      composite += "1>";
    } else {
      composite += "0>";
    }

    return composite;
  }
};

void writeHeaders(
  std::ofstream& benchmarkFile
) {
  // Write headers
  std::vector<std::string> headers {
    "N",
    "E",
    "Dlib",
    "Eigen",
    "EigenSIMD"
  };

  for(unsigned i = 0; i < 2; ++i) {
    benchmarkFile << "\"" << headers.at(i) << "\", ";
  }

  for(unsigned i = 2; i < headers.size(); ++i) {
    benchmarkFile << "\"" << headers.at(i) << "\", \"" << headers.at(i)
      << " sigma\"";
    if(i != headers.size() - 1) {
      benchmarkFile << ", ";
    }
  }

  benchmarkFile << nl;
}


enum class Algorithm {
  All,
  Dlib,
  EigenDouble,
  EigenFloat,
  EigenSIMDDouble,
  EigenSIMDFloat
};

void benchmark(
  const boost::filesystem::path& filePath,
  std::ofstream& benchmarkFile,
  const Algorithm algorithmChoice
) {
  using namespace molassembler;

  Molecule molecule = IO::read(
    filePath.string()
  );

  // Skip small molecules
  if(molecule.graph().N() < 10) {
    std::cout << "Skipping " << filePath.stem().string() << " since it's very small\n";
    return;
  }

  /* Embed, generate square bounds from distance bounds matrix and extract
   * chiral and dihedral constraints from spatial model
   */

  /* Timings */
  std::cout << std::fixed << std::setprecision(0);

  std::string nCount = "N = " + std::to_string(molecule.graph().N());

  std::cout
    << std::setw(6) << nCount
    << std::setw(10) << "Name"
    << std::setw(10) << "Rel. v"
    << std::setw(14) << "Successes"
    << std::setw(25) << "Time mu(sigma) / 1e-6s"
    << std::setw(25) << "Iter mu(sigma) / 1"
    << nl;

  std::vector<
    std::unique_ptr<TimingFunctor>
  > functors;

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::Dlib) {
    functors.emplace_back(
      std::make_unique<DlibFunctor>()
    );
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::EigenDouble) {
    functors.emplace_back(
      std::make_unique<
        EigenFunctor<4, double, false>
      >()
    );
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::EigenFloat) {
    functors.emplace_back(
      std::make_unique<
        EigenFunctor<4, float, false>
      >()
    );
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::EigenSIMDDouble) {
    functors.emplace_back(
      std::make_unique<
        EigenFunctor<4, double, true>
      >()
    );
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::EigenSIMDFloat) {
    functors.emplace_back(
      std::make_unique<
        EigenFunctor<4, float, true>
      >()
    );
  }

  auto results = timeFunctors<nExperiments>(molecule, functors);

  double smallestAverage = std::min_element(
    results.begin(),
    results.end(),
    [](const FunctorResults& a, const FunctorResults& b) -> bool {
      return a.timingAverage < b.timingAverage;
    }
  )->timingAverage;

  benchmarkFile
    << std::fixed << std::setprecision(0)
    << molecule.graph().N() << ", " << molecule.graph().B() << ", "
    << std::scientific << std::setprecision(6);

  for(unsigned i = 0; i < functors.size(); ++i) {
    const FunctorResults& functorResult = results.at(i);
    const auto& functorPtr = functors.at(i);

    std::string successCount = std::to_string(functorResult.count) + "/" + std::to_string(nExperiments);

    std::string timingStr = std::to_string(static_cast<int>(functorResult.timingAverage)) + "(" + std::to_string(static_cast<int>(functorResult.timingStddev)) + ")";

    std::string iterationStr = std::to_string(static_cast<int>(functorResult.iterationsAverage)) + "(" + std::to_string(static_cast<int>(functorResult.iterationsStddev)) + ")";

    std::cout
      << std::setw(16) << functorPtr->name()
      << std::setw(10) << (functorResult.timingAverage / smallestAverage)
      << std::setw(14) << successCount
      << std::setw(25) << timingStr
      << std::setw(25) << iterationStr
      << nl;

    benchmarkFile << functorResult.timingAverage << ", " << functorResult.timingStddev;

    if(i != functors.size() - 1) {
      benchmarkFile << ", ";
    } else {
      benchmarkFile << nl;
      std::cout << nl;
    }
  }
}

using namespace std::string_literals;
const std::string algorithmChoices =
  "  0 - All\n"
  "  1 - Dlib\n"
  "  2 - Eigen<4, dbl, 0>\n"
  "  3 - Eigen<4, flt, 0>\n"
  "  4 - Eigen<4, dbl, 1>\n"
  "  5 - Eigen<4, flt, 1>\n";

int main(int argc, char* argv[]) {
  using namespace molassembler;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("c", boost::program_options::value<unsigned>(), "Specify algorithm to benchmark")
    ("m", boost::program_options::value<std::string>(), "Path to MOLFiles to benchmark")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  // Program options
  if(options_variables_map.count("help") > 0) {
    std::cout << options_description << std::endl;
    return 0;
  }

  if(options_variables_map.count("m") == 0) {
    std::cout << "You have not specified any path to MOLFiles that could be used" << nl;
    return 0;
  }

  std::string molPath = options_variables_map["m"].as<std::string>();

  Algorithm choice = Algorithm::All;
  if(options_variables_map.count("c") > 0) {
    unsigned combination = options_variables_map["c"].as<unsigned>();

    if(combination > 2) {
      std::cout << "Specified algorithm is out of bounds. Valid choices are:" << nl
        << algorithmChoices;
      return 0;
    }

    choice = static_cast<Algorithm>(combination);
  }

  // Benchmark everything
  std::ofstream benchmarkFile ("refinement_timings.csv");
  writeHeaders(benchmarkFile);

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(molPath)
  ) {
    benchmark(currentFilePath, benchmarkFile, choice);
  }

  return 0;
}
