
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

#include "molassembler/IO.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Dense>

#include <chrono>
#include <iomanip>

using namespace Scine;
using namespace molassembler;

constexpr std::size_t nExperiments = 1000;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

struct TimingFunctor {
  virtual double value(
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    const Eigen::MatrixXd& positions,
    std::chrono::time_point<std::chrono::steady_clock>& start,
    std::chrono::time_point<std::chrono::steady_clock>& end
  ) = 0;

  virtual std::string name() = 0;
};

struct FunctorResults {
  unsigned count = 0;
  double timingAverage = 0;
  double timingStddev = 0;
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
    std::vector<double> results;
  };

  std::vector<Counter> counters(functors.size());

  for(unsigned n = 0; n < N; ++n) {
    // Prepare new data for the functors in every experiment
    DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

    const auto boundsList = spatialModel.makePairwiseBounds();

    const auto chiralConstraints = spatialModel.getChiralConstraints();
    const auto dihedralConstraints = spatialModel.getDihedralConstraints();

    DistanceGeometry::ExplicitGraph explicitGraph {
      molecule.graph().inner(),
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
      double result = functors.at(functorIndex)->value(
        squaredBounds,
        chiralConstraints,
        dihedralConstraints,
        embeddedPositions,
        start,
        end
      );

      counters.at(functorIndex).timings.push_back(duration_cast<nanoseconds>(end - start).count());
      counters.at(functorIndex).results.push_back(result);
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
      }

      return result;
    }
  );
}

struct DlibFunctor final : public TimingFunctor {
  double value (
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    const Eigen::MatrixXd& positions,
    std::chrono::time_point<std::chrono::steady_clock>& start,
    std::chrono::time_point<std::chrono::steady_clock>& end
  ) final {
    using namespace DistanceGeometry;

    Eigen::MatrixXd positionCopy = positions;

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

    double error = 0;

    start = std::chrono::steady_clock::now();
    for(unsigned i = 0; i < 10; ++i) {
      error += valueFunctor(dlibPositions);
      auto gradient = gradientFunctor(dlibPositions);
      error += gradient(0);
    }
    end = std::chrono::steady_clock::now();

    return error;
  }

  std::string name() final {
    return "Dlib";
  }
};

template<
  template<unsigned, typename, bool> class EigenRefinementType,
  unsigned dimensionality,
  typename FloatType,
  bool SIMD
>
double timeFunctionEvaluation(
  const Eigen::MatrixXd& squaredBounds,
  const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
  const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
  const Eigen::MatrixXd& positions,
  std::chrono::time_point<std::chrono::steady_clock>& start,
  std::chrono::time_point<std::chrono::steady_clock>& end
) {
  using NTPRefinementType = EigenRefinementType<dimensionality, FloatType, SIMD>;
  using PositionsType = typename NTPRefinementType::VectorType;

  Eigen::MatrixXd positionCopy = positions;

  /* Transfer positions into vector form */
  PositionsType transformedPositions = Eigen::Map<Eigen::VectorXd>(
    positionCopy.data(),
    positionCopy.cols() * positionCopy.rows()
  ).template cast<FloatType>();

  NTPRefinementType refinementFunctor {
    squaredBounds,
    chiralConstraints,
    dihedralConstraints
  };

  FloatType error = 0;
  PositionsType gradient(transformedPositions.size());
  gradient.setZero();

  start = std::chrono::steady_clock::now();
  for(unsigned i = 0; i < 10; ++i) {
    refinementFunctor(transformedPositions, error, gradient);
  }
  end = std::chrono::steady_clock::now();

  return error + gradient(0);
}

template<
  unsigned dimensionality,
  typename FloatType,
  bool SIMD
>
struct EigenFunctor final : public TimingFunctor {
  double value(
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    const Eigen::MatrixXd& positions,
    std::chrono::time_point<std::chrono::steady_clock>& start,
    std::chrono::time_point<std::chrono::steady_clock>& end
  ) final {
    return timeFunctionEvaluation<
      DistanceGeometry::EigenRefinementProblem,
      dimensionality,
      FloatType,
      SIMD
    >(
      squaredBounds,
      chiralConstraints,
      dihedralConstraints,
      positions,
      start,
      end
    );
  }

  std::string name() final {
    return DistanceGeometry::EigenRefinementProblem<dimensionality, FloatType, SIMD>::name();
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
    << std::setw(25) << "Time mu(sigma) / 1e-9s"
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

    std::cout
      << std::setw(16) << functorPtr->name()
      << std::setw(10) << (functorResult.timingAverage / smallestAverage)
      << std::setw(25) << timingStr
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
  "  1 - dlib\n"
  "  2 - Eigen<dimensionality=4, double, SIMD=false>\n"
  "  3 - Eigen<dimensionality=4, float, SIMD=false>\n"
  "  4 - Eigen<dimensionality=4, double, SIMD=true>\n"
  "  5 - Eigen<dimensionality4, float, SIMD=true>\n";

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
