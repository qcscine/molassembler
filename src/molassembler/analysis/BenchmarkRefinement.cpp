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
#include "molassembler/DistanceGeometry/EigenSIMDRefinement.h"
#include "molassembler/DistanceGeometry/DlibRefinement.h"
#include "molassembler/DistanceGeometry/DlibAdaptors.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"

#include "Utils/Optimizer/GradientBased/LBFGS.h"

#include "molassembler/IO.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Dense>

#include <chrono>
#include <iomanip>

using namespace Scine;
using namespace molassembler;

constexpr std::size_t nExperiments = 10;
constexpr std::size_t refinementStepLimit = 10000;
constexpr double refinementGradientTarget = 1e-5;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

template<typename TimingCallable, size_t N>
std::pair<double, double> timeFunctor(
  const Eigen::MatrixXd& squaredBounds,
  const std::vector<DistanceGeometry::ChiralityConstraint>& chiralConstraints,
  const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
  const Eigen::MatrixXd& positions
) {
  using namespace std::chrono;

  TimingCallable functor;

  time_point<steady_clock> start, end;
  std::array<double, N> timings;
  for(unsigned n = 0; n < N; ++n) {
    start = steady_clock::now();
    functor(
      squaredBounds,
      chiralConstraints,
      dihedralConstraints,
      positions
    );
    end = steady_clock::now();

    timings.at(n) = duration_cast<nanoseconds>(end - start).count();
  }

  auto average = temple::average(timings);

  return {
    average,
    temple::stddev(timings, average)
  };
}

struct DlibFunctor {
  boost::optional<AngstromWrapper> operator() (
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralityConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    Eigen::MatrixXd positions
  ) {
    using namespace DistanceGeometry;

    const unsigned N = positions.size() / 4;

    /* Transform positions to dlib space */
    ErrorFunctionValue::Vector dlibPositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
        Eigen::Map<Eigen::VectorXd>(
          positions.data(),
          positions.cols() * positions.rows()
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
      dlibAdaptors::IterationOrAllChiralitiesCorrectStrategy inversionStopStrategy {
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
      if(inversionStopStrategy.iterations >= refinementStepLimit) {
        return boost::none;
      }

      if(valueFunctor.proportionChiralConstraintsCorrectSign < 1.0) {
        return boost::none;
      }
    }

    dlibAdaptors::IterationOrGradientNormStopStrategy refinementStopStrategy {
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
    if(refinementStopStrategy.iterations >= refinementStepLimit) {
      return boost::none;
    }

    // Not all chirality constraints have the right sign
    if(valueFunctor.proportionChiralConstraintsCorrectSign < 1) {
      return boost::none;
    }

    return detail::convertToAngstromWrapper(dlibPositions);
  }
};

template<typename EigenRefinementType>
struct InversionOrIterLimitStop final : public Utils::GradientBasedCheck {
  const EigenRefinementType& refinementFunctorReference;

  InversionOrIterLimitStop(const EigenRefinementType& functor) : refinementFunctorReference(functor) {}

  virtual bool checkConvergence(
    const Eigen::VectorXd& /* parameters */,
    double /* value */,
    const Eigen::VectorXd& /* gradients */
  ) final {
    return refinementFunctorReference.proportionChiralConstraintsCorrectSign >= 1.0;
  }

  virtual bool checkMaxIterations(unsigned currentIteration) final {
    return currentIteration >= refinementStepLimit;
  }

  virtual void addSettingsDescriptors(Utils::UniversalSettings::DescriptorCollection& /* collection */) final {}
  virtual void applySettings(const Utils::Settings& /* s */) final {}
};

template<
  template<unsigned, typename> class EigenRefinementType,
  unsigned dimensionality,
  typename FloatType
>
boost::optional<AngstromWrapper> refine(
  const Eigen::MatrixXd& squaredBounds,
  const std::vector<DistanceGeometry::ChiralityConstraint>& chiralConstraints,
  const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
  Eigen::MatrixXd positions
) {
  using NTPRefinementType = EigenRefinementType<dimensionality, FloatType>;

  const unsigned N = positions.size() / dimensionality;
  /* Transfer positions into vector form */
  Eigen::VectorXd transformedPositions = Eigen::Map<Eigen::VectorXd>(
    positions.data(),
    positions.cols() * positions.rows()
  );

  NTPRefinementType refinementFunctor {
    squaredBounds,
    chiralConstraints,
    dihedralConstraints
  };

  double initiallyCorrectChiralConstraints = refinementFunctor.calculateProportionChiralConstraintsCorrectSign(transformedPositions);

  if(initiallyCorrectChiralConstraints < 0.5) {
    for(unsigned i = 0; i < N; ++i) {
      positions(dimensionality * i + 1) *= -1;
    }
    initiallyCorrectChiralConstraints = 1 - initiallyCorrectChiralConstraints;
  }

  Utils::LBFGS optimizer;

  /* First stage: Invert all chiral constraints */
  if(initiallyCorrectChiralConstraints < 1) {
    InversionOrIterLimitStop<NTPRefinementType> inversionChecker {refinementFunctor};
    const unsigned iterations = optimizer.optimize(
      transformedPositions,
      refinementFunctor,
      inversionChecker
    );

    if(iterations >= refinementStepLimit) {
      return boost::none;
    }

    if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1.0) {
      return boost::none;
    }
  }

  /* Second stage: Refine */
  refinementFunctor.compressFourthDimension = true;

  Utils::GradientBasedCheck gradientChecker;
  gradientChecker.maxIter = refinementStepLimit;
  gradientChecker.gradNorm = 1e-5;
  const unsigned iterations = optimizer.optimize(
    transformedPositions,
    refinementFunctor,
    gradientChecker
  );

  if(iterations >= refinementStepLimit) {
    return boost::none;
  }

  if(refinementFunctor.proportionChiralConstraintsCorrectSign < 1) {
    return boost::none;
  }

  Utils::PositionCollection finalPositions(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    finalPositions.row(i) = transformedPositions.segment<3>(dimensionality * i);
  }

  return AngstromWrapper {std::move(finalPositions), LengthUnit::Angstrom};
}

struct EigenFunctor {
  boost::optional<AngstromWrapper> operator() (
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralityConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    Eigen::MatrixXd positions
  ) {
    return refine<
      DistanceGeometry::EigenRefinementProblem,
      4,
      double
    >(
      squaredBounds,
      chiralConstraints,
      dihedralConstraints,
      std::move(positions)
    );
  }
};

struct EigenSIMDFunctor {
  boost::optional<AngstromWrapper> operator() (
    const Eigen::MatrixXd& squaredBounds,
    const std::vector<DistanceGeometry::ChiralityConstraint>& chiralConstraints,
    const std::vector<DistanceGeometry::DihedralConstraint>& dihedralConstraints,
    Eigen::MatrixXd positions
  ) {
    return refine<
      DistanceGeometry::EigenSIMDRefinementProblem,
      4,
      double
    >(
      squaredBounds,
      chiralConstraints,
      dihedralConstraints,
      std::move(positions)
    );
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
  Eigen,
  EigenSIMD
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

  DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

  const auto boundsList = spatialModel.makeBoundsList();

  const auto chiralConstraints = spatialModel.getChiralityConstraints();
  const auto dihedralConstraints = spatialModel.getDihedralConstraints();

  DistanceGeometry::ExplicitGraph explicitGraph {
    molecule,
    boundsList
  };

  auto distancesMatrixResult = explicitGraph.makeDistanceMatrix(randomnessEngine());
  if(!distancesMatrixResult) {
    std::cout << distancesMatrixResult.error().message() << "\n";
    return;
  }

  auto distanceBoundsResult = explicitGraph.makeDistanceBounds();
  if(!distanceBoundsResult) {
    std::cout << "Failure in distance bounds matrix construction: "
      << distanceBoundsResult.error().message() << "\n";
    return;
  }

  DistanceGeometry::DistanceBoundsMatrix distanceBounds {std::move(distanceBoundsResult.value())};

  Eigen::MatrixXd squaredBounds = distanceBounds.access().cwiseProduct(
    distanceBounds.access()
  );

  auto metricMatrix = DistanceGeometry::MetricMatrix(
    std::move(distancesMatrixResult.value())
  );

  auto embeddedPositions = metricMatrix.embed();

  /* Embed, generate square bounds from distance bounds matrix and extract
   * chiral and dihedral constraints from spatial model
   */

  /* Timings */
  std::cout << std::fixed << std::setprecision(0);
  std::cout << std::setw(14) << "N" << std::setw(8) << molecule.graph().N() << nl;

  std::vector<std::string> names;
  std::vector<double> times;
  std::vector<double> deviations;

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::Dlib) {
    auto timings = timeFunctor<
      DlibFunctor,
      nExperiments
    >(squaredBounds, chiralConstraints, dihedralConstraints, embeddedPositions);

    names.emplace_back("Dlib");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::Eigen) {
    auto timings = timeFunctor<
      EigenFunctor,
      nExperiments
    >(squaredBounds, chiralConstraints, dihedralConstraints, embeddedPositions);

    names.emplace_back("Eigen");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::EigenSIMD) {
    auto timings = timeFunctor<
      EigenSIMDFunctor,
      nExperiments
    >(squaredBounds, chiralConstraints, dihedralConstraints, embeddedPositions);

    names.emplace_back("EigenSIMD");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  auto smallest = *std::min_element(times.begin(), times.end());

  benchmarkFile
    << std::fixed << std::setprecision(0)
    << molecule.graph().N() << ", " << molecule.graph().B() << ", "
    << std::scientific << std::setprecision(6);

  for(unsigned i = 0; i < names.size(); ++i) {
    std::cout
      << std::setw(14) << names[i]
      << std::setw(8) << (times[i] / smallest)
      << std::setw(14) << times[i]
      << std::setw(14) << deviations[i]
      << nl;

    benchmarkFile << times[i] << ", " << deviations[i];

    if(i != names.size() - 1) {
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
  "  2 - Eigen\n"
  "  3 - Eigen with SIMD\n";

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
