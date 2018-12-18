/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "boost/program_options.hpp"

#include "gor1/Gor1.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/DistanceGeometry/Gor1.h"
#include "molassembler/DistanceGeometry/ImplicitGraphBoost.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/IO.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

#include <chrono>
#include <iomanip>

using namespace Scine;
using namespace molassembler;

constexpr size_t nExperiments = 10;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

template<typename TimingCallable, size_t N>
std::pair<double, double> timeFunctor(
  const Molecule& molecule,
  const DistanceGeometry::SpatialModel::BoundsList& bounds,
  DistanceGeometry::Partiality partiality
) {
  using namespace std::chrono;

  TimingCallable functor;

  time_point<steady_clock> start, end;
  std::array<double, N> timings;
  for(unsigned n = 0; n < N; ++n) {
    start = steady_clock::now();
    functor(molecule, bounds, partiality);
    end = steady_clock::now();

    timings.at(n) = duration_cast<nanoseconds>(end - start).count();
  }

  auto average = temple::average(timings);

  return {
    average,
    temple::stddev(timings, average)
  };
}

template<class Graph>
struct Gor1Functor {
  Eigen::MatrixXd operator() (
    const Molecule& molecule,
    const DistanceGeometry::SpatialModel::BoundsList& boundsList,
    DistanceGeometry::Partiality partiality
  ) {

    Graph graph {molecule, boundsList};
    if(auto distanceMatrixResult = graph.makeDistanceMatrix(partiality)) {
      return distanceMatrixResult.value();
    }

    std::cout << "Gor1 Functor distance matrix generation failed!";
    return {};
  }
};

struct DBM_FW_Functor {
  Eigen::MatrixXd operator() (
    const Molecule& molecule,
    const DistanceGeometry::SpatialModel::BoundsList& boundsList,
    DistanceGeometry::Partiality partiality
  ) {
    DistanceGeometry::DistanceBoundsMatrix bounds {molecule, boundsList};

    bounds.smooth();

    auto distancesMatrixResult = bounds.makeDistanceMatrix(partiality);
    if(!distancesMatrixResult) {
      std::cout << "DBM-FW distance matrix generation failed!";
    }

    return bounds.access();
  }
};

void writeHeaders(
  std::ofstream& benchmarkFile
) {
  // Write headers
  std::vector<std::string> headers {
    "N",
    "E",
    "Floyd-Warshall & DBM",
    "Gor & ExplicitGraph",
    "Gor & ImplicitGraph"
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
  ImplicitGor,
  ExplicitGor,
  MatrixFW
};

void benchmark(
  const boost::filesystem::path& filePath,
  std::ofstream& benchmarkFile,
  Algorithm algorithmChoice,
  DistanceGeometry::Partiality partiality
) {
  using namespace molassembler;

  Molecule sampleMol = IO::read(
    filePath.string()
  );

  DistanceGeometry::SpatialModel spatialModel {sampleMol, DistanceGeometry::Configuration {}};

  const auto boundsList = spatialModel.makeBoundsList();

  /*
   * Can calculate shortest paths with either:
   * - Floyd-Warshall in DistanceBoundsMatrix
   * - Bellman-Ford (like in ExplicitGraph)
   * - Gor1 (graph-independent implementation
   *
   * On either:
   * - BGL ExplicitGraph type (can access from ExplicitGraph)
   * - SPG type
   *
   * Time and compare correctness:
   * - DistanceBoundsMatrix  + Floyd-Warshall (currently ONLY correct impl)
   * - ExplicitGraph BGL graph + Bellman-Ford
   * - ExplicitGraph BGL graph + Gor1
   * - SPG w/out implicit   + Bellman-Ford
   * - SPG w/out implicit   + Gor1
   *
   * The bottom four combinations should yield equal shortest paths distances,
   * although they are not consistent with the triangle inequalities. Only DBM
   * + FW is _correct_ so far.
   */

  /* Timings */
  std::cout << std::fixed << std::setprecision(0);
  std::cout << std::setw(14) << "N" << std::setw(8) << sampleMol.graph().N() << nl;

  std::vector<std::string> names;
  std::vector<double> times;
  std::vector<double> deviations;

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::MatrixFW) {
    auto timings = timeFunctor<
      DBM_FW_Functor,
      nExperiments
    >(sampleMol, boundsList, partiality);

    names.emplace_back("DBM FW");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::ExplicitGor) {
    auto timings = timeFunctor<
      Gor1Functor<DistanceGeometry::ExplicitGraph>,
      nExperiments
    >(sampleMol, boundsList, partiality);

    names.emplace_back("Explicit Gor");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::ImplicitGor) {
    auto timings = timeFunctor<
      Gor1Functor<DistanceGeometry::ImplicitGraph>,
      nExperiments
    >(sampleMol, boundsList, partiality);

    names.emplace_back("Implicit Gor");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  auto smallest = *std::min_element(times.begin(), times.end());

  benchmarkFile
    << std::fixed << std::setprecision(0)
    << sampleMol.graph().N() << ", " << sampleMol.graph().B() << ", "
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
  "  0 - Matrix Floyd-Warshall\n"
  "  1 - Gor1 with ImplicitGraph\n"
  "  2 - Gor1 with ExplicitGraph\n";

const std::string partialityChoices =
  "  0 - Four-Atom Metrization\n"
  "  1 - 10% Metrization\n"
  "  2 - All (default)\n";

int main(int argc, char* argv[]) {
  using namespace molassembler;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("c", boost::program_options::value<unsigned>(), "Specify algorithm / graph combination to benchmark")
    ("p", boost::program_options::value<unsigned>(), "Specify partiality")
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
      std::cout << "Specified algorithm / graph combination is out of bounds. Valid choices are:" << nl
        << algorithmChoices;
      return 0;
    }

    choice = static_cast<Algorithm>(combination);
  }

  DistanceGeometry::Partiality partiality = DistanceGeometry::Partiality::All;
  if(options_variables_map.count("p") > 0) {
    unsigned index =  options_variables_map["p"].as<unsigned>();

    if(index > 2) {
      std::cout << "Specified metrization option is out of bounds. Valid choices are:" << nl
        << partialityChoices;
      return 0;
    }

    partiality = static_cast<DistanceGeometry::Partiality>(index);
  }

  // Benchmark everything
  std::ofstream benchmarkFile ("graph_timings.csv");
  writeHeaders(benchmarkFile);

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(molPath)
  ) {
    benchmark(currentFilePath, benchmarkFile, choice, partiality);
  }

  return 0;
}
