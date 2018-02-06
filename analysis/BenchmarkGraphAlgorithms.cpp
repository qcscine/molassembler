#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include "template_magic/Containers.h"
#include "DistanceGeometry/ImplicitGraphBoost.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/ExplicitGraph.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "template_magic/Numeric.h"
#include "IO.h"

#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "Graph/Gor1.h"
#include "DistanceGeometry/Gor1.h"

#include <chrono>

using namespace MoleculeManip;

constexpr size_t nExperiments = 10;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

template<typename TimingCallable, size_t N>
std::pair<double, double> timeFunctor(
  const Molecule& molecule,
  const DistanceGeometry::MoleculeSpatialModel::BoundList& boundList
) {
  using namespace std::chrono;

  TimingCallable functor;

  time_point<steady_clock> start, end;
  std::array<double, N> timings;
  for(unsigned n = 0; n < N; ++n) {
    start = steady_clock::now();
    functor(molecule, boundList);
    end = steady_clock::now();

    timings.at(n) = duration_cast<nanoseconds>(end - start).count();
  }

  auto average = TemplateMagic::average(timings);

  return {
    average,
    TemplateMagic::stddev(timings, average)
  };
}

template<class Graph>
struct Gor1Functor {
  Eigen::MatrixXd operator() (
    const Molecule& molecule,
    const DistanceGeometry::ImplicitGraph::BoundList& boundsList
  ) {

    Graph graph {molecule, boundsList};
    return graph.makeDistanceMatrix();
  }
};

struct DBM_FW_Functor {
  Eigen::MatrixXd operator() (
    const Molecule& molecule,
    const DistanceGeometry::ImplicitGraph::BoundList& boundsList
  ) {
    DistanceGeometry::DistanceBoundsMatrix bounds {molecule, boundsList};

    bounds.smooth();

    bounds.makeDistanceMatrix();

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
  IO::MOLFileHandler& molHandler
) {
  using namespace MoleculeManip;

  Molecule sampleMol = molHandler.readSingle(
    filePath.string()
  );

  DistanceGeometry::MoleculeSpatialModel spatialModel {
    sampleMol,
    DistanceGeometry::MoleculeSpatialModel::DistanceMethod::UFFLike
  };

  const auto boundsList = spatialModel.makeBoundList();

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
  std::cout << std::setw(14) << "N" << std::setw(8) << sampleMol.numAtoms() << nl;

  std::vector<std::string> names;
  std::vector<double> times;
  std::vector<double> deviations;

  bool dontBreak = algorithmChoice == Algorithm::All;
  switch(algorithmChoice) {
    case Algorithm::All:
    case Algorithm::MatrixFW: {
      auto timings = timeFunctor<
        DBM_FW_Functor,
        nExperiments
      >(sampleMol, boundsList);

      names.push_back("DBM FW");
      times.push_back(timings.first);
      deviations.push_back(timings.second);

      if(!dontBreak) {
        break;
      } else {
        [[fallthrough]];
      }
    };

    case Algorithm::ImplicitGor: {
      auto timings = timeFunctor<
        Gor1Functor<DistanceGeometry::ImplicitGraph>,
        nExperiments
      >(sampleMol, boundsList);

      names.push_back("Implicit Gor");
      times.push_back(timings.first);
      deviations.push_back(timings.second);

      if(!dontBreak) {
        break;
      } else {
        [[fallthrough]];
      }
    };

    case Algorithm::ExplicitGor: {
      auto timings = timeFunctor<
        Gor1Functor<DistanceGeometry::ExplicitGraph>,
        nExperiments
      >(sampleMol, boundsList);

      names.push_back("Explicit Gor");
      times.push_back(timings.first);
      deviations.push_back(timings.second);
    };
  }

  auto smallest = *std::min_element(times.begin(), times.end());

  benchmarkFile 
    << std::fixed << std::setprecision(0)
    << sampleMol.numAtoms() << ", " << sampleMol.numBonds() << ", "
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
  "  0 - Matrix Floyd-Warshall\n"s
  "  1 - Gor1 with ImplicitGraph\n"s
  "  2 - Gor1 with ExplicitGraph\n"s;

int main(int argc, char* argv[]) {
  using namespace MoleculeManip;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("c", boost::program_options::value<unsigned>(), "Specify algorithm / graph combination to benchmark")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);  


  // Program options
  if(options_variables_map.count("help")) {
    std::cout << options_description << std::endl;
    return 0;
  }

  Algorithm choice = Algorithm::All;
  if(options_variables_map.count("c")) {
    unsigned combination = options_variables_map["c"].as<unsigned>();

    if(combination > 2) {
      std::cout << "Specified algorithm / graph combination is out of bounds. Valid choices are:" << nl
        << algorithmChoices;
      return 0;
    }

    choice = static_cast<Algorithm>(combination);
  }

  // Benchmark everything
  IO::MOLFileHandler molHandler;
  std::ofstream benchmarkFile ("graph_timings.csv");
  writeHeaders(benchmarkFile);

  boost::filesystem::path filesPath("../tests/mol_files/ranking_tree_molecules");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    benchmark(currentFilePath, benchmarkFile, choice, molHandler);
  }

  return 0;
}
