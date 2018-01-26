#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include "template_magic/Containers.h"
#include "DistanceGeometry/ImplicitGraphBoost.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/ExplicitGraph.h"
#include "template_magic/Numeric.h"
#include "IO.h"

#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "Graph/Gor1.h"
#include "DistanceGeometry/Gor1.h"

#include <chrono>

/* TODO
 * Full algorithm execution needs to be compared, including graph setup!
 */

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
    "Gor & LG",
    "Gor & SPG"
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


void addAll(
  const boost::filesystem::path& filePath, 
  std::ofstream& benchmarkFile,
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

  DistanceGeometry::ExplicitGraph limits {sampleMol, boundsList};
  DistanceGeometry::DistanceBoundsMatrix boundsMatrix {sampleMol, boundsList};

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
  auto DBM_FW_timings = timeFunctor<
    DBM_FW_Functor,
    nExperiments
  >(sampleMol, boundsList);

  auto LG_Gor_timings = timeFunctor<
    Gor1Functor<DistanceGeometry::ExplicitGraph>,
    nExperiments
  >(sampleMol, boundsList);

  auto SPG_Gor_timings = timeFunctor<
    Gor1Functor<DistanceGeometry::ImplicitGraph>,
    nExperiments
  >(sampleMol, boundsList);

  auto smallest = std::min({
    DBM_FW_timings.first,
    LG_Gor_timings.first,
    SPG_Gor_timings.first
  });

  std::cout << std::fixed << std::setprecision(0);
  std::cout << std::setw(10) << "N" << std::setw(10) << sampleMol.numAtoms() << nl
    << std::setw(10) << "DBM FW" << std::setw(10) << DBM_FW_timings.first / smallest << nl
    << std::setw(10) << "LG Gor" << std::setw(10) << LG_Gor_timings.first / smallest << nl
    << std::setw(10) << "SPG Gor" << std::setw(10) << SPG_Gor_timings.first / smallest << nl << nl;

  /* Write results */
  benchmarkFile 
    << std::fixed << std::setprecision(0)
    << sampleMol.numAtoms() << ", " << sampleMol.numBonds() << ", "
    << std::scientific << std::setprecision(6)
    << DBM_FW_timings.first << ", "  << DBM_FW_timings.second << ", "
    << LG_Gor_timings.first << ", "  << LG_Gor_timings.second << ", "
    << SPG_Gor_timings.first << ", "  << SPG_Gor_timings.second << nl;
}

void DBM_FW(
  const boost::filesystem::path& filePath, 
  std::ofstream& benchmarkFile,
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

  auto DBM_FW_timings = timeFunctor<
    DBM_FW_Functor,
    nExperiments
  >(sampleMol, spatialModel.makeBoundList());

  benchmarkFile 
    << std::fixed << std::setprecision(0)
    << sampleMol.numAtoms() << ", " << sampleMol.numBonds() << ", "
    << std::scientific << std::setprecision(6)
    << DBM_FW_timings.first << ", "  << DBM_FW_timings.second << nl;
}

void SPG_Gor(
  const boost::filesystem::path& filePath, 
  std::ofstream& benchmarkFile,
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

  Gor1Functor<DistanceGeometry::ImplicitGraph> functor;

  auto bounds = functor(sampleMol, spatialModel.makeBoundList());

  /*auto SPG_Gor_timings = timeFunctor<
    Gor1Functor<MakeSPGFunctor>,
    nExperiments
  >(spatialModel);

  benchmarkFile 
    << std::fixed << std::setprecision(0)
    << sampleMol.numAtoms() << ", " << sampleMol.numBonds() << ", "
    << std::scientific << std::setprecision(6)
    << SPG_Gor_timings.first << ", "  << SPG_Gor_timings.second << nl;*/

  std::cout << bounds << nl;
}

using namespace std::string_literals;
const std::string algorithmChoices = 
  "  0 - Floyd-Warshall with DistanceBoundsMatrix\n"s
  "  1 - Bellman-Ford with ExplicitGraph\n"s
  "  2 - Gor1 with ExplicitGraph\n"s
  "  3 - Bellman-Ford with ImplicitGraph\n"s
  "  4 - Gor1 with ImplicitGraph\n"s;

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

  if(options_variables_map.count("c")) {
    unsigned combination = options_variables_map["c"].as<unsigned>();
    if(combination > 4) {
      std::cout << "Specified algorithm / graph combination is out of bounds. Valid choices are:" << nl
        << algorithmChoices;
      return 0;
    }

    IO::MOLFileHandler molHandler;
    std::ofstream benchmarkFile ("graph_timings.csv");
    writeHeaders(benchmarkFile);

    boost::filesystem::path filesPath("../tests/mol_files/ranking_tree_molecules");
    boost::filesystem::recursive_directory_iterator end;

    for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
      const boost::filesystem::path currentFilePath = *i;

      if(combination == 0) {
        DBM_FW(currentFilePath, benchmarkFile, molHandler);
      }  else if(combination == 4) {
        SPG_Gor(currentFilePath, benchmarkFile, molHandler);
      }
    }
  } else {
    // Benchmark everything
    IO::MOLFileHandler molHandler;
    std::ofstream benchmarkFile ("graph_timings.csv");
    writeHeaders(benchmarkFile);

    boost::filesystem::path filesPath("../tests/mol_files/ranking_tree_molecules");
    boost::filesystem::recursive_directory_iterator end;

    for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
      const boost::filesystem::path currentFilePath = *i;

      addAll(currentFilePath, benchmarkFile, molHandler);
    }
  }

  return 0;
}
