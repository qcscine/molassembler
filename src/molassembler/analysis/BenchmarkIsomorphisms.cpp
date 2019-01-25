/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
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
  const Molecule& a,
  const Molecule& b
) {
  using namespace std::chrono;

  TimingCallable functor;

  time_point<steady_clock> start, end;
  std::array<double, N> timings;
  for(unsigned n = 0; n < N; ++n) {
    start = steady_clock::now();
    functor(a, b);
    end = steady_clock::now();

    timings.at(n) = duration_cast<nanoseconds>(end - start).count();
  }

  auto average = temple::average(timings);

  return {
    average,
    temple::stddev(timings, average)
  };
}

struct BoostFunctor {
  bool operator() (const Molecule& a, const Molecule& b) {
    return a.modularCompare(
      b,
      temple::make_bitmask(AtomEnvironmentComponents::ElementTypes)
        | AtomEnvironmentComponents::BondOrders
        | AtomEnvironmentComponents::Symmetries
        | AtomEnvironmentComponents::Stereopermutations
    );
  }
};

struct TracesFunctor {
  bool operator() (const Molecule& a, const Molecule& b) {
    return a.trialModularCompare(
      b,
      temple::make_bitmask(AtomEnvironmentComponents::ElementTypes)
        | AtomEnvironmentComponents::BondOrders
        | AtomEnvironmentComponents::Symmetries
        | AtomEnvironmentComponents::Stereopermutations
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
    "Boost",
    "Traces"
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
  Boost,
  Traces
};

void benchmark(
  const boost::filesystem::path& filePath,
  std::ofstream& benchmarkFile,
  Algorithm algorithmChoice,
  bool canonicalize
) {

  Molecule molecule = IO::read(
    filePath.string()
  );

  // Skip small molecules
  if(molecule.graph().N() < 10) {
    std::cout << "Skipping " << filePath.stem().string() << " since it's very small\n";
    return;
  }

  // Generate an isomorphic molecule
  auto permutation = temple::iota<AtomIndex>(molecule.graph().N());
  temple::random::shuffle(permutation, randomnessEngine());
  Molecule isomorphic_molecule = molecule;
  isomorphic_molecule.applyPermutation(permutation);

  if(canonicalize) {
    molecule.canonicalize();
    isomorphic_molecule.canonicalize();
  }

  /* Timings */
  std::cout << std::fixed << std::setprecision(0);
  std::cout << std::setw(14) << "N" << std::setw(8) << molecule.graph().N() << nl;

  std::vector<std::string> names;
  std::vector<double> times;
  std::vector<double> deviations;

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::Boost) {
    auto timings = timeFunctor<
      BoostFunctor,
      nExperiments
    >(molecule, isomorphic_molecule);

    names.emplace_back("Boost");
    times.push_back(timings.first);
    deviations.push_back(timings.second);
  }

  if(algorithmChoice == Algorithm::All || algorithmChoice == Algorithm::Traces) {
    auto timings = timeFunctor<
      TracesFunctor,
      nExperiments
    >(molecule, isomorphic_molecule);

    names.emplace_back("Traces");
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

int main(int argc, char* argv[]) {

  bool canonicalize = false;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    (
      "canonicalize,c",
      boost::program_options::bool_switch(&canonicalize),
      "Canonicalize the molecules before isomorphism"
    )
    ("m", boost::program_options::value<std::string>(), "Path to MOLFiles to benchmark")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv).
    options(options_description).
    style(
      boost::program_options::command_line_style::unix_style
      | boost::program_options::command_line_style::allow_long_disguise
    ).run(),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);
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
  // Benchmark everything
  std::ofstream benchmarkFile ("isomorphism_timings.csv");
  writeHeaders(benchmarkFile);

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(molPath)
  ) {
    benchmark(currentFilePath, benchmarkFile, choice, canonicalize);
  }

  return 0;
}
