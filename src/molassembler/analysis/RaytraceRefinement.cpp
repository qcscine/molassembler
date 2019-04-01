/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "temple/Adaptors/Enumerate.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"
#include "molassembler/DistanceGeometry/DlibRefinement.h"
#include "molassembler/IO.h"
#include "molassembler/Log.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace Scine;
using namespace molassembler;

namespace detail {

inline std::string mapIndexToChar(const unsigned index) {
  std::string result {
    static_cast<char>(
      'A' + (
        (index / 25) % 25
      )
    )
  };

  result += static_cast<char>(
    'A' + (index % 25)
  );

  return result;
}

void writePOVFile(
  const Molecule& mol,
  const std::string& baseFilename,
  const unsigned structureIndex,
  const DistanceGeometry::RefinementStepData& stepData,
  const std::vector<DistanceGeometry::ChiralityConstraint>& constraints
) {
  const unsigned dimensionality = 4;
  assert(stepData.positions.size() % dimensionality == 0);
  const unsigned N = stepData.positions.size() / dimensionality;
  assert(N == mol.graph().N());

  /* Write the POV file for this step */
  std::stringstream filename;
  filename << baseFilename << "-"
    << std::setfill('0') << std::setw(3) << structureIndex
    << ".pov";

  std::ofstream outStream(
    filename.str()
  );

  outStream << "#version 3.7;\n"
    << "#include \"scene.inc\"\n\n";

  // Define atom names with positions
  for(unsigned i = 0; i < N; ++i) {
    outStream << "#declare " << detail::mapIndexToChar(i) <<  " = <";
    outStream << std::fixed << std::setprecision(4);
    outStream << stepData.positions(dimensionality * i) << ", ";
    outStream << stepData.positions(dimensionality * i + 1) << ", ";
    outStream << stepData.positions(dimensionality * i + 2) << ">;\n";
  }
  outStream << "\n";

  // Atoms
  for(unsigned i = 0; i < N; ++i) {
    double fourthDimAbs = std::fabs(stepData.positions(dimensionality * i + 3));
    if(fourthDimAbs < 1e-4) {
      outStream << "Atom(" << detail::mapIndexToChar(i) << ")\n";
    } else {
      outStream << "Atom4D(" << detail::mapIndexToChar(i) << ", "
        << fourthDimAbs << ")\n";
    }
  }
  outStream << "\n";

  // Bonds
  for(const BondIndex& edge : boost::make_iterator_range(mol.graph().bonds())) {
    outStream << "Bond("
      << detail::mapIndexToChar(edge.first) << ","
      << detail::mapIndexToChar(edge.second)
      << ")\n";
  }
  outStream << "\n";

  // Tetrahedra
  if(!constraints.empty()) {
    auto writePosition = [&stepData](
      std::ostream& os,
      const std::vector<AtomIndex>& indices
    ) -> std::ostream& {
      // Calculate the average position
      auto averagePosition = DistanceGeometry::ErrorFunctionValue::getAveragePos3D(
        stepData.positions,
        indices
      );

      os << "<" << averagePosition.x() << ","
        << averagePosition.y() << ","
        << averagePosition.z() << ">";

      return os;
    };

    for(const auto& chiralityConstraint : constraints) {
      outStream << "TetrahedronHighlight(";
      writePosition(outStream, chiralityConstraint.sites[0]);
      outStream << ", ";
      writePosition(outStream, chiralityConstraint.sites[1]);
      outStream << ", ";
      writePosition(outStream, chiralityConstraint.sites[2]);
      outStream << ", ";
      writePosition(outStream, chiralityConstraint.sites[3]);
      outStream << ")\n";
    }
    outStream << "\n";
  }

  // Gradients
  for(unsigned i = 0; i < N; i++) {
    if(
      dlib::length(
        dlib::rowm(
          stepData.gradient,
          dlib::range(dimensionality * i, dimensionality * i + 2)
        )
      ) > 1e-5
    ) {
      outStream << "GradientVector(" << detail::mapIndexToChar(i) <<", <"
        << std::fixed << std::setprecision(4)
        << (-stepData.gradient(dimensionality * i)) << ", "
        << (-stepData.gradient(dimensionality * i + 1)) << ", "
        << (-stepData.gradient(dimensionality * i + 2)) << ", "
        << ">)\n";
    }
  }

  outStream.close();
}

} // namespace detail

void writeDGPOVandProgressFiles(
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
      << dlib::length(refinementStep.gradient) << ","
      << static_cast<unsigned>(refinementStep.compress) << ","
      << refinementStep.proportionCorrectChiralityConstraints << "\n";
  }

  progressFile.close();

  const unsigned maxPOVFiles = 100;

  if(refinementData.steps.size() > maxPOVFiles) {
    // Determine 100 roughly equispaced conformations to write to POV files
    double stepLength = static_cast<double>(refinementData.steps.size()) / maxPOVFiles;
    auto listIter = refinementData.steps.begin();
    unsigned currentIndex = 0;
    for(unsigned i = 0; i < maxPOVFiles; ++i) {
      unsigned targetIndex = std::floor(i * stepLength);
      assert(targetIndex >= currentIndex && targetIndex < refinementData.steps.size());
      std::advance(listIter, targetIndex - currentIndex);
      currentIndex = targetIndex;

      detail::writePOVFile(
        mol,
        baseFilename,
        i,
        *listIter,
        refinementData.constraints
      );
    }
  } else {
    for(const auto enumPair : temple::adaptors::enumerate(refinementData.steps)) {
      detail::writePOVFile(
        mol,
        baseFilename,
        enumPair.index,
        enumPair.value,
        refinementData.constraints
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
      "from_file,f",
      boost::program_options::value<std::string>(),
      "Read molecule to generate from file"
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
      "contributions,c",
      boost::program_options::bool_switch(&showFinalContributions),
      "Show the final contributions to the refinement error functions"
    )
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

  // Manage the results
  if(options_variables_map.count("help") > 0) {
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

  Log::particulars.insert(Log::Particulars::DGStructureAcceptanceFailures);

  if(showFinalContributions) {
    Log::particulars.insert(Log::Particulars::DGFinalErrorContributions);
  }

  unsigned nSteps = 10000;
  if(options_variables_map.count("steps") > 0) {
    nSteps = options_variables_map["steps"].as<unsigned>();
  }

/* Generating work */
  // Generate from file
  if(options_variables_map.count("from_file") == 1) {
    auto filename = options_variables_map["from_file"].as<std::string>();

    if(!boost::filesystem::exists(filename)) {
      std::cout << "The specified file could not be found!" << std::endl;
      return 0;
    }

    auto mol = IO::read(filename);

    std::cout << mol << std::endl;

    boost::filesystem::path filepath {filename};
    std::string filestem = filepath.stem().string();

    std::ofstream graphFile(filestem +  "-graph.dot");
    graphFile << mol.dumpGraphviz();
    graphFile.close();

    DistanceGeometry::Configuration DGConfiguration;
    DGConfiguration.partiality = metrizationOption;
    DGConfiguration.refinementStepLimit = nSteps;

#ifndef NDEBUG
    auto debugData = DistanceGeometry::debugRefinement(
      mol,
      nStructures,
      DGConfiguration
    );

    for(const auto& enumPair : temple::adaptors::enumerate(debugData)) {
      const auto& structNum = enumPair.index;
      const auto& refinementData = enumPair.value;

      std::string baseName = filestem + "-"s + std::to_string(structNum);

      writeDGPOVandProgressFiles(
        mol,
        baseName,
        refinementData
      );

      IO::write(
        filestem + "-"s + std::to_string(structNum) + "-last.mol"s,
        mol,
        DistanceGeometry::detail::convertToAngstromWrapper(
          refinementData.steps.back().positions
        )
      );
    }

    auto failures = temple::sum(
      temple::map(
        debugData,
        [](const auto& refinementData) -> unsigned {
          return static_cast<unsigned>(refinementData.isFailure);
        }
      )
    );

    if(failures > 0) {
      std::cout << "WARNING: " << failures << " refinements failed.\n";
    }
#else
    auto conformers = DistanceGeometry::run(
      mol,
      nStructures,
      DGConfiguration
    );

    unsigned i = 0;
    unsigned failures = 0;
    for(const auto& conformerResult : conformers) {
      if(conformerResult) {
        IO::write(
          filestem + "-"s + std::to_string(i) + "-last.mol"s,
          mol,
          conformerResult.value()
        );
      } else {
        std::cout << "Conformer " << i << " failed: " << conformerResult.error().message() << "\n";
        ++failures;
      }

      ++i;
    }

    std::cout << "WARNING: " << failures << " refinement(s) failed.\n";
#endif
  }
}
