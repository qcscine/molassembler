#define BOOST_SYSTEM_NO_DEPRECATED
#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "temple/VectorView.h"
#include "temple/constexpr/ConsecutiveCompare.h"

#include "chemical_symmetries/DynamicProperties.h"

#include <iostream>
#include <fstream>
#include <iomanip>

const std::array<unsigned, 4> distortionColumns {{8, 8, 8, 30}};
const std::array<unsigned, 3> symmetryColumns {{5, 5, 25}};
const std::array<unsigned, 4> ambiguityColumns {{10, 25, 25, 4}};

std::ostream& nl(std::ostream& out) {
  out << '\n';
  return out;
}

std::string condense(const std::vector<unsigned>& indexVector) {
  using namespace std::string_literals;
  return "{"s + temple::condenseIterable(indexVector) + "}"s;
}

std::ostream& operator << (
  std::ostream& out,
  const Symmetry::properties::DistortionInfo& distortion
) {
  out << std::setw(distortionColumns[0]) << distortion.angularDistortion
    << std::setw(distortionColumns[1]) << distortion.chiralDistortion
    << std::setw(distortionColumns[2]) << (distortion.angularDistortion + distortion.chiralDistortion)
    << std::setw(distortionColumns[3]) << condense(distortion.indexMapping);

  return out;
}

void printMappingsHeader() {
  std::cout << std::setw(distortionColumns[0]) << "Angular"
    << std::setw(distortionColumns[1]) << "Chiral"
    << std::setw(distortionColumns[2]) << "Total"
    << std::setw(distortionColumns[3]) << "Mapping" << nl;
}

void printPermissibleSymmetries() {
  std::cout << std::setw(symmetryColumns[0]) << "Idx"
    << std::setw(symmetryColumns[1]) << "Size"
    << std::setw(symmetryColumns[2]) << "Name"
    << nl;

  for(unsigned i = 0; i < Symmetry::allNames.size(); i++) {
    std::cout << std::setw(symmetryColumns[0]) << i 
      << std::setw(symmetryColumns[1]) << Symmetry::size(Symmetry::allNames.at(i))
      << std::setw(symmetryColumns[2]) << Symmetry::name(Symmetry::allNames.at(i)) 
      << nl;
  }

  std::cout << std::endl;
}

void writeDistortions(
  const std::vector<Symmetry::properties::DistortionInfo>& distortions
) {
  std::cout << std::fixed << std::setprecision(2);

  auto sortedView = temple::sort(
    distortions,
    [](const auto& a, const auto& b) -> bool {
      return temple::consecutiveCompareSmaller(
        b.angularDistortion,
        a.angularDistortion,
        b.chiralDistortion,
        a.chiralDistortion
      );
    }
  );

  for(const auto& distortion : sortedView) {
    std::cout << distortion << nl;
  }
}

double calculateAmbiguity(
  const std::vector<Symmetry::properties::DistortionInfo>& distortions
) {
  /* Some measure between 0 and 1 that indicates how ambiguous choosing the
   * lowest mapping is.
   */

  auto sortByTotalView = temple::sort(
    distortions,
    [](const auto& a, const auto& b) -> bool {
      return (
        (a.angularDistortion + a.chiralDistortion)
        < (b.angularDistortion + b.chiralDistortion)
      );
    }
  );

  if(distortions.size() <= 1) {
    return 0;
  } 

  auto firstValue = (
    sortByTotalView.at(0).angularDistortion 
    + sortByTotalView.at(0).chiralDistortion
  );
  auto secondValue = (
    sortByTotalView.at(1).angularDistortion 
    + sortByTotalView.at(1).chiralDistortion
  );

  if(std::fabs(secondValue - firstValue) < 1e-10) {
    return 1;
  } 

  return firstValue / secondValue;
}

struct AmbiguityEntry {
  double ambiguity;
  Symmetry::Name source, target;
  boost::optional<unsigned> deletedIndex;

  AmbiguityEntry(
    const double& ambiguity,
    const Symmetry::Name& source,
    const Symmetry::Name& target,
    const boost::optional<unsigned>& deletedIndex = boost::none
  ) : ambiguity(ambiguity), source(source), target(target), deletedIndex(deletedIndex) {};
};

int main(int argc, char* argv[]) {
/* Set program options from command-line arguments */
  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("s", boost::program_options::value<unsigned>(), "Source symmetry index")
    ("t", boost::program_options::value<unsigned>(), "Target symmetry index")
    ("a", boost::program_options::value<bool>(), "Calculate ambiguity of lowest mappings")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);  

  // Manage the results
  if(options_variables_map.count("help")) {
    std::cout << options_description << nl
      << "The following symmetry indices are permissible:" << nl;
    printPermissibleSymmetries();
    return 0;
  }

  if(options_variables_map.count("s") && options_variables_map.count("t")) {
    unsigned sourceSymmetryArg = options_variables_map["s"].as<unsigned>();
    unsigned targetSymmetryArg = options_variables_map["t"].as<unsigned>();

    if(
      sourceSymmetryArg >= Symmetry::allNames.size()
      || targetSymmetryArg >= Symmetry::allNames.size()
    ) {
      std::cout << "Specified symmetry out of bounds. Valid symmetries:" << nl << nl;
      printPermissibleSymmetries();
      return 1;
    }

    if(sourceSymmetryArg == targetSymmetryArg) {
      std::cout << "The source and target symmetry may not be identical." << nl << nl;
      printPermissibleSymmetries();
      return 1;
    }

    Symmetry::Name sourceSymmetry(Symmetry::allNames.at(sourceSymmetryArg)),
                   targetSymmetry(Symmetry::allNames.at(targetSymmetryArg));

    int diff = (
      static_cast<int>(Symmetry::size(targetSymmetry))
      - static_cast<int>(Symmetry::size(sourceSymmetry))
    );

    if(std::abs(diff) > 1) {
      std::cout << "The selected symmetries must be adjacent in size!" << nl << nl;
      printPermissibleSymmetries();
      return 1;
    }

    if(diff == 1 || diff == 0) {
      auto distortions = Symmetry::properties::symmetryTransitionMappings(
        sourceSymmetry,
        targetSymmetry
      );

      printMappingsHeader();
      writeDistortions(distortions);
    } else {
      for(unsigned i = 0; i < Symmetry::size(sourceSymmetry); ++i) {
        auto distortions = Symmetry::properties::ligandLossTransitionMappings(
          sourceSymmetry,
          targetSymmetry,
          i
        );

        printMappingsHeader();
        writeDistortions(distortions);
        std::cout << nl << nl;
      }
    }
  }

  if(options_variables_map.count("a")) {
    std::vector<AmbiguityEntry> ambiguities;

    for(const auto& sourceSymmetry : Symmetry::allNames) {
      for(const auto& targetSymmetry : Symmetry::allNames) {
        if(sourceSymmetry == targetSymmetry) {
          // Skip identity mapping
          continue;
        }

        int diff = (
          static_cast<int>(Symmetry::size(targetSymmetry))
          - static_cast<int>(Symmetry::size(sourceSymmetry))
        );

        if(diff == 1 || diff == 0) {
          auto distortions = Symmetry::properties::symmetryTransitionMappings(
            sourceSymmetry,
            targetSymmetry
          );

          auto ambiguity = calculateAmbiguity(distortions);

          if(0.0 < ambiguity && ambiguity < 1.0) {
            ambiguities.emplace_back(
              ambiguity,
              sourceSymmetry,
              targetSymmetry
            );
          }
        } else if(diff == -1) {
          for(unsigned i = 0; i < Symmetry::size(sourceSymmetry); ++i) {
            auto distortions = Symmetry::properties::ligandLossTransitionMappings(
              sourceSymmetry,
              targetSymmetry,
              i
            );

            auto ambiguity = calculateAmbiguity(distortions);

            if(0.0 < ambiguity && ambiguity < 1.0) {
              ambiguities.emplace_back(
                ambiguity,
                sourceSymmetry,
                targetSymmetry,
                i
              );
            }
          }
        }
      }
    }

    auto sortedView = temple::sort(
      ambiguities,
      [](const auto& a, const auto& b) -> bool {
        return a.ambiguity < b.ambiguity;
      }
    );

    std::cout 
      << "Ambiguity tries to quantify how bad choosing the index mapping with the" << nl
      << "lowest total distortion is over considering the next best index mapping." << nl
      << "Zero indicates that the choice is unambiguous, one that the choice is " << nl
      << "completely ambiguous. Ambiguity values excluding zero and one are shown" << nl 
      << "(both are common). Idx is the index that is deleted when the target" << nl 
      << "symmetry is smaller than the source symmetry." << nl << nl;


    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(ambiguityColumns[0]) << "Ambiguity"
      << std::setw(ambiguityColumns[1]) << "Source"
      << std::setw(ambiguityColumns[2]) << "Target"
      << std::setw(ambiguityColumns[3]) << "Idx" << nl;

    for(const auto& entry : sortedView) {
      std::cout << std::setw(ambiguityColumns[0]) << entry.ambiguity
        << std::setw(ambiguityColumns[1]) << Symmetry::name(entry.source) 
        << std::setw(ambiguityColumns[2]) << Symmetry::name(entry.target)
        << std::setw(ambiguityColumns[3]);

      if(entry.deletedIndex) {
        std::cout << entry.deletedIndex.value();
      } else {
        std::cout << " ";
      }
      
      std::cout << nl;
    }
  }

  return 0;
}
