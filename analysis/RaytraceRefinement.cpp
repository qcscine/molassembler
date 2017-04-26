#include "DistanceGeometry/generateConformation.h"
#include "BoundsFromSymmetry.h"
#include "IO.h"
#include "Log.h"
#include "template_magic/Enumerate.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

inline char mapIndexToChar(const unsigned& index) {
  return 'A' + index;
}

int main() {
  Log::level = Log::Level::None;
  Log::particulars = {Log::Particulars::DGRefinementChiralityNumericalDebugInfo};

  for(const auto& symmetryName : Symmetry::allNames) {
    // Make a space-free string from the name
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Make a molecule and generate an ensemble
    auto mol = DGDBM::symmetricMolecule(symmetryName);

    auto debugData = detail::debugDistanceGeometry(
      mol,
      3,
      MetrizationOption::off,
      false
    );

    // Write the POV-Ray files
    for(const auto& enumPair : enumerate(debugData.refinements)) {
      const auto& structNum = enumPair.index;
      const auto& refinementData = enumPair.value;

      for(const auto& refinementEnumPair : enumerate(refinementData.steps)) {
        const auto& stepData = refinementEnumPair.value;

        std::stringstream filename;
        filename << Symmetry::name(symmetryName) << "-"
          << structNum << "-"
          << std::setfill('0') << std::setw(3) << refinementEnumPair.index
          << ".pov";

        std::ofstream outStream(
          filename.str()
        );

        outStream << "#version 3.7;\n"
          << "#include \"scene.inc\"\n\n";

        const unsigned dimensionality = 4;

        assert(stepData.positions.size() % dimensionality == 0);
        const unsigned N = stepData.positions.size() / dimensionality;

        // Define atom names with positions
        for(unsigned i = 0; i < N; i++) {
          outStream << "#declare " << mapIndexToChar(i) <<  " = <";
          outStream << std::fixed << std::setprecision(4);
          outStream << stepData.positions[dimensionality * i] << ", ";
          outStream << stepData.positions[dimensionality * i + 1] << ", ";
          outStream << stepData.positions[dimensionality * i + 2] << ">;\n";
        }
        outStream << "\n";

        // Atoms
        for(unsigned i = 0; i < N; i++) {
          outStream << "Atom4D(" << mapIndexToChar(i) << ", " 
            << stepData.positions[dimensionality * i + 3] << ")\n";
        }
        outStream << "\n";

        // Bonds
        for(const auto& edgePair : mol.getAdjacencyList().getEdges()) {
          outStream << "Bond(" 
            << mapIndexToChar(edgePair.first.first) << ","
            << mapIndexToChar(edgePair.first.second) 
            << ")\n";
        }
        outStream << "\n";

        // Tetrahedra
        if(refinementData.constraints.size() > 0) {
          for(const auto& chiralityConstraint : refinementData.constraints) {
            outStream << "TetrahedronHighlight("
              << mapIndexToChar(chiralityConstraint.indices[0]) << ", "
              << mapIndexToChar(chiralityConstraint.indices[1]) << ", "
              << mapIndexToChar(chiralityConstraint.indices[2]) << ", "
              << mapIndexToChar(chiralityConstraint.indices[3]) 
              << ")\n";
          }
          outStream << "\n";
        }

        // Gradients
        for(unsigned i = 0; i < N; i++) {
          outStream << "GradientVector(" << mapIndexToChar(i) <<", <"
            << std::fixed << std::setprecision(4)
            << (-stepData.gradient[3 * i]) << ", "
            << (-stepData.gradient[3 * i + 1]) << ", "
            << (-stepData.gradient[3 * i + 2]) << ", "
            << ">)\n";
        }

        outStream.close();
      }
    }
  }
}
