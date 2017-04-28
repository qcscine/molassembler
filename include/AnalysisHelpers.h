#ifndef INCLUDE_ANALYSIS_HELPERS_H
#define INCLUDE_ANALYSIS_HELPERS_H

#include "DistanceGeometry/generateConformation.h"
#include "template_magic/Enumerate.h"

namespace MoleculeManip {

namespace AnalysisHelpers {

namespace detail {

inline char mapIndexToChar(const unsigned& index) {
  return 'A' + index;
}

} // namespace detail

void writeDGPOVandProgressFiles(
  const Molecule& mol,
  const Symmetry::Name& symmetryName,
  const unsigned& structNum,
  const DistanceGeometry::detail::RefinementData& refinementData
) {
  const unsigned dimensionality = 4;

  // Make a space-free string from the name
  std::string spaceFreeName = Symmetry::name(symmetryName);
  std::replace(
    spaceFreeName.begin(),
    spaceFreeName.end(),
    ' ',
    '-'
  );

  // Write the POV-Ray files and progress files
  std::string progressFilename = spaceFreeName + "-"s 
    + std::to_string(structNum) + "-progress.csv"s;
  std::ofstream progressFile (progressFilename);

  progressFile << std::scientific;

  for(const auto& refinementEnumPair : enumerate(refinementData.steps)) {

    const auto& stepData = refinementEnumPair.value;
    assert(stepData.positions.size() % dimensionality == 0);
    const unsigned N = stepData.positions.size() / dimensionality;

    // Collect sum of absolute 4D values
    double totalAbs4D = 0;
    for(unsigned i = 0; i < N; i++) {
      totalAbs4D += std::fabs(stepData.positions[4 * i + 3]);
    }

    progressFile << stepData.error << "," 
      << stepData.gradient.norm() << "," 
      << static_cast<unsigned>(stepData.compress) << "," 
      << totalAbs4D << "\n";

    std::stringstream filename;
    filename << spaceFreeName << "-"
      << structNum << "-"
      << std::setfill('0') << std::setw(3) << refinementEnumPair.index
      << ".pov";

    std::ofstream outStream(
      filename.str()
    );

    outStream << "#version 3.7;\n"
      << "#include \"scene.inc\"\n\n";

    // Define atom names with positions
    for(unsigned i = 0; i < N; i++) {
      outStream << "#declare " << detail::mapIndexToChar(i) <<  " = <";
      outStream << std::fixed << std::setprecision(4);
      outStream << stepData.positions[dimensionality * i] << ", ";
      outStream << stepData.positions[dimensionality * i + 1] << ", ";
      outStream << stepData.positions[dimensionality * i + 2] << ">;\n";
    }
    outStream << "\n";

    // Atoms
    for(unsigned i = 0; i < N; i++) {
      outStream << "Atom4D(" << detail::mapIndexToChar(i) << ", " 
        << std::fabs(stepData.positions[dimensionality * i + 3]) << ")\n";
    }
    outStream << "\n";

    // Bonds
    for(const auto& edgePair : mol.getAdjacencyList().getEdges()) {
      outStream << "Bond(" 
        << detail::mapIndexToChar(edgePair.first.first) << ","
        << detail::mapIndexToChar(edgePair.first.second) 
        << ")\n";
    }
    outStream << "\n";

    // Tetrahedra
    if(!refinementData.constraints.empty()) {
      for(const auto& chiralityConstraint : refinementData.constraints) {
        outStream << "TetrahedronHighlight("
          << detail::mapIndexToChar(chiralityConstraint.indices[0]) << ", "
          << detail::mapIndexToChar(chiralityConstraint.indices[1]) << ", "
          << detail::mapIndexToChar(chiralityConstraint.indices[2]) << ", "
          << detail::mapIndexToChar(chiralityConstraint.indices[3]) 
          << ")\n";
      }
      outStream << "\n";
    }

    // Gradients
    for(unsigned i = 0; i < N; i++) {
      outStream << "GradientVector(" << detail::mapIndexToChar(i) <<", <"
        << std::fixed << std::setprecision(4)
        << (-stepData.gradient[3 * i]) << ", "
        << (-stepData.gradient[3 * i + 1]) << ", "
        << (-stepData.gradient[3 * i + 2]) << ", "
        << ">)\n";
    }

    outStream.close();

  }
}

} // namespace AnalysisHelpers

} // namespace MoleculeManip

#endif
