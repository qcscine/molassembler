#ifndef INCLUDE_ANALYSIS_HELPERS_H
#define INCLUDE_ANALYSIS_HELPERS_H

#include "DistanceGeometry/generateConformation.h"
#include "template_magic/Enumerate.h"

/*! @file
 *
 * Contains some helper functions to write DG step raytrace and analysis files.
 */

namespace MoleculeManip {

namespace AnalysisHelpers {

namespace detail {

inline std::string mapIndexToChar(const unsigned& index) {
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

} // namespace detail

void writeDGPOVandProgressFiles(
  const Molecule& mol,
  const std::string& baseFilename,
  const DistanceGeometry::detail::RefinementData& refinementData
) {
  const unsigned dimensionality = 4;

  // Write the POV-Ray files and progress files
  std::string progressFilename = baseFilename + "-progress.csv"s;
  std::ofstream progressFile (progressFilename);

  progressFile << std::scientific;

  for(const auto& refinementEnumPair : enumerate(refinementData.steps)) {

    const auto& stepData = refinementEnumPair.value;
    assert(stepData.positions.size() % dimensionality == 0);
    const unsigned N = stepData.positions.size() / dimensionality;

    progressFile 
      << stepData.distanceError << "," 
      << stepData.chiralError << "," 
      << stepData.fourthDimError << "," 
      << dlib::length(stepData.gradient) << "," 
      << static_cast<unsigned>(stepData.compress) << "," 
      << stepData.proportionCorrectChiralityConstraints << "\n";

    // Write the POV file for this step
    std::stringstream filename;
    filename << baseFilename << "-"
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
      outStream << stepData.positions(dimensionality * i) << ", ";
      outStream << stepData.positions(dimensionality * i + 1) << ", ";
      outStream << stepData.positions(dimensionality * i + 2) << ">;\n";
    }
    outStream << "\n";

    // Atoms
    for(unsigned i = 0; i < N; i++) {
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
    for(const auto& edgePair : mol.getEdges()) {
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
}

// Shortcut for simple symmetry molecules
void writeDGPOVandProgressFiles(
  const Molecule& mol,
  const Symmetry::Name& symmetryName,
  const unsigned& structNum,
  const DistanceGeometry::detail::RefinementData& refinementData
) {
  std::string baseName = Symmetry::spaceFreeName(symmetryName)
    + "-"s + std::to_string(structNum);

  writeDGPOVandProgressFiles(
    mol,
    baseName,
    refinementData
  );

}


} // namespace AnalysisHelpers

} // namespace MoleculeManip

#endif
