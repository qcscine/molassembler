// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_ANALYSIS_HELPERS_H
#define INCLUDE_MOLASSEMBLER_ANALYSIS_HELPERS_H

#include "boost/range/iterator_range_core.hpp"
#include "chemical_symmetries/Symmetries.h"
#include "temple/Adaptors/Enumerate.h"

#include "molassembler/DistanceGeometry/ConformerGeneration.h"
#include "molassembler/DistanceGeometry/RefinementProblem.h"

/*! @file
 *
 * @brief Helper functions to write DG step raytrace and analysis files.
 */

namespace molassembler {

namespace AnalysisHelpers {

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
      auto averagePosition = DistanceGeometry::errfDetail::getAveragePos3D(
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

// Shortcut for simple symmetry molecules
void writeDGPOVandProgressFiles(
  const Molecule& mol,
  const Symmetry::Name& symmetryName,
  const unsigned structNum,
  const DistanceGeometry::RefinementData& refinementData
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

} // namespace molassembler

#endif
