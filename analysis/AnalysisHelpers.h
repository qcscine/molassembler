#include <iomanip>
#include <fstream>
#include <cassert>

#include "Molecule.h"
#include "IO.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "DistanceGeometry/DGRefinementProblem.h"

Delib::PositionCollection toPositionCollection(
  const Eigen::VectorXd& vectorizedPositions,
  const MoleculeManip::DistanceGeometry::EmbeddingOption& embedding
) {
  auto dimensionality = static_cast<unsigned>(embedding);

  assert(vectorizedPositions.size() % dimensionality == 0);

  Delib::PositionCollection positions;

  for(unsigned i = 0; i < vectorizedPositions.size() / dimensionality; i++) {
    positions.push_back(
      Delib::Position( // converting constructor
        vectorizedPositions.segment<3>(dimensionality * i)
      )
    );
  }

  return positions;
}

void writeFile(
  const bool& optimized,
  const std::string symmetryString,
  const unsigned& structNum,
  const Eigen::VectorXd& vectorizedPositions,
  const MoleculeManip::DistanceGeometry::EmbeddingOption& embedding
) {
  using namespace std::string_literals;

  auto dimensionality = static_cast<unsigned>(embedding);

  assert(vectorizedPositions.size() % dimensionality == 0);

  unsigned N = vectorizedPositions.size() / dimensionality;

  std::vector<std::string> elementNames = {
    "Ru",
    "H",
    "F",
    "Cl",
    "Br",
    "I",
    "N",
    "C",
    "O",
    "S",
    "P"
  };

  std::string filename;

  if(optimized) filename = "opt-"s;
  else filename = "gen-"s;

  filename += symmetryString + std::to_string(structNum) + ".xyz"s;

  /* From here it's a simplified XYZ Writer */
  std::ofstream outStream(filename.c_str());
  outStream << std::setprecision(7) << std::fixed;
  outStream << N << std::endl
    << "Energy = " << std::endl;

  for(unsigned i = 0; i < N; i++) {
    if(elementNames[i].size() == 1) outStream << elementNames[i] << " ";
    else outStream << elementNames[i];

    outStream << std::setw(13) << vectorizedPositions[dimensionality * i + 0];
    outStream << std::setw(13) << vectorizedPositions[dimensionality * i + 1];
    outStream << std::setw(13) << vectorizedPositions[dimensionality * i + 2];

    if(i != N - 1) outStream << std::endl;
  }

  outStream.close();
}

/* TODO there is no underlying molecule data, we need it, maybe from the 
 * distance bounds matrix generating function?
 */
void writeMOLFile(
  const MoleculeManip::Molecule& molecule,
  const bool& optimized,
  const std::string& symmetryString,
  const unsigned& structNum,
  const Eigen::VectorXd& vectorizedPositions,
  const MoleculeManip::DistanceGeometry::EmbeddingOption& embedding
) {
  using namespace std::string_literals;

  auto moleculeCopy = molecule;

  auto dimensionality = static_cast<unsigned>(embedding);

  assert(vectorizedPositions.size() % dimensionality == 0);

  std::vector<Delib::ElementType> elementTypes = {
    Delib::ElementType::Ru,
    Delib::ElementType::H,
    Delib::ElementType::F,
    Delib::ElementType::Cl,
    Delib::ElementType::Br,
    Delib::ElementType::I,
    Delib::ElementType::N,
    Delib::ElementType::C,
    Delib::ElementType::O,
    Delib::ElementType::S,
    Delib::ElementType::P
  };

  for(unsigned i = 0; i < molecule.getNumAtoms(); i++) {
    moleculeCopy.changeElementType(
      i,
      elementTypes.at(i)
    );
  }

  std::string filename;

  if(optimized) filename = "opt-"s;
  else filename = "gen-"s;

  filename += symmetryString + std::to_string(structNum) + ".mol"s;

  MoleculeManip::IO::MOLFileHandler fileHandler;
  fileHandler.writeSingle(
    filename,
    moleculeCopy,
    toPositionCollection(
      vectorizedPositions,
      embedding
    )
  );
}

void writeErrorValues(
  const std::string& symmetryString,
  const std::vector<double> errorValues
) {
  using namespace std::string_literals;

  std::ofstream outStream(symmetryString+".csv"s);

  for(const auto& value: errorValues) {
    outStream << value << std::endl;
  }

  outStream.close();
}

char mapIndexToChar(const unsigned& index) {
  return 'A' + index;
}

void writePOVRayFile(
  const std::string& symmetryString,
  const unsigned& structNum,
  const unsigned& iterations,
  const MoleculeManip::DistanceGeometry::DGRefinementProblem<double>& problem,
  const MoleculeManip::Molecule& molecule,
  const cppoptlib::ConjugatedGradientDescentSolver<
    MoleculeManip::DistanceGeometry::DGRefinementProblem<double>
  >::TVector& positions,
  const cppoptlib::ConjugatedGradientDescentSolver<
    MoleculeManip::DistanceGeometry::DGRefinementProblem<double>
  >::StepResultType& stepResult,
  const MoleculeManip::DistanceGeometry::EmbeddingOption& embedding
) {
  using namespace std::string_literals;

  std::stringstream filename;
  filename << symmetryString << "-"
    << structNum << "-"
    << std::setfill('0') << std::setw(3) << iterations
    << ".pov";

  std::ofstream outStream(
    filename.str()
  );

  outStream << "#version 3.7;\n"
    << "#include \"scene.inc\"\n\n";

  auto dimensionality = static_cast<unsigned>(embedding);

  assert(positions.size() % dimensionality == 0);
  const unsigned N = positions.size() / dimensionality;

  // Define atom names with positions
  for(unsigned i = 0; i < N; i++) {
    outStream << "#declare " << mapIndexToChar(i) <<  " = <";
    outStream << std::fixed << std::setprecision(4);
    outStream << positions[dimensionality * i] << ", ";
    outStream << positions[dimensionality * i + 1] << ", ";
    outStream << positions[dimensionality * i + 2] << ">;\n";
  }
  outStream << "\n";

  // Atoms
  if(dimensionality == 3) {
    for(unsigned i = 0; i < N; i++) {
      outStream << "Atom(" << mapIndexToChar(i) << ")\n";
    }
  } else { // 4D
    for(unsigned i = 0; i < N; i++) {
      outStream << "Atom4D(" << mapIndexToChar(i) << ", " 
        << positions[dimensionality * i + 3] << ")\n";
    }
  }
  outStream << "\n";

  // Bonds
  for(const auto& edgePair : molecule.getAdjacencyList().getEdges()) {
    outStream << "Bond(" 
      << mapIndexToChar(edgePair.first.first) << ","
      << mapIndexToChar(edgePair.first.second) 
      << ")\n";
  }
  outStream << "\n";

  // Tetrahedra
  if(problem.constraints.size() > 0) {
    for(const auto& chiralityConstraint : problem.constraints) {
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
      << stepResult.negativeGradient[3 * i] << ", "
      << stepResult.negativeGradient[3 * i + 1] << ", "
      << stepResult.negativeGradient[3 * i + 2] 
      << ">)\n";
  }

  outStream.close();
}
