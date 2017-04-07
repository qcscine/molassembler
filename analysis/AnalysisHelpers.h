#include <iomanip>
#include <fstream>
#include <cassert>

#include "Molecule.h"
#include "IO.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "DistanceGeometry/DGRefinementProblem.h"

Delib::PositionCollection toPositionCollection(
  const Eigen::VectorXd& vectorizedPositions
) {
  assert(vectorizedPositions.size() % 3 == 0);
  Delib::PositionCollection positions;

  for(unsigned i = 0; i < vectorizedPositions.size() / 3; i++) {
    positions.push_back(
      Delib::Position( // converting constructor
        vectorizedPositions.segment<3>(3 * i)
      )
    );
  }

  return positions;
}

void writeFile(
  const bool& optimized,
  const std::string symmetryString,
  const unsigned& structNum,
  const Eigen::VectorXd& vectorizedPositions
) {
  using namespace std::string_literals;

  unsigned N = vectorizedPositions.size() / 3;
  assert(3 * N == vectorizedPositions.size());

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

    outStream << std::setw(13) << vectorizedPositions[3 * i + 0];
    outStream << std::setw(13) << vectorizedPositions[3 * i + 1];
    outStream << std::setw(13) << vectorizedPositions[3 * i + 2];

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
  const Eigen::VectorXd& vectorizedPositions
) {
  using namespace std::string_literals;

  auto moleculeCopy = molecule;

  unsigned N = vectorizedPositions.size() / 3;
  assert(3 * N == vectorizedPositions.size());

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
      vectorizedPositions
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
  >::StepResultType& stepResult
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

  outStream << "#include \"a_geo_tetra1.inc\"\n\n";

  assert(positions.size() % 3 == 0);
  const unsigned N = positions.size() / 3;

  // Define atom names with positions
  for(unsigned i = 0; i < N; i++) {
    outStream << "#declare " << mapIndexToChar(i) <<  " = <";
    outStream << std::fixed << std::setprecision(4);
    outStream << positions[3 * i] << ", ";
    outStream << positions[3 * i + 1] << ", ";
    outStream << positions[3 * i + 2] << ">;\n";
  }
  outStream << "\n";

  // Atoms
  for(unsigned i = 0; i < N; i++) {
    outStream << "object { sphere{" << mapIndexToChar(i) 
      << ", 0.1} pigment {color myGray }}\n";
  }
  outStream << "\n";

  // Bonds
  for(const auto& edgePair : molecule.getAdjacencyList().getEdges()) {
    outStream << "object { cylinder{" 
      << mapIndexToChar(edgePair.first.first) << ","
      << mapIndexToChar(edgePair.first.second) 
      << ", 0.05} pigment { color myGray } }\n";
  }
  outStream << "\n";

  // Tetrahedra
  if(problem.constraints.size() > 0) {
    for(const auto& chiralityConstraint : problem.constraints) {
      outStream << "object{ Tetrahedron_by_Corners("
        << mapIndexToChar(std::get<0>(chiralityConstraint)) << ", "
        << mapIndexToChar(std::get<1>(chiralityConstraint)) << ", "
        << mapIndexToChar(std::get<2>(chiralityConstraint)) << ", "
        << mapIndexToChar(std::get<3>(chiralityConstraint)) 
        << ", Rl, Rp, 0) pigment{ color steelBlue }}\n";
    }
    outStream << "\n";
  }

  // Gradients
  for(unsigned i = 0; i < N; i++) {
    outStream << "object{ Vector(" << mapIndexToChar(i) <<", " 
      << mapIndexToChar(i) << " + <"
      << std::fixed << std::setprecision(4)
      << stepResult.negativeGradient[3 * i] << ", "
      << stepResult.negativeGradient[3 * i + 1] << ", "
      << stepResult.negativeGradient[3 * i + 2] 
      << ">, 0.05) pigment{ color tomato}}\n";
  }

  outStream.close();
}
