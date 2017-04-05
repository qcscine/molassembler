#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"

#include "BoundsFromSymmetry.h"
#include "IO.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

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
  const Molecule& molecule,
  const bool& optimized,
  const std::string& symmetryString,
  const unsigned& structNum,
  const Eigen::VectorXd& vectorizedPositions
) {
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

  IO::MOLFileHandler fileHandler;
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
  std::ofstream outStream(symmetryString+".csv"s);

  for(const auto& value: errorValues) {
    outStream << value << std::endl;
  }

  outStream.close();
}


/* NOTES to current state of optimized output structures
 *
 * - Initial problem is gone, was an issue with the generation of the test
 *   matrices
 * - Found one bug in the refinement problem (a missing square), heavily
 *   decreases optimization time
 * - Found another bug that caused generated structures to be considerably
 *   expanded, although well reflective of the overall symmetries (a missing
 *   sqrt in the embedding procedure).
 * - Now high-symmetry geometries are maybe somewhat over-constrained, leading 
 *   to messy structures. Avenues to try:
 *
 *   - Add more variance to 1-3 distance bounds
 *     -> did nothing
 *
 *   - Add more sanity tests to the various symmetries' angle functions
 *     -> Found mistakes in the angle functions of pent-bipy and sq-antiprism, 
 *        fixing the issue
 *        
 * - Refinement cannot reach a global minimum in almost all cases. Probably due
 *   to triangle inequality bound violations. Fix by introducing metrization.
 */

int main() {

  // Make a problem solver
  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGSolver;

  // Set stop criteria
  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  //stopCriteria.iterations = 1000;
  stopCriteria.fDelta = 1e-5;
  DGSolver.setStopCriteria(stopCriteria);

  const unsigned nStructures = 100;

  for(const auto& symmetryName : Symmetry::allNames) {
    // Make a space-free string from the name
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Generate distance bounds
    auto simpleMol = DGDBM::symmetricMolecule(symmetryName);
    auto distanceBoundsMatrix = simpleMol.getDistanceBoundsMatrix();

    /*std::cout << "Sample distances matrix for symmetry '" 
      << Symmetry::name(symmetryName) << std::endl
      << distanceBoundsMatrix.generateDistanceMatrix(
        MetrizationOption::off
      ) << std::endl;*/

    std::vector<double> refinedErrorValues;

    for(unsigned structNum = 0; structNum < nStructures; structNum++) {

      // Calculate metric matrix from selected distances
      MetricMatrix metricMatrix(
        distanceBoundsMatrix.generateDistanceMatrix() // with metrization!
      );

      // Embed
      auto embeddedPositions = metricMatrix.embed(EmbeddingOption::threeDimensional);

      // Vectorize
      Eigen::VectorXd vectorizedPositions(
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      );

      // Write the unoptimized result to a file
      writeMOLFile(
        simpleMol,
        false,
        spaceFreeName,
        structNum,
        vectorizedPositions
      );

      // Create the RefinementProblem
      DGRefinementProblem<double> problem(
        std::vector<ChiralityConstraint>({}), // no chirality constraints
        distanceBoundsMatrix
      );

      // Run the minimization
      DGSolver.minimize(problem, vectorizedPositions);

      // Save end value in Problem
      refinedErrorValues.emplace_back(
        problem.value(
          vectorizedPositions
        )
      );

      // Write the result to a file
      writeMOLFile(
        simpleMol,
        true,
        spaceFreeName,
        structNum,
        vectorizedPositions
      );
    }

    // Write refined error values to file
    writeErrorValues(
      spaceFreeName,
      refinedErrorValues
    );
  }
}
