/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "molassembler/IO.h"
#include "molassembler/Interpret.h"
#include "molassembler/Graph.h"
#include "molassembler/Molecule.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

void readTriples(const std::string& file, Scine::Utils::BondOrderCollection& bondOrders) {
  std::ifstream in(file);
  std::string line;

  unsigned numBondOrders = 0;
  unsigned i,j;
  double value;
  if(in.is_open()) {
    while(true) {
      in >> i >> j >> value;
      if(in.eof() || in.fail()) {
        break;
      }

      bondOrders.setOrder(i, j, value);
      ++numBondOrders;
    }
  }

  std::cout << "Read " << numBondOrders << " bond orders from '" << file << "'\n";
}

Eigen::MatrixXd readCSV(std::string file, int rows, int cols) {
  std::ifstream in(file);

  std::string line;

  int row = 0;
  int col = 0;

  Eigen::MatrixXd res(rows, cols);

  if (in.is_open()) {
    while (std::getline(in, line)) {
      char *ptr = const_cast<char *>(line.c_str());
      int len = line.length();

      col = 0;

      char *start = ptr;
      for (int i = 0; i < len; i++) {

        if (ptr[i] == ',') {
          res(row, col++) = atof(start);
          start = ptr + i + 1;
        }
      }
      res(row, col) = atof(start);

      row++;
    }

    in.close();
  }

  return res;
}

int main(int argc, char* argv[]) {
  using namespace Scine;
  using namespace molassembler;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("m", boost::program_options::value<std::string>(), "Molecule data file (.xyz, .mol, etc)")
    ("bo-csv", boost::program_options::value<std::string>(), "Matrix CSV file to read as bond orders")
    ("bo-triples", boost::program_options::value<std::string>(), "Matrix triples list file to read as bond orders")
    ("threshold", boost::program_options::value<double>(), "Manual threshold for bond stereopermutator interpretation")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  if(options_variables_map.count("m") > 0) {
    boost::filesystem::path filepath {
      options_variables_map["m"].as<std::string>()
    };

    if(!boost::filesystem::exists(filepath)) {
      std::cout << "The molecular file does not exist" << std::endl;
      return 1;
    }

    // This can throw in lots of cases
    auto readData = Utils::ChemicalFileHandler::read(filepath.string());
    const auto& atomCollection = readData.first;
    auto& bondOrders = readData.second;

    if(bondOrders.empty()) {
      if(options_variables_map.count("bo-csv") == 0 && options_variables_map.count("bo-triples") == 0) {
        std::cout << "The molecule data file did not contain bond order information and you did not supply a file to read as bond orders\n";
        return 1;
      }

      const int N = atomCollection.size();
      bondOrders.resize(N);

      if(options_variables_map.count("bo-csv") == 1) {
        auto fileMatrix = readCSV(
          options_variables_map["bo-csv"].as<std::string>(),
          N,
          N
        );
        for(Eigen::Index i = 0; i < N; ++i) {
          for(Eigen::Index j = i + 1; j < N; ++j) {
            if(fileMatrix(i, j) > 0) {
              bondOrders.setOrder(i, j, fileMatrix(i, j));
            }
          }
        }
      } else if(options_variables_map.count("bo-triples") == 1) {
        readTriples(
          options_variables_map["bo-triples"].as<std::string>(),
          bondOrders
        );
      }
    }

    double defaultThreshold = 1.4;
    if(options_variables_map.count("threshold") == 1) {
      defaultThreshold = options_variables_map["threshold"].as<double>();
    }

    auto interpretation = interpret::molecules(atomCollection, bondOrders, interpret::BondDiscretizationOption::RoundToNearest, defaultThreshold);
    auto positions = interpret::applyInterpretationMap(interpretation.componentMap, atomCollection);

    if(interpretation.molecules.size() == 1) {
      std::ofstream graphFile("interpreted.dot");
      graphFile << interpretation.molecules.front().dumpGraphviz();
      graphFile.close();

      io::write("interpreted.json", interpretation.molecules.front());
      io::write("interpreted.mol", interpretation.molecules.front(), positions.front().getPositions());
    } else {
      for(unsigned i = 0; i < interpretation.molecules.size(); ++i) {
        const auto& mol = interpretation.molecules[i];
        std::string filebase = "interpreted-" + std::to_string(i);

        std::ofstream graphFile(filebase + ".dot");
        graphFile << mol.dumpGraphviz();
        graphFile.close();

        io::write(filebase + ".json", mol);
        io::write(filebase + ".mol", mol, positions.at(i).getPositions());
      }
    }

    for(const auto& mol : interpretation.molecules) {
      std::cout << "Interpreted molecule N=" << mol.graph().N() << ": "
        << mol << "\n\n";
    }

  } else {
    std::cout << options_description << std::endl;
  }

  return 0;
}
