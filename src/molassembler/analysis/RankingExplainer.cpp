/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"

#include "molassembler/Cycles.h"
#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/IO.h"
#include "molassembler/Log.h"
#include "molassembler/Molecule/RankingTree.h"

#include <random>

using namespace Scine;
using namespace molassembler;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

void writeExpandedTree(
  const std::string& fileName,
  const AtomIndex expandOnIndex
) {
  auto molecule = IO::read(
    "../tests/mol_files/ranking_tree_molecules/"s
    + fileName
  );

  auto expandedTree = RankingTree(
    molecule.graph(),
    molecule.graph().cycles(),
    molecule.stereopermutators(),
    molecule.dumpGraphviz(),
    expandOnIndex
  );

  std::ofstream dotFile(fileName + ".dot");
  dotFile << expandedTree.dumpGraphviz();
  dotFile.close();
}

int main(int argc, char* argv[]) {
#ifdef NDEBUG
  std::cout << "This analysis binary requires a debug build of the library.";
  return 0;
#endif

  using namespace std::string_literals;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help,h", "Produce help message")
    (
      "file,f",
      boost::program_options::value<std::string>(),
      "Read molecule to generate from file"
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

  if(options_variables_map.count("help") > 0) {
    std::cout << options_description << nl;
    return 0;
  }

  if(options_variables_map.count("file") > 0) {
    // Set log particulars for debug information
    Log::level = Log::Level::Debug;
    Log::particulars.insert(Log::Particulars::RankingTreeDebugInfo);

    auto filename = options_variables_map["file"].as<std::string>();

    if(!boost::filesystem::exists(filename)) {
      std::cout << "The specified file could not be found!" << nl;
      return 1;
    }

    // This triggers all debug messages during tree instantiations and ranking
    auto mol = IO::read(filename);

    std::cout << mol << std::endl;

    auto filepath = boost::filesystem::path(filename);
    std::string folderName = "ranking-"s + filepath.stem().string();

    // Create a new ranking-results folder if it doesn't exist
    if(!boost::filesystem::is_directory(folderName)) {
      boost::filesystem::create_directory(folderName);
    }

    const boost::regex fileFilterRegex {R"(ranking-tree-[0-9]{1,}-[0-9]\.dot)"};

    boost::filesystem::directory_iterator endIter;

    std::map<
      unsigned, // step index
      unsigned // number of graphs - 1
    > numGraphsMap;

    for(
      boost::filesystem::directory_iterator iter {"."};
      iter != endIter;
      ++iter
    ) {
      if(!boost::filesystem::is_regular_file(iter->status())) {
        continue;
      }

      // Move any generated ranking-tree files to the new folder
      boost::smatch what;

      if(!boost::regex_match(iter -> path().filename().string(), what, fileFilterRegex)) {
        continue;
      }

      boost::filesystem::path newPath(folderName);
      newPath /= iter->path().filename();

      auto splat = StdlibTypeAlgorithms::split(iter->path().filename().string(), '-');
      auto step = std::stoul(splat.at(2));
      auto graphIndex = std::stoul(splat.at(3));

      // Collect information on how many graphs were generated for each step
      if(numGraphsMap.count(step) == 0) {
        numGraphsMap.emplace(
          step,
          graphIndex
        );
      } else if(numGraphsMap.at(step) < graphIndex) {
        numGraphsMap.at(step) = graphIndex;
      }

      // Move the file
      boost::filesystem::rename(
        iter->path(),
        newPath
      );

    }

    // Write a bash file for the generation of combined graphs
    std::ofstream bashFile(folderName + "/create_graphs.sh"s);

    /* Explanation of combined bash command
     * 1. Concatenate streams of separately layouted graphviz graphs, in order
     * 2. Pipe that into gvpack, which combines the graphs into an ordered
     *    array (by which is read in first), into columns with only one row
     *    (redirect stderr to /dev/null, since the warning of node names being
     *    adapted is unneeded)
     * 3. Create the final layout of the combined graph using neato, preserving
     *    position attributes from gvpack (-n2) into an SVG file
     * 4. Write to a iteration-compatibly delimited (i.e. step 0 is 000)
     *    filename
     */

    for(const auto& iterPair : numGraphsMap) {
      const auto& step = iterPair.first;
      const auto& highestGraphIndex = iterPair.second;

      // The first graph is the mol-graph, it needs to be laid out with neato
      bashFile << "cat <(neato ranking-tree-" << step << "-0.dot) ";

      for(unsigned i = 1; i <= highestGraphIndex; ++i) {
        bashFile << "<(dot ranking-tree-" << step << "-" << i << ".dot) ";
      }

      bashFile << "| gvpack -array_uc1 -Glayout=neato 2>/dev/null "
        << "| neato -n2 -Tsvg "
        << "> ranking-tree-"
        << std::setw(3) << std::setfill('0') << step
        << ".svg" << nl;

    }

    bashFile.close();
  }

  return 0;
}
