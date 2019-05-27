/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE GraphAlgorithmTestModule
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/range/iterator_range_core.hpp"
#include "boost/test/unit_test.hpp"

#include "molassembler/Cycles.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/IO.h"
#include "temple/Stringify.h"
#include "temple/TinySet.h"

#include "boost/graph/graphviz.hpp"
#include "molassembler/Molecule/MolGraphWriter.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"

#include <iostream>

using namespace Scine;
using namespace molassembler;

inline std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

struct LinkTestData {
  using LinksType = std::set<
    std::pair<AtomIndex, AtomIndex>
  >;

  AtomIndex source;

  LinksType expectedLinks;

  LinkTestData() = default;
  LinkTestData(
    AtomIndex passSource,
    LinksType passLinks
  ) : source {passSource},
      expectedLinks {std::move(passLinks)}
  {}

};

const std::map<std::string, LinkTestData> linkTestData {
  {
    "02_adjusted",
    {
      0,
      {
        {1, 5},
        {18, 22}
      }
    }
  },
  {
    "03",
    {
      0,
      {
        {1, 2}
      }
    }
  },
  {
    "04",
    {
      0,
      {
        {1, 2},
        {1, 3},
        {2, 3}
      }
    }
  },
  {
    "06",
    {
      0,
      {
        {1, 2},
        {3, 4}
      }
    }
  },
  {
    "07",
    {
      0,
      {
        {1, 11},
        {1, 23},
        {11, 17},
        {17, 23}
      }
    }
  },
  {
    "09",
    {
      0,
      {
        {1, 2},
        {1, 3},
        {2, 45},
        {3, 45}
      }
    }
  },
  {
    "11",
    {
      0,
      {
        {1, 2},
        {1, 44},
        {2, 43},
        {43, 44}
      }
    }
  },
  {
    "12",
    {
      0,
      {
        {1, 2},
        {1, 3},
        {2, 3}
      }
    }
  },
  {
    "13",
    {
      0,
      {
        {2, 3},
        {2, 4},
        {7, 8},
        {7, 9},
        {8, 9}
      }
    }
  },
  {
    "15",
    {
      0,
      {
        {2, 3},
        {2, 23},
        {4, 5},
        {4, 52},
        {5, 52}
      }
    }
  },
  {
    "16",
    {
      0,
      {
        {3, 4},
        {3, 23},
        {5, 6},
        {6, 23}
      }
    }
  },
  {
    "18",
    {
      0,
      {
        {7, 14},
        {7, 16},
        {7, 20}
      }
    }
  },
  {
    "19",
    {
      0,
      {
        {3, 4},
        {3, 23},
        {5, 6},
        {6, 23}
      }
    }
  },
  {
    "20",
    {
      0,
      {
        {3, 4},
        {4, 21},
        {1, 21},
        {1, 2}
      }
    }
  },
  {
    "21",
    {
      0,
      {
        {1, 3},
        {1, 24},
        {3, 25},
        {24, 25}
      }
    }
  },
  {
    "23",
    {
      0,
      {
        {1, 2},
        {1, 4},
        {2, 3},
        {3, 4}
      }
    }
  },
  {
    "25",
    {
      0,
      {
        {7, 14},
        {7, 16},
        {7, 20}
      }
    }
  },
  {
    "27",
    {
      0,
      {
        {2, 3},
        {2, 5},
        {3, 4},
        {4, 5}
      }
    }
  },
  {
    "29",
    {
      0,
      {
        {1, 2},
        {3, 4}
      }
    }
  },
  {
    "30",
    {
      0,
      {
        {2, 3},
        {2, 5},
        {3, 4},
        {4, 5}
      }
    }
  },
  {
    "31",
    {
      0,
      {
        {5, 6},
        {7, 8}
      }
    }
  },
  {
    "37",
    {
      0,
      {
        {3, 4},
        {3, 23},
        {5, 6},
        {6, 23}
      }
    }
  },
  {
    "38",
    {
      0,
      {
        {2, 5}
      }
    }
  },
  {
    "39",
    {
      0,
      {
        {3, 32}
      }
    }
  },
  {
    "41",
    {
      0,
      {
        {3, 18}
      }
    }
  },
  {
    "42",
    {
      0,
      {
        {7, 8},
        {7, 10},
        {8, 9},
        {9, 10}
      }
    }
  },
  {
    "44",
    {
      0,
      {
        {7, 14},
        {7, 26},
        {14, 20},
        {20, 26}
      }
    }
  },
  {
    "Co(ox)3",
    {
      0,
      {
        {1, 6},
        {7, 12},
        {13, 18}
      }
    }
  }
};

bool testSubstituentLinks(const boost::filesystem::path& filePath) {
  auto mol = IO::read(filePath.string());

  std::cout << mol << nl;

  const auto& relevantData = linkTestData.at(filePath.stem().string());

  auto ligands = GraphAlgorithms::ligandSiteGroups(
    mol.graph().inner(),
    relevantData.source
  );

  auto links = GraphAlgorithms::siteLinks(
    mol.graph().inner(),
    relevantData.source,
    ligands,
    {}
  );

  using LigandLevelSet = std::set<
    std::pair<unsigned, unsigned>
  >;

  LigandLevelSet condensedCalculated;
  for(const auto& linkData : links) {
    condensedCalculated.insert(linkData.indexPair);
  }

  std::map<AtomIndex, unsigned> indexToLigandMap;
  for(unsigned i = 0; i < ligands.size(); ++i) {
    for(const auto& ligandIndex : ligands.at(i)) {
      indexToLigandMap.emplace(
        ligandIndex,
        i
      );
    }
  }

  LigandLevelSet condensedExpected;
  for(const auto& linkData : relevantData.expectedLinks) {
    auto aLigand = indexToLigandMap.at(linkData.first);
    auto bLigand = indexToLigandMap.at(linkData.second);

    condensedExpected.emplace(
      std::min(aLigand, bLigand),
      std::max(aLigand, bLigand)
    );
  }

  if(condensedExpected != condensedCalculated) {
    std::cout << "Links test fails for " << filePath.stem().string() << "." << nl
      << "Expected: "
      << temple::stringify(condensedExpected) << ", got "
      << temple::stringify(condensedCalculated)
      << ". From:\n";

    for(const auto& link : links) {
      std::cout << temple::stringify(link.indexPair) << ": "
        << temple::stringify(link.cycleSequence) << "\n";
    }

    return false;
  }

  return true;
}

BOOST_AUTO_TEST_CASE(siteLinksTests) {
  boost::filesystem::path filesPath("inorganics/multidentate");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    std::cout << "Stem: " << currentFilePath.stem().string() << nl;
    std::cout << "Extension: " << currentFilePath.extension().string() << nl;

    if(
      currentFilePath.extension() == ".mol"
      && linkTestData.count(currentFilePath.stem().string()) > 0
    ) {
      BOOST_CHECK(
        testSubstituentLinks(currentFilePath)
      );
    }
  }
}

using BondListType = std::vector<
  std::array<AtomIndex, 2>
>;

const std::map<std::string, BondListType> hapticTestData {
  {
    "05",
    {
      {0, 24},
      {0, 26},
      {0, 27},
      {0, 28},
      {0, 29},
      {0, 31}
    }
  },
  {
    "10",
    {
      {0, 7},
      {0, 8}
    }
  },
  {
    "14",
    {
      {0, 7},
      {0, 14}
    }
  },
  {
    "32",
    {
      {0, 4},
      {0, 5},
      {0, 6},
      {0, 7}
    }
  },
};

using LigandsList = std::vector<
  std::vector<AtomIndex>
>;

// Pre-sorted ligands in ascending size and index
const std::map<std::string, LigandsList> hapticLigandsData {
  {
    "05",
    {
       {1},
       {2},
       {24, 26, 27, 28, 29, 31}
    }
  },
  {
    "10",
    {
       {15},
       {17},
       {19},
       {7, 8}
    }
  },
  {
    "14",
    {
       {1},
       {2},
       {7, 14}
    }
  },
  {
    "32",
    {
       {11},
       {12},
       {13},
       {4, 5, 6, 7}
    }
  },
};

bool testHapticBonds(const boost::filesystem::path& filePath) {
  auto rawData = Utils::ChemicalFileHandler::read(filePath.string());

  const unsigned N = rawData.first.size();
  InnerGraph graph (N);

  // Basically a UFF bond discretization scheme
  for(unsigned i = 0; i < N; ++i) {
    for(unsigned j = i + 1; j < N; ++j) {
      double bondOrder = rawData.second.getOrder(i, j);

      if(bondOrder > 0.5) {
        auto bond = static_cast<BondType>(
          std::round(bondOrder) - 1
        );

        if(bondOrder > 6.5) {
          bond = BondType::Sextuple;
        }

        graph.addEdge(i, j, bond);
      }
    }

    // copy element type data too
    graph.elementType(i) = rawData.first.getElement(i);
  }

  GraphAlgorithms::updateEtaBonds(graph);

  const auto& relevantData = hapticTestData.at(filePath.stem().string());

  // Test Eta edge classification correctness
  bool pass = true;
  temple::TinyUnorderedSet<InnerGraph::Edge> expectedEtaBonds;
  for(const auto& expectedEtaBond : relevantData) {
    AtomIndex i = expectedEtaBond.front(),
                  j = expectedEtaBond.back();
    auto edgeOption = graph.edgeOption(i, j);

    if(!edgeOption) {
      std::cout << "The expected eta bond " << i << " - " << j
        << " is not even present in the graph!" << std::endl;
      pass = false;
    }

    expectedEtaBonds.insert(edgeOption.value());
  }

  for(
    const InnerGraph::Edge edge :
    boost::make_iterator_range(
      graph.edges()
    )
  ) {
    AtomIndex i = graph.source(edge);
    AtomIndex j = graph.target(edge);
    BondType bty = graph.bondType(edge);
    bool expectEtaBond = expectedEtaBonds.count(edge) > 0;

    if(expectEtaBond && bty != BondType::Eta) {
      std::cout << "The expected eta bond " << i << " - " << j
        << " is not an Eta bond type, instead has underlying enum value "
        << static_cast<std::underlying_type<BondType>::type>(bty)
        << std::endl;
      pass = false;
    }

    if(!expectEtaBond && bty == BondType::Eta) {
      std::cout << "Did not expect an eta bond for " << i << " - " << j
        << ", but is one in the graph!" << std::endl;
      pass = false;
    }
  }

  // Test ligands classification
  auto ligands = GraphAlgorithms::ligandSiteGroups(graph, 0);

  /* Sort the ligands by size and then by individual atom indices so that we can
   * compare lexicograhically
   */
  std::stable_sort(
    std::begin(ligands),
    std::end(ligands),
    [&](const auto& a, const auto& b) -> bool {
      return a.size() < b.size();
    }
  );

  for(auto& ligand : ligands) {
    std::sort(
      std::begin(ligand),
      std::end(ligand),
      std::less<>()
    );
  }

  auto& expectedLigands = hapticLigandsData.at(filePath.stem().string());
  if(ligands != expectedLigands) {
    pass = false;
    std::cout << "Ligands at the central atom do not match expectation.\n"
      << "Expected: " << temple::stringify(expectedLigands) << "\n"
      << "Got: " << temple::stringify(ligands) << "\n\n";
  }

  if(!pass) {
    StereopermutatorList emptyStereopermutatorList;
    MolGraphWriter propertyWriter(&graph, &emptyStereopermutatorList);

    std::ofstream outFile(filePath.stem().string() + ".dot");
    boost::write_graphviz(
      outFile,
      graph.bgl(),
      propertyWriter,
      propertyWriter,
      propertyWriter
    );

    outFile.close();
  }

  return pass;
}

BOOST_AUTO_TEST_CASE(hapticGraphsTests) {
  boost::filesystem::path filesPath("inorganics/haptic");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    std::cout << "Stem: '" << currentFilePath.stem().string() << "'" << nl;
    std::cout << "Extension: '" << currentFilePath.extension().string() << "'" << nl;

    if(
      currentFilePath.extension() == ".mol"
      && hapticTestData.count(currentFilePath.stem().string()) > 0
    ) {
      BOOST_CHECK(
        testHapticBonds(currentFilePath)
      );
    }
  }
}
