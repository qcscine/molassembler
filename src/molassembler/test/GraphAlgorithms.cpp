#define BOOST_TEST_MODULE GraphAlgorithmTestModule
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "GraphAlgorithms.h"
#include "IO.h"
#include "temple/Stringify.h"

#include "boost/graph/graphviz.hpp"
#include "MolGraphWriter.h"

#include <iostream>

using namespace molassembler;

inline std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

struct LinkTestData {
  using LinksType = std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >;

  AtomIndexType source;

  LinksType expectedLinks;

  LinkTestData() = default;
  LinkTestData(
    AtomIndexType passSource,
    const LinksType& passLinks
  ) : source {passSource},
      expectedLinks {passLinks}
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
    mol.getGraph(),
    relevantData.source
  );

  auto links = GraphAlgorithms::substituentLinks(
    mol.getGraph(),
    mol.getCycleData(),
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

  std::map<AtomIndexType, unsigned> indexToLigandMap;
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

BOOST_AUTO_TEST_CASE(substituentLinksTests) {
  boost::filesystem::path filesPath("test_files/inorganics/multidentate");
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
  std::array<AtomIndexType, 2>
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
  std::vector<AtomIndexType>
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

bool testHapticBonds(boost::filesystem::path filePath) {
  IO::MOLFileHandler molHandler;
  auto rawData = molHandler.read(filePath.string());

  const unsigned N = rawData.elements.size();
  GraphType graph (N);

  // Basically a UFF bond discretization scheme
  for(unsigned i = 0; i < N; ++i) {
    for(unsigned j = i + 1; j < N; ++j) {
      double bondOrder = rawData.bondOrders.getOrder(i, j);

      if(bondOrder > 0.5) {
        BondType bond = static_cast<BondType>(
          std::round(bondOrder) - 1
        );

        if(bondOrder > 6.5) {
          bond = BondType::Sextuple;
        }

        auto edgeAddPair = boost::add_edge(i, j, graph);
        graph[edgeAddPair.first].bondType = bond;
      }
    }

    // copy element type data too
    graph[i].elementType = rawData.elements.at(i);
  }

  GraphAlgorithms::findAndSetEtaBonds(graph);

  const auto& relevantData = hapticTestData.at(filePath.stem().string());

  // Test Eta edge classification correctness
  bool pass = true;
  temple::TinyUnorderedSet<GraphType::edge_descriptor> expectedEtaBonds;
  for(const auto& expectedEtaBond : relevantData) {
    AtomIndexType i = expectedEtaBond.front(),
                  j = expectedEtaBond.back();
    auto edgePair = boost::edge(i, j, graph);

    if(!edgePair.second) {
      std::cout << "The expected eta bond " << i << " - " << j
        << " is not even present in the graph!" << std::endl;
      pass = false;
    }

    expectedEtaBonds.insert(edgePair.first);
  }

  for(
    const auto edge :
    RangeForTemporary<GraphType::edge_iterator>(
      boost::edges(graph)
    )
  ) {
    AtomIndexType i = boost::source(edge, graph);
    AtomIndexType j = boost::target(edge, graph);
    BondType bty = graph[edge].bondType;
    bool expectEtaBond = expectedEtaBonds.count(edge);

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
    MolGraphWriter propertyWriter(&graph);

    std::ofstream outFile(filePath.stem().string() + ".dot");
    boost::write_graphviz(
      outFile,
      graph,
      propertyWriter,
      propertyWriter,
      propertyWriter
    );

    outFile.close();
  }

  return pass;
}

BOOST_AUTO_TEST_CASE(hapticGraphsTests) {
  boost::filesystem::path filesPath("test_files/inorganics/haptic");
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
