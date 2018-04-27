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

  auto links = GraphAlgorithms::substituentLinks(
    mol.getGraph(),
    mol.getCycleData(),
    relevantData.source,
    mol.getAdjacencies(relevantData.source)
  );

  LinkTestData::LinksType condensed;
  for(const auto& linkData : links) {
    condensed.insert(linkData.indexPair);
  }

  if(condensed != relevantData.expectedLinks) {
    std::cout << "Links test fails for " << filePath.stem().string() << "." << nl
      << "Expected: "
      << temple::stringify(relevantData.expectedLinks) << ", got "
      << temple::stringify(condensed)
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
      {0, 29},
      {0, 31},
      {0, 33}
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

bool testHapticBonds(boost::filesystem::path filePath) {
  IO::MOLFileHandler molHandler;
  auto rawData = molHandler.read(filePath.string());

  unsigned N = rawData.atoms.size();
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
    graph[i].elementType = rawData.atoms.getElement(i);
  }

  graph = GraphAlgorithms::findAndSetEtaBonds(std::move(graph));

  const auto& relevantData = hapticTestData.at(filePath.stem().string());

  /* TODO do this differently. Check ALL bonds, and make sure only those
   * expected were turned into eta bonds. If findAndSetEtaBonds just made
   * everything an eta bond this would pass
   */

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

  GraphType::edge_iterator iter, end;
  std::tie(iter, end) = boost::edges(graph);

  for(; iter != end; ++iter) {
    AtomIndexType i = boost::source(*iter, graph);
    AtomIndexType j = boost::target(*iter, graph);
    BondType bty = graph[*iter].bondType;
    bool expectEtaBond = expectedEtaBonds.count(*iter);

    if(expectEtaBond && bty != BondType::Eta) {
      std::cout << "The expected eta bond " << i << " - " << j
        << " is not an Eta bond type, but has underlying enum value "
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
