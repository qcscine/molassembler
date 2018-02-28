#define BOOST_TEST_MODULE GraphAlgorithmTestModule
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "GraphAlgorithms.h"
#include "IO.h"
#include "temple/Stringify.h"

using namespace MoleculeManip;

inline std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

struct TestingData {
  using LinksType = RankingInformation::LinksType;

  AtomIndexType source;

  LinksType expectedLinks;

  TestingData() = default;
  TestingData(
    AtomIndexType passSource,
    const LinksType& passLinks
  ) : source {passSource},
      expectedLinks {passLinks}
  {}

};

const std::map<
  std::string,
  TestingData
> testData {
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
  IO::MOLFileHandler molHandler;

  auto mol = molHandler.readSingle(filePath.string());

  std::cout << mol << nl;

  const auto& relevantData = testData.at(filePath.stem().string());

  auto rankedData = mol.rankPriority(relevantData.source);

  if(rankedData.linkedPairs != relevantData.expectedLinks) {
    std::cout << "Links test fails for " << filePath.stem().string() << "." << nl
      << "Expected: "
      << temple::stringify(relevantData.expectedLinks) << ", got "
      << temple::stringify(rankedData.linkedPairs)
      << ".";

    return false;
  }

  return true;
}

BOOST_AUTO_TEST_CASE(substituentLinks) {
  boost::filesystem::path filesPath("../tests/mol_files/inorganics/multidentate");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    std::cout << "Stem: " << currentFilePath.stem().string() << nl;
    std::cout << "Extension: " << currentFilePath.extension().string() << nl;

    if(
      currentFilePath.extension() == ".mol"
      && testData.count(currentFilePath.stem().string()) > 0
    ) {
      BOOST_CHECK(
        testSubstituentLinks(currentFilePath)
      );
    }
  }
}
