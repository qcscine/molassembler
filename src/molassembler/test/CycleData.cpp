#define BOOST_TEST_MODULE RingDecomposition
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "temple/Stringify.h"
#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"

#include "IO.h"
#include <iostream>

/* TODO add further tests for free functions
 */

struct ExpectationData {
  std::vector<unsigned> cycleSizes;

  ExpectationData(std::vector<unsigned>&& passCycleSizes) : cycleSizes {std::move(passCycleSizes)} {
    temple::sort(cycleSizes);
  }
};

std::map<std::string, ExpectationData> decompositionData {
  {
    "fenestrane-3",
    {
      {3, 3, 3, 3}
    }
  },
  {
    "cyclooctatriene-fused-w-cyclopropane",
    {
      {3, 8}
    }
  },
  {
    "rotanes-5",
    {
      {3, 3, 3, 3, 5, 3}
    }
  },
  {
    "strained-db-aromatic-multicycles-1",
    {
      {5, 6, 6, 5, 6, 6}
    }
  },
  {
    "quadricyclane",
    {
      {3, 3, 4, 5, 5}
    }
  },
  {
    "aza-spiro-cyclo-3-3-with-cyclohexanyls",
    {
      {6, 6, 3, 3}
    }
  }
};

void readAndDecompose(const boost::filesystem::path& filePath) {
  using namespace molassembler;

  // Read the file
  auto mol = IO::read(filePath.string());

  Cycles cycles {mol.getGraph()};

  std::vector<unsigned> cycleSizes;
  for(const auto cyclePtr : cycles) {
    cycleSizes.emplace_back(Cycles::size(cyclePtr));
  }

  BOOST_CHECK(!cycleSizes.empty());

  auto findIter = decompositionData.find(filePath.stem().string());

  if(findIter != decompositionData.end()) {
    temple::sort(cycleSizes);

    BOOST_CHECK_MESSAGE(
      cycleSizes == findIter->second.cycleSizes,
      "Expected cycle sizes " << temple::condenseIterable(findIter->second.cycleSizes)
        << ", but got " << temple::condenseIterable(cycleSizes) << " for "
        << filePath.stem().string()
    );

    unsigned cycleSizeThreshold = temple::Math::ceil(temple::average(cycleSizes));

    for(
      const auto cyclePtr :
      cycles.iterate(Cycles::predicates::SizeLessThan(cycleSizeThreshold))
    ) {
      BOOST_CHECK_MESSAGE(
        Cycles::size(cyclePtr) < cycleSizeThreshold,
        "Expected edge set to have size less than " << cycleSizeThreshold
          << ", but got one with size " << Cycles::size(cyclePtr) << ": "
          << temple::stringify(Cycles::edgeVertices(cyclePtr))
      );
    }
  }

  std::cout << "'" << filePath.stem().string() << "' -> "
    << temple::condenseIterable(cycleSizes)
    << std::endl;
}

BOOST_AUTO_TEST_CASE(ringDecomposition) {
  boost::filesystem::path filesPath("test_files/strained_organic_molecules");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    BOOST_CHECK_NO_THROW(
      readAndDecompose(currentFilePath)
    );
  }
}
