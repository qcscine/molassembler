/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/range/iterator_range_core.hpp"
#include "boost/test/unit_test.hpp"

#include "temple/Stringify.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"

#include "molassembler/IO.h"
#include "molassembler/Cycles.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/Molecule.h"
#include <iostream>
#include <array>
#include <iterator>

/* TODO add further tests for free functions
 */

using namespace Scine;
using namespace molassembler;


struct ExpectationData {
  std::vector<unsigned> cycleSizes;

  ExpectationData(std::vector<unsigned>&& passCycleSizes)
    : cycleSizes(temple::sort(passCycleSizes)) {}
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
  auto findIter = decompositionData.find(filePath.stem().string());

  if(findIter != decompositionData.end()) {
    // Read the file
    auto mol = io::read(filePath.string());

    const Cycles& cycles = mol.graph().cycles();

    auto cycleSizes = temple::map(
      cycles,
      [](const auto& cycleEdges) -> unsigned {
        return cycleEdges.size();
      }
    );

    BOOST_CHECK(!cycleSizes.empty());

    temple::inplace::sort(cycleSizes);

    BOOST_CHECK_MESSAGE(
      cycleSizes == findIter->second.cycleSizes,
      "Expected cycle sizes " << temple::condense(findIter->second.cycleSizes)
        << ", but got " << temple::condense(cycleSizes) << " for "
        << filePath.stem().string()
    );
  }
}

BOOST_AUTO_TEST_CASE(ringDecomposition) {
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("strained_organic_molecules")
  ) {
    BOOST_CHECK_NO_THROW(
      readAndDecompose(currentFilePath)
    );
  }
}

BOOST_AUTO_TEST_CASE(cycleIterators) {
  std::vector<
    std::tuple<std::string, AtomIndex, unsigned>
  > tests {
    {"strained_organic_molecules/fenestrane-4.mol", 0ul, 4u}
  };

  std::string file;
  AtomIndex i;
  unsigned relevantCycleCount;
  for(const auto& testTuple : tests) {
    std::tie(file, i, relevantCycleCount) = testTuple;

    Molecule mol;

    BOOST_REQUIRE_NO_THROW(mol = io::read(file));

    const Cycles& cycles = mol.graph().cycles();

    auto iteratorPair = cycles.containing(i);

    BOOST_REQUIRE(iteratorPair.first != iteratorPair.second);

    const unsigned iteratorDistance = std::distance(iteratorPair.first, iteratorPair.second);

    BOOST_CHECK_MESSAGE(
      iteratorDistance == relevantCycleCount,
      "Iterator distance for " << file << " is " << iteratorDistance << ", not " << relevantCycleCount
    );
  }
}
