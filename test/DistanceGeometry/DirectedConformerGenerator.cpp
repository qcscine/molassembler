/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "Molassembler/DirectedConformerGenerator.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/IO.h"

#include "Fixtures.h"

#include "Utils/Typenames.h"

#include "Molassembler/Temple/Invoke.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/constexpr/Math.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std::string_literals;
using namespace Scine;
using namespace Molassembler;

namespace {

void executeTest(
  const std::string& filename,
  const unsigned numConsideredBonds,
  const unsigned idealEnsembleSize
) {
  auto mol = IO::read(filename);
  DirectedConformerGenerator generator(mol);

  BOOST_CHECK_MESSAGE(
    generator.bondList().size() == numConsideredBonds,
    "Bond list yielded by generator does not have expected size. Expected "
    << numConsideredBonds << " for " << filename << ", got "
    << generator.bondList().size() << " instead."
  );

  BOOST_CHECK_MESSAGE(
    generator.idealEnsembleSize() == idealEnsembleSize,
    "Generator ideal ensemble size does not yield expected number of "
    "conformers. Expected " << idealEnsembleSize << " for " << filename
      << ", got " << generator.idealEnsembleSize() << " instead."
  );

  // If there are
  if(idealEnsembleSize == 0) {
    return;
  }

  // Make a stricter configuration. 2000 steps should be enough, even for testosterone
  const auto fitting = BondStereopermutator::FittingMode::Nearest;
  DistanceGeometry::Configuration configuration {};
  configuration.refinementStepLimit = 2000;

  /* Ensure we can make generate all conformers we have hypothesized exist */
  const unsigned maxTries = 5;
  while(generator.decisionListSetSize() != generator.idealEnsembleSize()) {
    auto newDecisionList = generator.generateNewDecisionList();
    bool success = false;
    for(unsigned attempt = 0; attempt < maxTries; ++attempt) {
      const auto positionResult = generator.generateRandomConformation(newDecisionList, configuration, fitting);
      if(positionResult) {
        success = true;
        break;
      }

      std::cout << "Conformer generation failure: "
        << positionResult.error().message() << "\n";
    }

    BOOST_CHECK_MESSAGE(
      success,
      "Could not generate " << filename << " conformer w/ decision list: "
        << Temple::stringify(newDecisionList) << " in " << maxTries << " attempts"
    );
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(DirectedConformerGeneration, *boost::unit_test::label("DG")) {
  std::vector<
    std::tuple<std::string, unsigned, unsigned>
  > testCases {
    {"directed_conformer_generation/butane.mol", 1, 3},
    {"directed_conformer_generation/pentane.mol", 2, 9},
#ifdef NDEBUG
    {"directed_conformer_generation/caffeine.mol", 3, 8},
    {"isomorphisms/testosterone.mol", 1, 3},
#endif
  };

  for(const auto& tup : testCases) {
    Temple::invoke(executeTest, tup);
  }
}

BOOST_FIXTURE_TEST_CASE(DirectedConfGenHomomorphicSwap, LowTemperatureFixture, *boost::unit_test::label("DG")) {
  auto mol = IO::Experimental::parseSmilesSingleMolecule("CCN");
  DirectedConformerGenerator generator {mol};
  BOOST_REQUIRE_MESSAGE(
    !generator.bondList().empty(),
    "No interesting bonds in CCN at low temperature!"
  );

  while(generator.decisionListSetSize() != generator.idealEnsembleSize()) {
    auto newDecisionList = generator.generateNewDecisionList();
    auto conformerResult = generator.generateRandomConformation(newDecisionList);
    if(conformerResult) {
      auto reinterpretedDecisionlist = generator.getDecisionList(conformerResult.value());
      BOOST_CHECK_MESSAGE(
        newDecisionList == reinterpretedDecisionlist,
        "Reinterpreted decision list is not the same as the one used to generate!"
      );

      /* Speculatively swap a pair of hydrogen atoms at one side of the
       * interesting bond
       */
      const BondIndex interestingBond = generator.bondList().front();

      std::vector<AtomIndex> hydrogens;
      for(const AtomIndex v : mol.graph().adjacents(interestingBond.first)) {
        if(mol.graph().elementType(v) == Utils::ElementType::H) {
          hydrogens.push_back(v);
        }
      }

      BOOST_REQUIRE(hydrogens.size() >= 2);
      const AtomIndex i = hydrogens.at(0);
      const AtomIndex j = hydrogens.at(1);

      Utils::PositionCollection& pos = conformerResult.value();
      pos.row(i).swap(pos.row(j));

      reinterpretedDecisionlist = generator.getDecisionList(pos);
      BOOST_CHECK_MESSAGE(
        newDecisionList == reinterpretedDecisionlist,
        "Swap of hydrogen atom positions yields different decision list!"
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(DirConfGenRelabelerBins, *boost::unit_test::label("DG")) {
  using Relabeler = DirectedConformerGenerator::Relabeler;

  std::vector<double> observedDihedrals {{
    Temple::Math::toRadians(1.0),
    Temple::Math::toRadians(2.0),
    Temple::Math::toRadians(2.4),
    Temple::Math::toRadians(10.0),
    Temple::Math::toRadians(40.0),
    Temple::Math::toRadians(-93.0),
    Temple::Math::toRadians(-95.0)
  }};

  BOOST_CHECK_EQUAL(Relabeler::densityBins(observedDihedrals, 10 * M_PI / 180.0).size(), 3);

  BOOST_CHECK_EQUAL(Relabeler::densityBins(std::vector<double>(1, 3.0), 1.0).size(), 1);

  BOOST_CHECK_EQUAL(Relabeler::densityBins(std::vector<double>(2, 3.0), 1.0).size(), 1);

  BOOST_CHECK_EQUAL(
    Relabeler::densityBins(std::vector<double> {{-M_PI/2, M_PI/2}}, M_PI).size(),
    1
  );

  BOOST_CHECK_EQUAL(
    Relabeler::densityBins(std::vector<double> {{-3 * M_PI / 4, 3 * M_PI / 4}}, 2 * M_PI / 3).size(),
    1
  );

  const auto periodicInputDihedrals = std::vector<double> {{
    Temple::Math::toRadians(1.0),
    Temple::Math::toRadians(119.0)
  }};

  const auto expectedPeriodicOutput = Relabeler::Intervals {
    Relabeler::Interval {
      Temple::Math::toRadians(119.0),
      Temple::Math::toRadians(1.0)
    }
  };

  BOOST_CHECK(
    Relabeler::densityBins(
      periodicInputDihedrals,
      Temple::Math::toRadians(5.0),
      3
    ) == expectedPeriodicOutput
  );
}

BOOST_AUTO_TEST_CASE(DirConfGenRelabeler, *boost::unit_test::label("DG")) {
  auto mol = IO::Experimental::parseSmilesSingleMolecule("CCC(=O)O");
  auto generator = DirectedConformerGenerator(mol);

  std::vector<Utils::PositionCollection> conformers;
  generator.enumerateRandom([&](const auto& /* decList */, const auto& conf) {
    conformers.push_back(conf);
  });

  BOOST_REQUIRE_EQUAL(generator.idealEnsembleSize(), 12);
  BOOST_REQUIRE_MESSAGE(
    conformers.size() >= 10,
    "Require at least 10 conformers to properly test relabeler, got " << conformers.size()
  );

  auto relabeler = generator.relabeler();
  for(const auto& pos : conformers) {
    relabeler.add(pos);
  }

  auto bins = relabeler.bins();
  auto binIndices = relabeler.binIndices(bins);
  auto midpoints = relabeler.binMidpointIntegers(binIndices, bins);
}
