/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "boost/test/unit_test.hpp"

#include "molassembler/Conformers.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "molassembler/Options.h"

#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Stringify.h"

using namespace Scine::Molassembler;

BOOST_AUTO_TEST_CASE(ReproducibleConformers) {
  const unsigned seed = 6564;
  auto& prng = randomnessEngine();

  Molecule mol = IO::read("stereocenter_detection_molecules/RSs-halogenated-propane.mol");
  prng.seed(seed);
  const auto prngStatePriorGeneration = randomnessEngine();
  const auto a = generateRandomConformation(mol);
  const auto prngStatePostGeneration = randomnessEngine();

  prng.seed(seed);
  BOOST_CHECK(randomnessEngine() == prngStatePriorGeneration);
  const auto b = generateRandomConformation(mol);
  BOOST_CHECK(randomnessEngine() == prngStatePostGeneration);

  BOOST_REQUIRE_MESSAGE(
    not (a.has_value() xor b.has_value()),
    "Molecule generation yielded different success/failure states for repeated use of same prng seed"
  );

  BOOST_REQUIRE_MESSAGE(
    a.has_value() && b.has_value(),
    "Molecule generation failed to yield two conformers for RSs-halogenated-propane"
  );

  BOOST_CHECK_MESSAGE(
    a.value().isApprox(b.value(), 1e-3),
    "The two sets of positions returned by repeated use of the same prng are not fully identical."
    << " Difference norm is " << (a.value() - b.value()).norm()
  );
}

BOOST_AUTO_TEST_CASE(ReproducibleEnsembles) {
  const unsigned seed = 6564;
  const unsigned ensembleSize = 10;
  auto& prng = randomnessEngine();

  Molecule mol = IO::read("stereocenter_detection_molecules/RSs-halogenated-propane.mol");
  prng.seed(seed);
  const auto prngStatePriorGeneration = randomnessEngine();
  const auto a = generateRandomEnsemble(mol, ensembleSize);
  const auto prngStatePostGeneration = randomnessEngine();

  prng.seed(seed);
  BOOST_CHECK(randomnessEngine() == prngStatePriorGeneration);
  const auto b = generateRandomEnsemble(mol, ensembleSize);
  BOOST_CHECK(randomnessEngine() == prngStatePostGeneration);

  auto writeEnsemble = [](const auto& resultVector) {
    for(const auto& result : resultVector) {
      if(result) {
        std::cout << result.value().transpose() << "\n\n";
      } else {
        std::cout << "Failure: " << result.error().message() << "\n\n";
      }
    }
  };

  /* This is a little different... If ensemble generation is run in parallel,
   * we do not guarantee ordering, I think. We can't really order the matrices,
   * so we need to do full searching to ensure there is a bijective mapping
   * between a and b.
   */
  std::vector<bool> matchedBs(10, false);
  for(unsigned i = 0; i < ensembleSize; ++i) {
    bool foundMatch = false;

    for(unsigned j = 0; j < ensembleSize; ++j) {
      // Skip already matched results in b
      if(matchedBs.at(j)) {
        continue;
      }

      // Same success/failure state
      if(a.at(i).has_value() == b.at(j).has_value()) {
        if(a.at(i).has_value()) {
          // Both are successes
          if(a.at(i).value().isApprox(b.at(j).value(), 1e-3)) {
            // Both have the same positions
            matchedBs.at(j) = true;
            foundMatch = true;
            break;
          }
        } else {
          // Both are failures
          matchedBs.at(j) = true;
          foundMatch = true;
          break;
        }
      }
    }

    if(!foundMatch) {
      std::cout << "Could not find a match for conformer #" << i << "\n";
      writeEnsemble(a);
      writeEnsemble(b);
      break;
    }
  }

  bool matchedAll = Temple::all_of(matchedBs);

  BOOST_CHECK_MESSAGE(
    matchedAll,
    "Not all conformers could be matched between two re-seeded ensemble generations"
  );
}
