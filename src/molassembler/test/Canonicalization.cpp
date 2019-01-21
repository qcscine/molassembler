/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include "temple/Random.h"
#include "temple/Stringify.h"

#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Graph/Canonicalization.h"

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(moleculeCanonicalizationTests) {
  // Get two isomorphic molecules, calculate the canonical labelling, check that the molecules are 1:1 the same
  using namespace std::string_literals;

  boost::filesystem::path directoryBase("isomorphisms");

  const boost::regex isomorphismFileRegex {R"(.+_isomorphism.mol)"};
  const boost::regex removeRegex {R"(_isomorphism.mol)"};

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    boost::smatch what;

    if(!boost::regex_match(currentFilePath.filename().string(), what, isomorphismFileRegex)) {
      continue;
    }

    auto originalFilePath = currentFilePath.parent_path() / (
      boost::regex_replace(currentFilePath.filename().string(), removeRegex, "") + ".mol"
    );

    if(!boost::filesystem::exists(originalFilePath)) {
      std::cout << "There is no matching file " << originalFilePath
        << " to " << currentFilePath << std::endl;
      continue;
    }

    Molecule a = IO::read(originalFilePath.string());
    Molecule b = IO::read(currentFilePath.string());

    BOOST_REQUIRE_MESSAGE(
      a.graph().N() == b.graph().N(),
      "Purportedly isomorphic molecules a and b do not have the same amount of atoms"
    );

    const AtomIndex N = a.graph().N();

    const auto comparisonBitmask = temple::make_bitmask(AtomEnvironmentComponents::ElementTypes)
      | AtomEnvironmentComponents::BondOrders
      | AtomEnvironmentComponents::Symmetries
      | AtomEnvironmentComponents::Stereopermutations;

    auto aHashes = hashes::generate(a.graph().inner(), a.stereopermutators(), comparisonBitmask);
    auto bHashes = hashes::generate(b.graph().inner(), b.stereopermutators(), comparisonBitmask);

    auto aAutomorphism = canonicalAutomorphism(a, aHashes);
    auto bAutomorphism = canonicalAutomorphism(b, bHashes);

    bool allPass = true;

    for(AtomIndex i = 0; i < N; ++i) {
      if(aHashes.at(aAutomorphism.at(i)) != bHashes.at(bAutomorphism.at(i))) {
        std::cout << "Mismatch at old index " << i << "\n";
        allPass = false;
      }
    }

    if(!allPass) {
      std::cout << temple::stringify(aAutomorphism) << "\n" << temple::stringify(bAutomorphism) << "\n";
    }

    BOOST_CHECK_MESSAGE(
      allPass,
      "Could not match up canonical automorphisms for "
        << originalFilePath.string() << " and " << currentFilePath.string()
    );
  }
}
