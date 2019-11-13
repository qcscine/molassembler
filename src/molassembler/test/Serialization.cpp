/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "temple/Random.h"
#include "temple/Stringify.h"

#include "molassembler/IO.h"
#include "molassembler/IO/Base64.h"
#include "molassembler/Interpret.h"
#include "molassembler/Molecule.h"
#include "molassembler/Options.h"
#include "molassembler/Serialization.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(Base64Reversibility) {
  // Fuzz the encode/decode pair

  const unsigned N = 100;
  for(unsigned i = 0; i < N; ++i) {
    unsigned messageLength = temple::random::getSingle<unsigned>(90, 110, randomnessEngine());

    auto sample = temple::random::getN<std::uint8_t>(
      std::numeric_limits<std::uint8_t>::min(),
      std::numeric_limits<std::uint8_t>::max(),
      messageLength,
      randomnessEngine()
    );

    std::string encoded = base64::encode(sample);
    auto decoded = base64::decode(encoded);

    BOOST_CHECK_MESSAGE(
      decoded == sample,
      "Encode / decode pair failed for message of length " << messageLength
        << ": {" << temple::condense(sample) << "}.\n"
        << "Encoded : {" << encoded << "}\n"
        << "Decoded : {" << temple::condense(decoded) << "}"
    );
  }
}

BOOST_AUTO_TEST_CASE(MoleculeSerializationReversibility) {
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ranking_tree_molecules")
  ) {
    auto molecule = IO::read(currentFilePath.string());

    std::string json = JSONSerialization(molecule);
    Molecule decoded = JSONSerialization(json);

    BOOST_CHECK_MESSAGE(
      decoded == molecule,
      "JSON serialization / deserialization failed! The molecule " << currentFilePath.string() << "before and after a serialization/deserialization pair is no longer the same!\n"
        << "JSON representation of original molecule: " << json << "\n"
        << "JSON representation of decoded molecule: "
        << JSONSerialization(json).operator std::string()
    );
  }
}

// After canonicalization, serializations of identical molecules must be identical
BOOST_AUTO_TEST_CASE(MoleculeCanonicalSerialization) {
  boost::filesystem::path directoryBase("isomorphisms");

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    auto readData = Utils::ChemicalFileHandler::read(currentFilePath.string());
    auto permutedData = IO::shuffle(readData.first, readData.second);

    auto interpretSingle = [](const Utils::AtomCollection& ac, const Utils::BondOrderCollection& boc) -> Molecule {
      interpret::MoleculesResult interpretation;
      if(boc.empty()) {
        // Unfortunately, the file type does not include bond order information
        interpretation = interpret::molecules(ac, interpret::BondDiscretizationOption::RoundToNearest);
      } else {
        interpretation = interpret::molecules(ac, boc, interpret::BondDiscretizationOption::RoundToNearest);
      }

      if(interpretation.molecules.size() > 1) {
        throw std::runtime_error(
          std::string("File is not a single molecule, but contains ")
            + std::to_string(interpretation.molecules.size())
            + " components."
        );
      }

      return interpretation.molecules.front();
    };

    Molecule a = interpretSingle(readData.first, readData.second);
    Molecule b = interpretSingle(std::get<0>(permutedData), std::get<1>(permutedData));

    // Canonicalize both molecules
    a.canonicalize();
    b.canonicalize();

    // Get standardized JSON strings for both
    std::string aSerialization = JSONSerialization(a).standardize();
    std::string bSerialization = JSONSerialization(b).standardize();

    // Ensure canonicalized JSON representations are lexicographically equal
    BOOST_CHECK_MESSAGE(
      aSerialization == bSerialization,
      "After canonicalization, JSON representations of " << currentFilePath << " are not identical:\n\n"
      << aSerialization << "\n\n" << bSerialization << "\n\n\n"
    );
  }
}
