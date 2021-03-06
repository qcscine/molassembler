/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"

#include "Molassembler/IO.h"
#include "Molassembler/IO/Base64.h"
#include "Molassembler/Interpret.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Options.h"
#include "Molassembler/Serialization.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"

using namespace Scine;
using namespace Molassembler;

BOOST_AUTO_TEST_CASE(Base64Reversibility, *boost::unit_test::label("Molassembler")) {
  // Fuzz the encode/decode pair

  const unsigned N = 100;
  for(unsigned i = 0; i < N; ++i) {
    unsigned messageLength = Temple::Random::getSingle<unsigned>(90, 110, randomnessEngine());

    auto sample = Temple::Random::getN<std::uint8_t>(
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
        << ": {" << Temple::condense(sample) << "}.\n"
        << "Encoded : {" << encoded << "}\n"
        << "Decoded : {" << Temple::condense(decoded) << "}"
    );
  }
}

BOOST_AUTO_TEST_CASE(MoleculeSerializationReversibility, *boost::unit_test::label("Molassembler")) {
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ranking_tree_molecules")
  ) {
    auto molecule = IO::read(currentFilePath.string());

    std::string json = JsonSerialization(molecule);
    Molecule decoded = JsonSerialization(json);

    BOOST_CHECK_MESSAGE(
      decoded == molecule,
      "JSON serialization / deserialization failed! The molecule " << currentFilePath.string() << "before and after a serialization/deserialization pair is no longer the same!\n"
        << "JSON representation of original molecule: " << json << "\n"
        << "JSON representation of decoded molecule: "
        << JsonSerialization(json).operator std::string()
    );
  }
}

// After canonicalization, serializations of identical molecules must be identical
BOOST_AUTO_TEST_CASE(MoleculeCanonicalSerialization, *boost::unit_test::label("Molassembler")) {
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
      Interpret::MoleculesResult interpretation;
      if(boc.empty()) {
        // Unfortunately, the file type does not include bond order information
        interpretation = Interpret::molecules(ac, Interpret::BondDiscretizationOption::RoundToNearest);
      } else {
        interpretation = Interpret::molecules(ac, boc, Interpret::BondDiscretizationOption::RoundToNearest);
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
    std::string aSerialization = JsonSerialization(a).standardize();
    std::string bSerialization = JsonSerialization(b).standardize();

    // Ensure canonicalized JSON representations are lexicographically equal
    BOOST_CHECK_MESSAGE(
      aSerialization == bSerialization,
      "After canonicalization, JSON representations of " << currentFilePath << " are not identical:\n\n"
      << aSerialization << "\n\n" << bSerialization << "\n\n\n"
    );
  }
}
