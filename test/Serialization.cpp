/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
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
#include "nlohmann/json.hpp"

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

// Comparison of graphs with different Molassembler versions.
BOOST_AUTO_TEST_CASE(StringGraphComparison, *boost::unit_test::label("Molassembler")) {
  std::string graphA = "pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==";
  std::string graphB = "pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA=";
  std::string graphC = "pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAgA=;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwECAA==";
  std::string graphD = "pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA=;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==";
  std::string graphDInvers = "pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==;pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA=";
  std::string graphDDuplicate = "pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==";
  std::string graphDDuplicate2 = "pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==;pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA=";
  std::string graphDDuplicate3 = "pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==;pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA=;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA==";
  std::string invalidGraph = "some;stuff;may;break;...";

  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphA, graphA));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphA, graphB));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphA, graphC));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphA, graphD));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphB, graphC));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphB, graphD));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphC, graphC));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphC, graphD));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphD, graphC));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphDInvers, graphDInvers));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphD, graphDInvers));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphDInvers, graphD));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphDInvers, graphB));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphDInvers, graphDDuplicate));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphDDuplicate, graphDInvers));
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphDDuplicate, graphDDuplicate));
  // A;A;B vs A;B;A
  BOOST_CHECK(JsonSerialization::base64EqualMolecules(graphDDuplicate2, graphDDuplicate3));
  BOOST_CHECK(not JsonSerialization::base64EqualMolecules(graphDDuplicate, graphDDuplicate3));
  BOOST_CHECK_THROW(JsonSerialization::base64EqualMolecules(invalidGraph, invalidGraph), nlohmann::detail::parse_error);
}

BOOST_AUTO_TEST_CASE(DecisionListComparison, *boost::unit_test::label("Molassembler")) {
  std::string listA = "(52, 57, 63, 1):(54, 60, 65, 3)";
  std::string listB = "(1, 7, 12, 1):(-165, -159, -154, 1)";
  std::string listC = ";(70, 76, 81, 2)";
  std::string listD = "(1, 7, 12, 1):(-165, -159, -154, 1);(70, 76, 81, 2)";
  std::string listE = "(1, 7, 12, 1):(-165, -159, -154, 1);;(70, 76, 81, 2)";
  std::string listF = "(173, 178, -177, 1)";
  std::string listG = "(174, 179, -176, 1)";
  std::string listH = "(174,179,-177,1)";
  std::string listI = "(177,-179,-174,1)";
  std::string invalidList = "(52, 57, 1):(54, 60, 65, 3)";
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listA, listA));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listB, listB));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listA, listB));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listB, listA));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listC, listC));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listB, listD));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listD, listD));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listE, listD));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listF, listG));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listG, listF));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listH, listI));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listI, listH));
  BOOST_CHECK_THROW(
    JsonSerialization::equalDecisionLists(invalidList, listA), std::runtime_error);

  // Check close intervals.
  std::string listAa = "(49, 54, 59, 1):(59, 64, 69, 3)";
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listAa, listA));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listA, listAa));

  // Check changing symmetry number.
  std::string listAb = "(49, 54, 59, 2):(59, 64, 69, 3)";
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listAb, listAa));

  // Check symmetry shifted angles
  std::string listAc = "(229, 234, 239, 2):(59, 64, 69, 3)";
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listAb, listAc));
  std::string listAd = "(-131, -126, -121, 2):(179, 184, 189, 3)";
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listAb, listAd));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listAc, listAd));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listAa, listAd));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listAa, listAc));

  // Check permutations of the decision list.
  std::string listEScrambledA = ";(1, 7, 12, 1):(-165, -159, -154, 1);(70, 76, 81, 2)";
  std::string listEScrambledB = ";(70, 76, 81, 2);(1, 7, 12, 1):(-165, -159, -154, 1)";
  std::string listEScrambledC = "(70, 76, 81, 2);(1, 7, 12, 1):(-165, -159, -154, 1);";
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listE, listEScrambledA));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listE, listEScrambledB));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listE, listEScrambledC));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listEScrambledA, listEScrambledA));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listEScrambledA, listEScrambledC));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listEScrambledA, listEScrambledB));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listEScrambledB, listEScrambledC));
  BOOST_CHECK(JsonSerialization::equalDecisionLists(listEScrambledC, listEScrambledB));

  // Check duplicated molecules, e.g. A;A vs A;B
  std::string duplicatedE = ";(1, 7, 12, 1):(-165, -159, -154, 1);(1, 7, 12, 1):(-165, -159, -154, 1)";
  BOOST_CHECK(JsonSerialization::equalDecisionLists(duplicatedE, duplicatedE));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(listEScrambledA, duplicatedE));
  BOOST_CHECK(not JsonSerialization::equalDecisionLists(duplicatedE, listEScrambledA));
}
