#define BOOST_TEST_MODULE SerializationTestModule
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "Serialization.h"
#include "IO.h"
#include "Base64.h"
#include "temple/Random.h"

using namespace molassembler;

BOOST_AUTO_TEST_CASE(base64Tests) {
  // Fuzz the encode/decode pair

  const unsigned N = 100;
  for(unsigned i = 0; i < N; ++i) {
    unsigned messageLength = temple::random.getSingle<unsigned>(90, 110);

    auto sample = temple::random.getN<std::uint8_t>(
      std::numeric_limits<std::uint8_t>::min(),
      std::numeric_limits<std::uint8_t>::max(),
      messageLength
    );

    std::string encoded = base64::encode(sample);
    auto decoded = base64::decode(encoded);

    BOOST_CHECK_MESSAGE(
      decoded == sample,
      "Encode / decode pair failed for message of length " << messageLength
        << ": {" << temple::condenseIterable(sample) << "}.\n"
        << "Encoded : {" << encoded << "}\n"
        << "Decoded : {" << temple::condenseIterable(decoded) << "}"
    );
  }
}

const std::string directoryPrefix = "test_files/ranking_tree_molecules/";

BOOST_AUTO_TEST_CASE(moleculeSerializationTests) {
  boost::filesystem::path filesPath(directoryPrefix);
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    Molecule molecule = IO::read(
      currentFilePath.string()
    );

    std::string json = toJSON(molecule);
    Molecule decoded = fromJSON(json);

    BOOST_CHECK_MESSAGE(
      decoded == molecule,
      "JSON serialization / deserialization failed!\n"
        << "JSON representation of original molecule: " << json << "\n"
        << "JSON representation of decoded molecule: " << toJSON(decoded)
    );
  }
}
