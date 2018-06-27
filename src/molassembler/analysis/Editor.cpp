#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "Conformers.h"
#include "IO.h"

using namespace std::string_literals;
using namespace molassembler;

const std::string useExplanation =
  "Combination examples:\n"
  "- Add an atom to an existing position: i, (e or z) & b\n"
  "- Add a bond between existing atom indices: i, j & b\n"
  "- Change a bond type: i, j & b\n"
  "- Change an element type: i & (e or z)\n"
  "- Change a particular geometry: i & s\n"
  "- Remove an atom: i & r\n"
  "- Remove a bond: i, j & r\n";

BondType interpretBondType(const unsigned index) {
  if(index > static_cast<unsigned>(BondType::Eta)) {
    throw std::logic_error("Index specified is not a valid bond type");
  }

  return static_cast<BondType>(index);
}

Delib::ElementType interpretElementType(const std::string& str) {
  try {
    Delib::ElementType e = Delib::ElementInfo::elementTypeForSymbol(str);
    return e;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    throw;
  }
}

const std::map<std::string, std::string> shortLongOptionMap {
  {"M", "molecule"},
  {"O", "output"},
  {"i", "position"},
  {"j", "auxiliary"},
  {"s", "symmetry"},
  {"e", "element_type"},
  {"z", "atomic_number"},
  {"b", "bond_type"}
};

std::string l(const std::string& shortOptionName) {
  return shortLongOptionMap.at(shortOptionName);
}

int main(int argc, char* argv[]) {
/* Set program options from command-line arguments */
  // Defaults
  bool removeFlag = false;

  // Set up option parsing
  boost::program_options::positional_options_description positional_description;
  positional_description.add("molecule,M", -1);

  boost::program_options::options_description mandatory_description("Required arguments");
  mandatory_description.add_options()
    ("molecule,M", boost::program_options::value<std::string>()->required(), "Molecule file to read in")
    ("output,O", boost::program_options::value<std::string>()->required(), "The name of the output file.")
  ;

  boost::program_options::options_description modification_description("Molecule modification arguments");
  modification_description.add_options()
    ("position,i", boost::program_options::value<unsigned>(), "Specify an index at which a modification should be performed")
    ("auxiliary,j", boost::program_options::value<unsigned>(), "Specify the auxiliary index at which a modification should be performed")
    ("symmetry,s", boost::program_options::value<std::string>(), "Specify a symmetry by its full string name")
    ("element_type,e", boost::program_options::value<std::string>(), "Specify an element type by a string")
    ("atomic_number,z", boost::program_options::value<unsigned>(), "Specify an element type by Z")
    ("bond_type,b", boost::program_options::value<unsigned>(), "Specify a bond type by index in BondType enum (0 - single, 1 - double, 2 - triple, ...)")
    ("remove,r", boost::program_options::bool_switch(&removeFlag), "Specify that the modification is a removal")
  ;

  boost::program_options::options_description allOptions("Generic options");
  allOptions.add_options()("help", "Produce help message");
  allOptions.add(mandatory_description).add(modification_description);

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv).
    options(allOptions).
    positional(positional_description).
    style(
      boost::program_options::command_line_style::unix_style
      | boost::program_options::command_line_style::allow_long_disguise
    ).run(),
    options_variables_map
  );

  // Manage the results
  if(options_variables_map.count("help") > 0) {
    std::cout << allOptions << "\n" << useExplanation;
    return 0;
  }

  // Check for required arguments
  boost::program_options::notify(options_variables_map);

  // Check if a file was specified
  if(options_variables_map.count(l("M")) == 0) {
    std::cout << "No file specified to read in." << std::endl;
    std::cout << allOptions << "\n" << useExplanation;
    return 1;
  }

  // Ensure the input file exists
  auto inputFilename = options_variables_map[l("M")].as<std::string>();
  if(!boost::filesystem::exists(inputFilename)) {
    std::cout << "The file '" << inputFilename << "' does not exist." << std::endl;
    return 1;
  }

  // Check that an output filename was specified
  if(options_variables_map.count(l("O")) == 0) {
    std::cout << "No output filename specified." << std::endl;
    std::cout << allOptions << "\n" << useExplanation;
    return 1;
  }

  std::string outputFilename = options_variables_map[l("O")].as<std::string>();
  boost::filesystem::path outputPath {outputFilename};

  if(
    (std::set<std::string> {".mol", ".xyz", ".masm"}).count(
      outputPath.extension().string()
    ) == 0
  ) {
    std::cout << "Unrecognized output file extension. Can only handle .mol, .xyz or .masm"
      << std::endl;
    return 1;
  }

  // Shorthand to check if an option is set
  auto o = [&](const std::string& shortOptionName) -> bool {
    return options_variables_map.count(l(shortOptionName)) == 1;
  };

  if(
    !(
      o("i") && temple::Math::XOR(
        (o("e") || o("z")) && !o("j") && !o("s") && !removeFlag,
        removeFlag && !(o("e") || o("z")) && !o("b") && !o("s"),
        !(o("e") || o("z")) && !removeFlag && temple::Math::XOR(
          o("j") && o("b") && !o("s"),
          !o("j") && !o("b") && o("s")
        )
      )
    )
  ) {
    std::cout << "Your combination of command-line arguments is not sensible.\n"
      << "Perhaps you are trying to do multiple edits in one step?\n"
      << allOptions << "\n" << useExplanation;

    return 1;
  }

  Molecule mol = IO::read(inputFilename);

  std::cout << "Molecule before edit:\n" << mol << std::endl;

  if(o("i")) {
    AtomIndexType i = options_variables_map[l("i")].as<unsigned>();

    if(o("b") && (o("e") || o("z"))) {
      BondType b = interpretBondType(
        options_variables_map[l("b")].as<unsigned>()
      );

      Delib::ElementType e;
      if(o("e")) {
        e = interpretElementType(
          options_variables_map[l("e")].as<std::string>()
        );
      } else {
        e = static_cast<Delib::ElementType>(
          options_variables_map[l("z")].as<unsigned>()
        );
      }

      mol.addAtom(e, i, b);
    }

    if(!o("b") && (o("e") || o("z"))) {
      Delib::ElementType e;
      if(o("e")) {
        e = interpretElementType(
          options_variables_map[l("e")].as<std::string>()
        );
      } else {
        e = static_cast<Delib::ElementType>(
          options_variables_map[l("z")].as<unsigned>()
        );
      }

      mol.setElementType(i, e);
    }

    if(o("j") && o("b")) {
      AtomIndexType j = options_variables_map[l("j")].as<unsigned>();
      BondType b = interpretBondType(
        options_variables_map[l("b")].as<unsigned>()
      );

      if(mol.isAdjacent(i, j)) {
        std::cout << "The indices " << i << " and " << j
          << " are adjacent. Changing their bond type." << std::endl;
        mol.setBondType(i, j, b);
      } else {
        std::cout << "The indices " << i << " and " << j
          << " are not adjacent. Adding a bond." << std::endl;
        mol.addBond(i, j, b);
      }
    }

    if(removeFlag && !o("j")) {
      mol.removeAtom(i);
    }

    if(removeFlag && o("j")) {
      AtomIndexType j = options_variables_map[l("j")].as<unsigned>();
      mol.removeBond(i, j);
    }

    if(o("s")) {
      // Try to interpret s as a symmetry name
      Symmetry::Name newGeometry;

      try {
        newGeometry = Symmetry::nameFromString(
          options_variables_map[l("s")].as<std::string>()
        );
      } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        std::cout << "Existing geometry names:\n";
        for(const auto name : Symmetry::allNames) {
          std::cout << "- " << Symmetry::name(name) << "\n";
        }

        return 1;
      }

      mol.setGeometryAtAtom(i, newGeometry);
    }
  }

  std::cout << "Molecule after edit: " << mol << std::endl;

  if(outputPath.extension() == ".mol" || outputPath.extension() == ".xyz") {
    // Generate a conformation for the new molecule
    auto conformerGenerationResult = generateConformation(mol);
    if(conformerGenerationResult) {
      IO::write(
        outputFilename,
        mol,
        conformerGenerationResult.value()
      );
    } else {
      std::cout << "No conformation could be generated for the resulting molecule." << std::endl;
      return 1;
    }
  } else if(outputPath.extension() == ".masm") {
    IO::write(outputFilename, mol);
  }

  return 0;
}
