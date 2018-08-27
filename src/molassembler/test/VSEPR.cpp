#define BOOST_TEST_MODULE VSEPRTestsModule
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "molassembler/Modeling/LocalGeometryModel.h"
#include "chemical_symmetries/Symmetries.h"

using namespace molassembler;
using namespace molassembler::LocalGeometry;

using TestCaseType = std::tuple<
  std::string, // a Name for the current compound
  Delib::ElementType, // The central atom type
  unsigned, // The number of bonding sites
  std::vector<BindingSiteInformation>, // a list of ligand types
  int // charge centered on the central atom
>;

// Runs a specific test case type
void testVSEPR(
  const Symmetry::Name& expectedSymmetryName,
  const std::vector<TestCaseType>& testCases
) {
  std::string complexName;
  unsigned nSites;
  int charge;
  std::vector<BindingSiteInformation> ligands;
  Delib::ElementType centerAtomType;

  for(const auto& testCase : testCases) {

    std::tie(complexName, centerAtomType, nSites, ligands, charge) = testCase;

    boost::optional<Symmetry::Name> symmetryName;

    try {
      symmetryName = vsepr(
        centerAtomType,
        nSites,
        ligands,
        charge
      );
    } catch(std::exception& e) {
      std::cout << "Caught an exception when trying "
        << complexName << "." << std::endl;
      throw e;
    }

    BOOST_CHECK_MESSAGE(
      symmetryName,
      "VSEPR could not find a symmetry for this case"
    );

    if(symmetryName) {
      BOOST_CHECK_MESSAGE(
        symmetryName.value() == expectedSymmetryName,
        complexName << " has been determined as "
          << Symmetry::name(symmetryName.value())
          << ", expected: "
          << Symmetry::name(expectedSymmetryName)
      );
    }
  }
}

// Helper function to compose ligand situations
auto makeLigand(
  const unsigned L,
  const unsigned X,
  const Delib::ElementType& type,
  const BondType& bty
) {
  return BindingSiteInformation {
    L,
    X,
    {type},
    bty
  };
}

// Helper function to compose ligand situations
std::vector<BindingSiteInformation> repeat(
  const BindingSiteInformation& ligand,
  const unsigned N
) {
  return std::vector<BindingSiteInformation> (N, ligand);
}

// Helper function to compose ligand situations
std::vector<BindingSiteInformation> merge(
  const std::vector<BindingSiteInformation>& a,
  const std::vector<BindingSiteInformation>& b
) {
  std::vector<BindingSiteInformation> ret = a;
  std::copy(
    b.begin(),
    b.end(),
    std::back_inserter(ret)
  );

  return ret;
}

BOOST_AUTO_TEST_CASE( VSEPRTests ) {
  testVSEPR( // AX2E0
    Symmetry::Name::Linear,
    {
      TestCaseType {
        "CO2",
        Delib::ElementType::C,
        2,
        repeat(
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          2
        ),
        0
      },
      TestCaseType {
        "BeF2",
        Delib::ElementType::Be,
        2,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          2
        ),
        0
      },
      TestCaseType {
        "C≡C-H",
        Delib::ElementType::C,
        2,
        {
          makeLigand(0, 3, Delib::ElementType::C, BondType::Triple),
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single)
        },
        0
      }
    }
  );

  testVSEPR( // AX3E0
    Symmetry::Name::TrigonalPlanar,
    {
      TestCaseType {
        "BCl3",
        Delib::ElementType::B,
        3,
        repeat(
          makeLigand(0, 1, Delib::ElementType::Cl, BondType::Single),
          3
        ),
        0
      },
      TestCaseType {
        "H2C=O",
        Delib::ElementType::C,
        3,
        {
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
        },
        0
      },
      TestCaseType {
        "SO3",
        Delib::ElementType::S,
        3,
        repeat(
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          3
        ),
        0
      }
    }
  );

  testVSEPR( // AX2E1 (V-shape / Bent)
    Symmetry::Name::Bent,
    {
      TestCaseType {
        "SO2",
        Delib::ElementType::S,
        2,
        repeat(
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          2
        ),
        0
      },
      TestCaseType {
        "CH2··",
        Delib::ElementType::C,
        2,
        repeat(
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          2
        ),
        0
      },
      TestCaseType {
        "NH2^+",
        Delib::ElementType::N,
        2,
        repeat(
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          2
        ),
        +1
      },
      TestCaseType {
        /* Lewis structure is
         *
         *  O = N⁺ – O⁻  < – >  ⁻O – N⁺ = O
         *
         */
        "NO2·",
        Delib::ElementType::N,
        2,
        {
          makeLigand(0, 1, Delib::ElementType::O, BondType::Single),
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
        },
        +1
      }
    }
  );

  testVSEPR( // AX4E0
    Symmetry::Name::Tetrahedral,
    {
      TestCaseType {
        "CH4",
        Delib::ElementType::C,
        4,
        repeat(
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          4
        ),
        0
      },
      TestCaseType {
        "XeO4",
        Delib::ElementType::Xe,
        4,
        repeat(
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          4
        ),
        0
      },
      TestCaseType {
        "ClO4-",
        Delib::ElementType::Cl,
        4,
        {
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          makeLigand(0, 1, Delib::ElementType::O, BondType::Single)
        },
        0 // charge is not on central chlorine atom!
      },
      TestCaseType {
        "SNF3",
        Delib::ElementType::S,
        4,
        {
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          makeLigand(0, 3, Delib::ElementType::N, BondType::Triple)
        },
        0
      }
    }
  );

  testVSEPR( // AX3E1
    Symmetry::Name::CutTetrahedral,
    {
      TestCaseType {
        "NH3",
        Delib::ElementType::N,
        3,
        repeat(
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          3
        ),
        0
      },
      TestCaseType {
        "XeO3",
        Delib::ElementType::Xe,
        3,
        repeat(
          makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
          3
        ),
        0
      }
    }
  );

  testVSEPR( // AX2E2
    Symmetry::Name::Bent,
    {
      TestCaseType {
        "OH2",
        Delib::ElementType::O,
        2,
        repeat(
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          2
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::TrigonalBiPyramidal,
    {
      TestCaseType {
        "PCl5",
        Delib::ElementType::P,
        5,
        repeat(
          makeLigand(0, 1, Delib::ElementType::Cl, BondType::Single),
          5
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::Seesaw,
    {
      TestCaseType {
        "SF4",
        Delib::ElementType::S,
        4,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          4
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::TShaped,
    {
      TestCaseType {
        "ClF3",
        Delib::ElementType::Cl,
        3,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          3
        ),
        0
      }
    }
  );

  testVSEPR( // AX2E3
    Symmetry::Name::Linear,
    {
      TestCaseType {
        "I3-",
        Delib::ElementType::I,
        2,
        {
          makeLigand(0, 1, Delib::ElementType::I, BondType::Single),
          makeLigand(1, 0, Delib::ElementType::I, BondType::Single)
        },
        -1
      },
      TestCaseType {
        "XeF2",
        Delib::ElementType::Xe,
        2,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          2
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::Octahedral,
    {
      TestCaseType {
        "SF6",
        Delib::ElementType::S,
        6,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          6
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::SquarePyramidal,
    {
      TestCaseType {
        "BrF5",
        Delib::ElementType::Br,
        5,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          5
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::SquarePlanar,
    {
      TestCaseType {
        "XeF4",
        Delib::ElementType::Xe,
        4,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          4
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::PentagonalBiPyramidal,
    {
      TestCaseType {
        "IF7",
        Delib::ElementType::I,
        7,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          7
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Name::PentagonalPyramidal,
    {
      TestCaseType {
        "XeOF5^-",
        Delib::ElementType::Xe,
        6,
        merge(
          repeat(
            makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
            5
          ),
          repeat(
            makeLigand(0, 2, Delib::ElementType::O, BondType::Double),
            1
          )
        ),
        -1
      }
    }
  );

  testVSEPR(
    Symmetry::Name::PentagonalPlanar,
    {
      TestCaseType {
        "XeF5-",
        Delib::ElementType::Xe,
        5,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          5
        ),
        -1
      }
    }
  );

  testVSEPR(
    Symmetry::Name::SquareAntiPrismatic,
    {
      TestCaseType {
        "XeF8^2-",
        Delib::ElementType::Xe,
        8,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          8
        ),
        0 // charge is not located on Xe!
      }, TestCaseType {
        "IF8^-",
        Delib::ElementType::I,
        8,
        repeat(
          makeLigand(0, 1, Delib::ElementType::F, BondType::Single),
          8
        ),
        -1
      }
    }
  );

  // dative bonding
  testVSEPR(
    Symmetry::Name::Tetrahedral,
    {
      TestCaseType {
        "H3N -> BF3, on N",
        Delib::ElementType::N,
        4,
        {
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 0, Delib::ElementType::B, BondType::Single)
        },
        +1
      },
      TestCaseType {
        "H3N -> BF3, on B",
        Delib::ElementType::B,
        4,
        {
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(0, 1, Delib::ElementType::H, BondType::Single),
          makeLigand(1, 0, Delib::ElementType::N, BondType::Single)
        },
        -1
      }
    }
  );
}
