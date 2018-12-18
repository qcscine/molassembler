/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "molassembler/Modeling/LocalGeometryModel.h"
#include "chemical_symmetries/Symmetries.h"

using namespace Scine;
using namespace molassembler;
using namespace LocalGeometry;

using TestCaseType = std::tuple<
  std::string, // a Name for the current compound
  Scine::Utils::ElementType, // The central atom type
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
  Scine::Utils::ElementType centerAtomType;

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
  const Scine::Utils::ElementType& type,
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
  using Element = Scine::Utils::ElementType;

  testVSEPR( // AX2E0
    Symmetry::Name::Linear,
    {
      TestCaseType {
        "CO2",
        Element::C,
        2,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          2
        ),
        0
      },
      TestCaseType {
        "BeF2",
        Element::Be,
        2,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          2
        ),
        0
      },
      TestCaseType {
        "C≡C-H",
        Element::C,
        2,
        {
          makeLigand(0, 3, Element::C, BondType::Triple),
          makeLigand(0, 1, Element::H, BondType::Single)
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
        Element::B,
        3,
        repeat(
          makeLigand(0, 1, Element::Cl, BondType::Single),
          3
        ),
        0
      },
      TestCaseType {
        "H2C=O",
        Element::C,
        3,
        {
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 2, Element::O, BondType::Double),
        },
        0
      },
      TestCaseType {
        "SO3",
        Element::S,
        3,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
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
        Element::S,
        2,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          2
        ),
        0
      },
      TestCaseType {
        "CH2··",
        Element::C,
        2,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          2
        ),
        0
      },
      TestCaseType {
        "NH2^+",
        Element::N,
        2,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
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
        Element::N,
        2,
        {
          makeLigand(0, 1, Element::O, BondType::Single),
          makeLigand(0, 2, Element::O, BondType::Double),
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
        Element::C,
        4,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          4
        ),
        0
      },
      TestCaseType {
        "XeO4",
        Element::Xe,
        4,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          4
        ),
        0
      },
      TestCaseType {
        "ClO4-",
        Element::Cl,
        4,
        {
          makeLigand(0, 2, Element::O, BondType::Double),
          makeLigand(0, 2, Element::O, BondType::Double),
          makeLigand(0, 2, Element::O, BondType::Double),
          makeLigand(0, 1, Element::O, BondType::Single)
        },
        0 // charge is not on central chlorine atom!
      },
      TestCaseType {
        "N≡SF3",
        Element::S,
        4,
        {
          makeLigand(0, 1, Element::F, BondType::Single),
          makeLigand(0, 1, Element::F, BondType::Single),
          makeLigand(0, 1, Element::F, BondType::Single),
          makeLigand(0, 3, Element::N, BondType::Triple)
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
        Element::N,
        3,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          3
        ),
        0
      },
      TestCaseType {
        "XeO3",
        Element::Xe,
        3,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
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
        Element::O,
        2,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
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
        Element::P,
        5,
        repeat(
          makeLigand(0, 1, Element::Cl, BondType::Single),
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
        Element::S,
        4,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::Cl,
        3,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::I,
        2,
        {
          makeLigand(0, 1, Element::I, BondType::Single),
          makeLigand(1, 0, Element::I, BondType::Single)
        },
        -1
      },
      TestCaseType {
        "XeF2",
        Element::Xe,
        2,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::S,
        6,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::Br,
        5,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::Xe,
        4,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::I,
        7,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::Xe,
        6,
        merge(
          repeat(
            makeLigand(0, 1, Element::F, BondType::Single),
            5
          ),
          repeat(
            makeLigand(0, 2, Element::O, BondType::Double),
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
        Element::Xe,
        5,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::Xe,
        8,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          8
        ),
        0 // charge is not located on Xe!
      }, TestCaseType {
        "IF8^-",
        Element::I,
        8,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
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
        Element::N,
        4,
        {
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 0, Element::B, BondType::Single)
        },
        +1
      },
      TestCaseType {
        "H3N -> BF3, on B",
        Element::B,
        4,
        {
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(0, 1, Element::H, BondType::Single),
          makeLigand(1, 0, Element::N, BondType::Single)
        },
        -1
      }
    }
  );
}
