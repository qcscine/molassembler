#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include "VSEPR.h"

using namespace LocalGeometry;

using TestCaseType = std::tuple<
  std::string, // a Name for the current compound
  Delib::ElementType, // The central atom type
  unsigned, // The number of bonding sites
  std::vector<Model::LigandType>, // a list of ligand types
  int // charge centered on the central atom
>;

// Runs a specific test case type
template<typename ModelClass>
void testModel(
  const Symmetry::Name& expectedSymmetryName,
  const std::vector<TestCaseType>& testCases
) {
  std::string complexName;
  unsigned nSites;
  int charge;
  std::vector<Model::LigandType> ligands;
  Delib::ElementType centerAtomType;

  for(const auto& testCase : testCases) {

    std::tie(complexName, centerAtomType, nSites, ligands, charge) = testCase;

    Symmetry::Name symmetryName;

    try {
      symmetryName = ModelClass::determineGeometry(
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
      symmetryName == expectedSymmetryName,
      complexName.c_str() << " has been determined as " 
        << Symmetry::name(symmetryName).c_str()
        << ", expected: "
        << Symmetry::name(expectedSymmetryName).c_str()
    );

  }
}

// Helper function to compose ligand situations
auto makeLigand(
  const unsigned& L,
  const unsigned& X,
  const Delib::ElementType& type,
  const BondType& bty
) {
  return Model::LigandType {
    L,
    X, 
    {
      {
        type,
        bty
      }
    }
  };
}

// Helper function to compose ligand situations
std::vector<Model::LigandType> repeat(
  const Model::LigandType& ligand,
  const unsigned& N
) {
  return std::vector<Model::LigandType> (N, ligand);
}

// Helper function to compose ligand situations
std::vector<Model::LigandType> merge(
  const std::vector<Model::LigandType>& a,
  const std::vector<Model::LigandType>& b
) {
  std::vector<Model::LigandType> ret = a;
  std::copy(
    b.begin(),
    b.end(),
    std::back_inserter(ret)
  );

  return ret;
}

BOOST_AUTO_TEST_CASE( VSEPRTests ) {
  testModel<VSEPR>( // AX2E0
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

  testModel<VSEPR>( // AX3E0
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

  testModel<VSEPR>( // AX2E1 (V-shape / Bent)
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

  testModel<VSEPR>( // AX4E0
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

  testModel<VSEPR>( // AX3E1
    Symmetry::Name::TrigonalPyramidal,
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

  testModel<VSEPR>( // AX2E2
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>( // AX2E3
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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

  testModel<VSEPR>(
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
  testModel<VSEPR>(
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
