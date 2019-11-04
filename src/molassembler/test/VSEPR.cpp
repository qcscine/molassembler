/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "molassembler/Modeling/LocalGeometryModel.h"
#include "shapes/Data.h"

using namespace Scine;
using namespace molassembler;
using namespace LocalGeometry;

using TestCaseType = std::tuple<
  std::string, // a Name for the current compound
  Scine::Utils::ElementType, // The central atom type
  std::vector<BindingSiteInformation>, // a list of ligand types
  int // charge centered on the central atom
>;

// Runs a specific test case type
void testVSEPR(
  const Symmetry::Shape& expectedShapeName,
  const std::vector<TestCaseType>& testCases
) {
  std::string complexName;
  int charge;
  std::vector<BindingSiteInformation> ligands;
  Scine::Utils::ElementType centerAtomType;

  for(const auto& testCase : testCases) {

    std::tie(complexName, centerAtomType, ligands, charge) = testCase;

    boost::optional<Symmetry::Shape> shape;

    try {
      shape = vsepr(
        centerAtomType,
        ligands,
        charge
      );
    } catch(std::exception& e) {
      std::cout << "Caught an exception when trying "
        << complexName << "." << std::endl;
      throw e;
    }

    BOOST_CHECK_MESSAGE(
      shape,
      "VSEPR could not find a shape for this case"
    );

    if(shape) {
      BOOST_CHECK_MESSAGE(
        shape.value() == expectedShapeName,
        complexName << " has been determined as "
          << Symmetry::name(shape.value())
          << ", expected: "
          << Symmetry::name(expectedShapeName)
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
    Symmetry::Shape::Line,
    {
      TestCaseType {
        "CO2",
        Element::C,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          2
        ),
        0
      },
      TestCaseType {
        "BeF2",
        Element::Be,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          2
        ),
        0
      },
      TestCaseType {
        "C≡C-H",
        Element::C,
        {
          makeLigand(0, 3, Element::C, BondType::Triple),
          makeLigand(0, 1, Element::H, BondType::Single)
        },
        0
      }
    }
  );

  testVSEPR( // AX3E0
    Symmetry::Shape::EquilateralTriangle,
    {
      TestCaseType {
        "BCl3",
        Element::B,
        repeat(
          makeLigand(0, 1, Element::Cl, BondType::Single),
          3
        ),
        0
      },
      TestCaseType {
        "H2C=O",
        Element::C,
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
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          3
        ),
        0
      }
    }
  );

  testVSEPR( // AX2E1 (V-shape / Bent)
    Symmetry::Shape::Bent,
    {
      TestCaseType {
        "SO2",
        Element::S,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          2
        ),
        0
      },
      TestCaseType {
        "CH2··",
        Element::C,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          2
        ),
        0
      },
      TestCaseType {
        "NH2^+",
        Element::N,
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
        {
          makeLigand(0, 1, Element::O, BondType::Single),
          makeLigand(0, 2, Element::O, BondType::Double),
        },
        +1
      }
    }
  );

  testVSEPR( // AX4E0
    Symmetry::Shape::Tetrahedron,
    {
      TestCaseType {
        "CH4",
        Element::C,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          4
        ),
        0
      },
      TestCaseType {
        "XeO4",
        Element::Xe,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          4
        ),
        0
      },
      TestCaseType {
        "ClO4-",
        Element::Cl,
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
    Symmetry::Shape::VacantTetrahedron,
    {
      TestCaseType {
        "NH3",
        Element::N,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          3
        ),
        0
      },
      TestCaseType {
        "XeO3",
        Element::Xe,
        repeat(
          makeLigand(0, 2, Element::O, BondType::Double),
          3
        ),
        0
      }
    }
  );

  testVSEPR( // AX2E2
    Symmetry::Shape::Bent,
    {
      TestCaseType {
        "OH2",
        Element::O,
        repeat(
          makeLigand(0, 1, Element::H, BondType::Single),
          2
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::TrigonalBipyramid,
    {
      TestCaseType {
        "PCl5",
        Element::P,
        repeat(
          makeLigand(0, 1, Element::Cl, BondType::Single),
          5
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::Seesaw,
    {
      TestCaseType {
        "SF4",
        Element::S,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          4
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::T,
    {
      TestCaseType {
        "ClF3",
        Element::Cl,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          3
        ),
        0
      }
    }
  );

  testVSEPR( // AX2E3
    Symmetry::Shape::Line,
    {
      TestCaseType {
        "I3-",
        Element::I,
        {
          makeLigand(0, 1, Element::I, BondType::Single),
          makeLigand(1, 0, Element::I, BondType::Single)
        },
        -1
      },
      TestCaseType {
        "XeF2",
        Element::Xe,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          2
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::Octahedron,
    {
      TestCaseType {
        "SF6",
        Element::S,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          6
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::SquarePyramid,
    {
      TestCaseType {
        "BrF5",
        Element::Br,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          5
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::Square,
    {
      TestCaseType {
        "XeF4",
        Element::Xe,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          4
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::PentagonalBipyramid,
    {
      TestCaseType {
        "IF7",
        Element::I,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          7
        ),
        0
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::PentagonalPyramid,
    {
      TestCaseType {
        "XeOF5^-",
        Element::Xe,
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
    Symmetry::Shape::Pentagon,
    {
      TestCaseType {
        "XeF5-",
        Element::Xe,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          5
        ),
        -1
      }
    }
  );

  testVSEPR(
    Symmetry::Shape::SquareAntiprism,
    {
      TestCaseType {
        "XeF8^2-",
        Element::Xe,
        repeat(
          makeLigand(0, 1, Element::F, BondType::Single),
          8
        ),
        0 // charge is not located on Xe!
      }, TestCaseType {
        "IF8^-",
        Element::I,
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
    Symmetry::Shape::Tetrahedron,
    {
      TestCaseType {
        "H3N -> BF3, on N",
        Element::N,
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
