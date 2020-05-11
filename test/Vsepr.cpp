/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Modeling/ShapeInference.h"
#include "molassembler/Shapes/Data.h"

#include "boost/optional.hpp"

#include <iostream>

using namespace Scine;
using namespace Molassembler;
using namespace ShapeInference;

using TestCaseType = std::tuple<
  std::string, // a Name for the current compound
  Utils::ElementType, // The central atom type
  std::vector<BindingSite>, // a list of ligand types
  int // charge centered on the central atom
>;

// Runs a specific test case type
void testVsepr(
  const Shapes::Shape& expectedShapeName,
  const std::vector<TestCaseType>& testCases
) {
  std::string complexName;
  int charge;
  std::vector<BindingSite> ligands;
  Utils::ElementType centerAtomType;

  for(const auto& testCase : testCases) {

    std::tie(complexName, centerAtomType, ligands, charge) = testCase;

    boost::optional<Shapes::Shape> shape;

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
          << Shapes::name(shape.value())
          << ", expected: "
          << Shapes::name(expectedShapeName)
      );
    }
  }
}

// Helper function to compose ligand situations
auto makeLigand(
  const unsigned L,
  const unsigned X,
  const Utils::ElementType& type,
  const BondType& bty
) {
  return BindingSite {
    L,
    X,
    {type},
    bty
  };
}

// Helper function to compose ligand situations
std::vector<BindingSite> repeat(
  const BindingSite& ligand,
  const unsigned N
) {
  return std::vector<BindingSite> (N, ligand);
}

// Helper function to compose ligand situations
std::vector<BindingSite> merge(
  const std::vector<BindingSite>& a,
  const std::vector<BindingSite>& b
) {
  std::vector<BindingSite> ret = a;
  std::copy(
    b.begin(),
    b.end(),
    std::back_inserter(ret)
  );

  return ret;
}

BOOST_AUTO_TEST_CASE(VseprTests) {
  using Element = Utils::ElementType;

  testVsepr( // AX2E0
    Shapes::Shape::Line,
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

  testVsepr( // AX3E0
    Shapes::Shape::EquilateralTriangle,
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

  testVsepr( // AX2E1 (V-shape / Bent)
    Shapes::Shape::Bent,
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

  testVsepr( // AX4E0
    Shapes::Shape::Tetrahedron,
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

  testVsepr( // AX3E1
    Shapes::Shape::VacantTetrahedron,
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

  testVsepr( // AX2E2
    Shapes::Shape::Bent,
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

  testVsepr(
    Shapes::Shape::TrigonalBipyramid,
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

  testVsepr(
    Shapes::Shape::Seesaw,
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

  testVsepr(
    Shapes::Shape::T,
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

  testVsepr( // AX2E3
    Shapes::Shape::Line,
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

  testVsepr(
    Shapes::Shape::Octahedron,
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

  testVsepr(
    Shapes::Shape::SquarePyramid,
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

  testVsepr(
    Shapes::Shape::Square,
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

  testVsepr(
    Shapes::Shape::PentagonalBipyramid,
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

  testVsepr(
    Shapes::Shape::PentagonalPyramid,
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

  testVsepr(
    Shapes::Shape::Pentagon,
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

  testVsepr(
    Shapes::Shape::SquareAntiprism,
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
  testVsepr(
    Shapes::Shape::Tetrahedron,
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
