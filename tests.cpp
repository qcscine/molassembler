#define BOOST_TEST_MODULE SymmetryTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Symmetries.h"
#include <set>
#include <iostream>
#include <numeric>
#include <Eigen/Geometry>
#include "constexprProperties.h"

#include "Properties.h"

using namespace Symmetry;

auto makeCoordinateGetter(const Name& name) {
  return [&](const unsigned& index) -> Eigen::Vector3d {
    return symmetryData().at(name).coordinates.at(index);
  };
}

std::vector<unsigned> rotate(
  const std::vector<unsigned>& toRotate,
  const std::vector<unsigned>& rotationVector
) {
  std::vector<unsigned> rotated (toRotate.size()); 

  for(unsigned i = 0; i < toRotate.size(); i++) {
    rotated[i] = toRotate[
      rotationVector[i]
    ];
  }

  return rotated;
}

BOOST_AUTO_TEST_CASE(symmetryDataConstructedCorrectly) {
  BOOST_TEST_REQUIRE(Symmetry::symmetryData().size() == Symmetry::nSymmetries);

  BOOST_REQUIRE(
    TemplateMagic::all_of(
      TemplateMagic::map(
        Symmetry::allNames,
        [&](const auto& symmetryName) -> bool {
          return Symmetry::symmetryData().count(symmetryName) == 1;
        }
      )
    )
  );
}

BOOST_AUTO_TEST_CASE(angleFuntionsInSequence) {
  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::zipMap(
        Symmetry::data::angleFunctions,
        std::vector<Symmetry::data::AngleFunctionPtr> {{
          &Symmetry::data::Linear::angleFunction,
          &Symmetry::data::Bent::angleFunction,
          &Symmetry::data::TrigonalPlanar::angleFunction, // 3
          &Symmetry::data::TrigonalPyramidal::angleFunction,
          &Symmetry::data::TShaped::angleFunction,
          &Symmetry::data::Tetrahedral::angleFunction, // 4
          &Symmetry::data::SquarePlanar::angleFunction,
          &Symmetry::data::Seesaw::angleFunction,
          &Symmetry::data::SquarePyramidal::angleFunction, // 5
          &Symmetry::data::TrigonalBiPyramidal::angleFunction,
          &Symmetry::data::PentagonalPlanar::angleFunction,
          &Symmetry::data::Octahedral::angleFunction, // 6
          &Symmetry::data::TrigonalPrismatic::angleFunction,
          &Symmetry::data::PentagonalPyramidal::angleFunction,
          &Symmetry::data::PentagonalBiPyramidal::angleFunction, // 7
          &Symmetry::data::SquareAntiPrismatic::angleFunction // 8
        }},
        [](const auto& aPtr, const auto& bPtr) -> bool {
          return aPtr == bPtr;
        }
      )
    )
  );
}

BOOST_AUTO_TEST_CASE( correctRotationVectorSize ) {
  // every rotation vector size must equal size of symmetry
  for(const auto& name : allNames) {
    for(const auto& rotationVector : rotations(name)) {
      BOOST_CHECK(rotationVector.size() == size(name));
    }
  }

}

BOOST_AUTO_TEST_CASE( rotationVectorSanityTests ) {

  // every rotation may have every number 0 -> (size of symmetry - 1) only once
  for(const auto& name : allNames) {
    std::set<unsigned> members;
    for(unsigned i = 0; i < size(name); i++) {
      members.insert(members.end(), i);
    }

    for(const auto& rotationVector : rotations(name)) {
      std::set<unsigned> converted {
        rotationVector.begin(),
        rotationVector.end()
      };

      BOOST_CHECK(converted.size() == size(name)); // no duplicates
      
      BOOST_CHECK(
        std::accumulate(
          rotationVector.begin(),
          rotationVector.end(),
          true,
          [&members](const bool& carry, const unsigned& rotationElement) {
            return carry && members.count(rotationElement) == 1;
          }
        )
      );
    }
  }

  /* every rotation must return to the original after a finite number of
   * applications
   */
  unsigned maxIter = 100;
  for(const auto& name : allNames) {
    std::vector<unsigned> initialConfiguration (
      size(name),
      0
    );

    std::iota(
      initialConfiguration.begin(),
      initialConfiguration.end(),
      0
    );

    for(const auto& rotationVector : rotations(name)) {
      // copy in from initial
      auto configuration = initialConfiguration;

      bool pass = false;
      for(unsigned N = 0; N < maxIter; N++) {
        configuration = rotate(configuration, rotationVector);
        if(configuration == initialConfiguration) {
          pass = true;
          break;
        }
      }

      BOOST_CHECK(pass);
    }
  }

}

BOOST_AUTO_TEST_CASE( angleFunctionInputSymmetry ) {

  // every angle function must be symmetrical on input of valid unsigned indices
  for(const auto& symmetryName: allNames) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(symmetryName) && passesAll; i++) {
      for(unsigned j = i + 1; j < size(symmetryName); j++) {
        if(angleFunction(symmetryName)(i, j) != angleFunction(symmetryName)(j, i)) {
          passesAll = false;
          std::cout << name(symmetryName) 
            << " is not symmetrical w.r.t. input indices: falsified by (" 
            << i << ", " << j <<") -> (" << angleFunction(symmetryName)(i, j) 
            << ", " << angleFunction(symmetryName)(j, i) << ")." << std::endl;
          break;
        }
      }
    }

    BOOST_CHECK(passesAll);
  }

}

BOOST_AUTO_TEST_CASE( angleFunctionZeroForIdenticalInput) {

  // every angle function must return 0 for identical indices
  for(const auto& symmetryName: allNames) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(symmetryName); i++) {
      if(angleFunction(symmetryName)(i, i) != 0) {
        passesAll = false;
        std::cout << name(symmetryName)
          << "'s angle function does not return zero for identical indices ("
          << i << ", " << i << ")." << std::endl;
        break;
      }
    }

    BOOST_CHECK(passesAll);
  }
}

BOOST_AUTO_TEST_CASE(anglesWithinRadiansBounds) {
  for(const auto& symmetryName : allNames) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(symmetryName); i++) {
      for(unsigned j = 0; j < size(symmetryName); j++) {
        if(
          !(
            0 <= angleFunction(symmetryName)(i, j) 
          ) || !(
            angleFunction(symmetryName)(i, j) <= M_PI
          )
        ) {
          passesAll = false;
          std::cout << name(symmetryName)
            << "'s angle function is not within radians bounds for indices ("
            << i << ", " << j << ") -> " << angleFunction(symmetryName)(i, j)
            << std::endl;
          break;
        }
      }
    }

    BOOST_CHECK(passesAll);
  }
}

BOOST_AUTO_TEST_CASE( rightAmountOfCoordinates) {
  // every information must have the right amount of coordinates
  for(const auto& symmetryName: allNames) {
    BOOST_CHECK(
      symmetryData().at(symmetryName).coordinates.size() == 
      symmetryData().at(symmetryName).size
    );
  }
}

BOOST_AUTO_TEST_CASE( allCoordinateVectorsLengthOne) {
  for(const auto& symmetryName: allNames) {
    bool all_pass = true;

    for(const auto& coordinate: symmetryData().at(symmetryName).coordinates) {
      if(coordinate.norm() - 1 > 1e10) {
        all_pass = false;
        break;
      }
    }

    BOOST_CHECK(all_pass);
  }
}

BOOST_AUTO_TEST_CASE( anglesMatchCoordinates) {

  /* The results of the angle functions ought to match the geometries specified
   * by the coordinates
   */

  for(const auto& symmetryName: allNames) {
    auto getCoordinates =  [&](const unsigned& index) -> Eigen::Vector3d {
      return symmetryData().at(symmetryName).coordinates.at(index);
    };

    bool all_pass = true;

    for(unsigned i = 0; i < size(symmetryName); i++) {
      for(unsigned j = i + 1; j < size(symmetryName); j++) {
        auto angleInCoordinates = std::acos(
          getCoordinates(i).dot(
            getCoordinates(j)
          ) / (
            getCoordinates(i).norm() * getCoordinates(j).norm()
          )
        );

        auto angleDifference = angleInCoordinates - angleFunction(symmetryName)(i, j);

        // Tolerate only one degree difference
        if(std::fabs(angleDifference) > 1) {
          all_pass = false;

          std::cout << name(symmetryName) 
            << ": angleFunction != angles from coordinates ("  
            << i << ", " << j << "): " << angleDifference 
            << ", angleFunction = " << angleFunction(symmetryName)(i, j) 
            << ", angle from coordinates = " << angleInCoordinates << std::endl;
        }
      }
    }

    BOOST_CHECK(all_pass);
  }
}

BOOST_AUTO_TEST_CASE( allTetrahedraPositive) {
  /* Checks if sequence that tetrahedra are defined in leads to a positive 
   * volume when calculated via
   *
   *  (1 - 4) dot [ (2 - 4) x (3 - 4) ]
   *
   */
  for(const auto& symmetryName: allNames) {
    auto getCoordinates = [&](const boost::optional<unsigned>& indexOption) -> Eigen::Vector3d {
      if(indexOption) {
        return symmetryData().at(symmetryName).coordinates.at(indexOption.value());
      } else {
        return {0, 0, 0};
      }
    };

    bool all_pass = true;

    for(const auto& tetrahedron: tetrahedra(symmetryName)) {

      double tetrahedronVolume = (
        getCoordinates(tetrahedron[0]) - getCoordinates(tetrahedron[3])
      ).dot(
        (
          getCoordinates(tetrahedron[1]) - getCoordinates(tetrahedron[3])
        ).cross(
          getCoordinates(tetrahedron[2]) - getCoordinates(tetrahedron[3])
        )
      );

      if(tetrahedronVolume < 0) {
        all_pass = false;
        std::cout << name(symmetryName) << ": Tetrahedron {";

        for(unsigned i = 0; i < 4; i++) {
          if(tetrahedron[i]) std::cout << tetrahedron[i].value();
          else std::cout << "C";

          if(i != 3) std::cout << ", ";
        }

        std::cout << "} has negative volume (" << tetrahedronVolume << ")." 
          << std::endl;
      }
    }

    BOOST_CHECK(all_pass);
  }
}

BOOST_AUTO_TEST_CASE( tetrahedraDefinitionIndicesUnique ) {
  for(const auto& symmetryName : allNames) {
    for(const auto& tetrahedron : tetrahedra(symmetryName)) {
      bool containsAnEmptyOption = false;

      for(const auto& edgeOption : tetrahedron) {
        if(!edgeOption) {
          containsAnEmptyOption = true;
          break;
        }
      }

      std::set<unsigned> indices;

      for(const auto& edgeOption : tetrahedron) {
        if(edgeOption) {
          indices.insert(edgeOption.value());
        }
      }

      BOOST_CHECK(indices.size() + static_cast<unsigned>(containsAnEmptyOption) == 4);
    }
  }
}

BOOST_AUTO_TEST_CASE(smallestAngleValue) {
  BOOST_CHECK(0 < smallestAngle && smallestAngle < M_PI);
}

/*BOOST_AUTO_TEST_CASE(writeAllCoordinates) {
  for(const auto& symmetryName : allNames) {
    std::cout << Symmetry::name(symmetryName) << "\n";
    std::cout << "ConstexprSymmetryInfo::CoordinatesType {{\n";
    for(const auto& coordinate : symmetryData().at(symmetryName).coordinates) {
      std::cout << "  {" << coordinate[0]
        << ", " << coordinate[1]
        << ", " << coordinate[2]
        << "}\n";
    }
    std::cout << "}}\n\n";
  }
}*/

template<class SymmetryClassFrom, class SymmetryClassTo>
struct LigandGainTest {
  static constexpr bool initialize() {
    constexpr auto constexprReturn = ligandGainMappings<
      SymmetryClassFrom,
      SymmetryClassTo
    >();



    return true;
  }

  static constexpr bool value = initialize();
};

template<class SymmetryClass>
struct RotationGenerationTest {
  static bool value() {

    // This is a DynamicSet of SymmetryClass-sized Arrays
    auto constexprRotations = generateAllRotations<SymmetryClass>(
      iota<SymmetryClass::size>()
    );

    // This is a std::set of SymmetryClass-sized std::vectors
    auto dynamicRotations = properties::generateAllRotations(
      SymmetryClass::name,
      detail::iota<unsigned>(SymmetryClass::size)
    );

    auto convertedRotations = constexprRotations.mapToSTL(
      [&](const auto& indexList) -> std::vector<unsigned> {
        return {
          indexList.begin(),
          indexList.end()
        };
      }
    );

    if(convertedRotations.size() != constexprRotations.size()) {
      std::cout << "In symmetry " << SymmetryClass::stringName << ", "
        << "constexpr rotations set reports " << constexprRotations.size()
        << " elements but the STL mapped variant has only " 
        << convertedRotations.size() << " elements!" << std::endl;
    }

    bool pass = true;

    // Size mismatch
    if(convertedRotations.size() != dynamicRotations.size()) {
      pass = false;
    } else {
      pass = (
        TemplateMagic::setDifference(
          convertedRotations,
          dynamicRotations
        ).size() == 0
      );
    }

    if(!pass) {
      std::cout << "Rotation generation differs for "
        << SymmetryClass::stringName 
        << " symmetry: Sizes of generated sets are different. "
        << "constexpr - " << convertedRotations.size() << " != "
        << dynamicRotations.size() << " - dynamic" << std::endl;
      std::cout << " Maximum #rotations: " << maxRotations<SymmetryClass>() 
        << std::endl;

      std::cout << " Converted constexpr:" << std::endl;
      for(const auto& element : convertedRotations) {
        std::cout << " {" << TemplateMagic::condenseIterable(element) 
          << "}\n";
      }

      std::cout << " Dynamic:" << std::endl;
      for(const auto& element : dynamicRotations) {
        std::cout << " {" << TemplateMagic::condenseIterable(element) 
          << "}\n";
      }
    } 

    return pass;
  }
/* Previously, when the interface with which this is used (unpackToFunction) was
 * unable to cope with both value data members and value function members 
 * equally, this was of the form::
 *
 *   static bool initialize() {
 *     ...
 *     return value;
 *   }
 *   static bool value = initialize();
 *
 * Although this seems equivalent, it really isn't, since initialize() is called
 * at static initialization time instead of at first use as when value is 
 * a function. Since, in this case, the value function depends on another static
 * value (generateAllRotations -> symmetryData), the value-initialize variant
 * leads to a static initialization fiasco, where it is unclear whether the
 * value data member or symmetryData is initialized first.
 *
 * Only in the case of static constexpr is the value-initialize variant 
 * semantically equivalent to the data member variant.
 */
};

BOOST_AUTO_TEST_CASE(constexprProperties) {
  // First the smaller cases
  BOOST_CHECK_MESSAGE(
    TemplateMagic::all_of(
      ConstexprMagic::TupleType::map<
        Symmetry::data::allSymmetryDataTypes,
        RotationGenerationTest
      >()
    ),
    "There is a discrepancy between constexpr and dynamic rotation generation"
  );

  constexpr auto mappings = ligandGainMappings<data::SquarePlanar, data::SquarePyramidal>();
  std::cout << "Linear to TShaped mappings: angular = " 
    << mappings.angleDistortion
    << ", chiral = " << mappings.chiralDistortion
    << ", multiplicity = " << mappings.mappings.size()
    << std::endl;

  for(const auto& indexMapping : mappings.mappings) {
    std::cout << TemplateMagic::condenseIterable(indexMapping) << std::endl;
  }
}
