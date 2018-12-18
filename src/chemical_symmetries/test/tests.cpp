// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#define BOOST_TEST_MODULE SymmetryTests

#include <boost/test/unit_test.hpp>
#include <Eigen/Geometry>

#include "temple/Adaptors/Zip.h"
#include "temple/constexpr/Numeric.h"
#include "temple/constexpr/ToSTL.h"
#include "temple/constexpr/TupleTypePairs.h"
#include "temple/Functional.h"
#include "temple/SetAlgorithms.h"
#include "temple/Stringify.h"
#include "temple/Functor.h"

#include "chemical_symmetries/Symmetries.h"
#include "chemical_symmetries/Properties.h"

#include <set>
#include <iostream>
#include <numeric>

using namespace Scine;
using namespace Symmetry;

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
    temple::all_of(
      Symmetry::allNames,
      [&](const auto& symmetryName) -> bool {
        return Symmetry::symmetryData().count(symmetryName) == 1;
      }
    )
  );
}

BOOST_AUTO_TEST_CASE(angleFuntionsInSequence) {
  BOOST_CHECK(
    temple::all_of(
      temple::adaptors::zip(
        Symmetry::data::angleFunctions,
        std::vector<Symmetry::data::AngleFunctionPtr> {{
          &Symmetry::data::Linear::angleFunction,
          &Symmetry::data::Bent::angleFunction,
          &Symmetry::data::TrigonalPlanar::angleFunction, // 3
          &Symmetry::data::CutTetrahedral::angleFunction,
          &Symmetry::data::TShaped::angleFunction,
          &Symmetry::data::Tetrahedral::angleFunction, // 4
          &Symmetry::data::SquarePlanar::angleFunction,
          &Symmetry::data::Seesaw::angleFunction,
          &Symmetry::data::TrigonalPyramidal::angleFunction,
          &Symmetry::data::SquarePyramidal::angleFunction, // 5
          &Symmetry::data::TrigonalBiPyramidal::angleFunction,
          &Symmetry::data::PentagonalPlanar::angleFunction,
          &Symmetry::data::Octahedral::angleFunction, // 6
          &Symmetry::data::TrigonalPrismatic::angleFunction,
          &Symmetry::data::PentagonalPyramidal::angleFunction,
          &Symmetry::data::PentagonalBiPyramidal::angleFunction, // 7
          &Symmetry::data::SquareAntiPrismatic::angleFunction // 8
        }}
      ),
      [](const auto& aPtr, const auto& bPtr) -> bool {
        return aPtr == bPtr;
      }
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
          [&members](const bool carry, const unsigned rotationElement) {
            return carry && (members.count(rotationElement) == 1);
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
    auto getCoordinates =  [&](const unsigned index) -> Eigen::Vector3d {
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
      }

      return {0, 0, 0};
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
          if(tetrahedron[i]) {
            std::cout << tetrahedron[i].value();
          } else {
            std::cout << "C";
          }

          if(i != 3) {
            std::cout << ", ";
          }
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

BOOST_AUTO_TEST_CASE(smallestAngleValueCorrect) {
  const double comparisonSmallestAngle = temple::min(
    temple::map(
      allNames,
      [](const Name& symmetryName) -> double {
        double symmetrySmallestAngle = angleFunction(symmetryName)(0, 1);

        for(unsigned i = 0; i < size(symmetryName); i++) {
          for(unsigned j = i + 1; j < size(symmetryName); j++) {
            double angle = angleFunction(symmetryName)(i, j);
            if(angle < symmetrySmallestAngle) {
              symmetrySmallestAngle = angle;
            }
          }
        }

        return symmetrySmallestAngle;
      }
    )
  );

  BOOST_CHECK(0 < smallestAngle && smallestAngle < M_PI);
  BOOST_CHECK_MESSAGE(
    std::fabs(
      smallestAngle - comparisonSmallestAngle
    ) < 1e-4,
    "The constant smallest angle set by the library is NOT the smallest "
    << "returned angle within the library. Current value of smallestAngle: "
    << smallestAngle
    << ", true smallest angle:" << comparisonSmallestAngle
  );
}

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
/* NOTE: can refactor out doLigandGainTestIfAdjacent with a simple if-constexpr
 * in C++17
 */
template<class SymmetryClassFrom, class SymmetryClassTo>
std::enable_if_t<
  (
    SymmetryClassFrom::size + 1 == SymmetryClassTo::size
    || SymmetryClassFrom::size == SymmetryClassTo::size
  ),
  bool
> doLigandGainTestIfAdjacent() {
  /* Struct:
   * .mappings - dynamic array of fixed-size index mappings
   * .angularDistortion, .chiralDistortion - doubles
   */
  auto constexprMappings = allMappings.at(
    static_cast<unsigned>(SymmetryClassFrom::name),
    static_cast<unsigned>(SymmetryClassTo::name)
  ).value();
  /* Vector of structs:
   * .indexMapping - vector containing the index mapping
   * .totalDistortion, .chiralDistortion - doubles
   */
  auto dynamicMappings = properties::selectBestTransitionMappings(
    properties::symmetryTransitionMappings(
      SymmetryClassFrom::name,
      SymmetryClassTo::name
    )
  );

  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    properties::floatingPointEqualityThreshold
  };

  if(
    comparator.isUnequal(
      dynamicMappings.angularDistortion,
      constexprMappings.angularDistortion
    ) || comparator.isUnequal(
      dynamicMappings.chiralDistortion,
      constexprMappings.chiralDistortion
    )
  ) {
    return false;
  }

  // Do a full set comparison
  auto convertedMappings = temple::map_stl(
    temple::toSTL(constexprMappings.mappings),
    [&](const auto& indexList) -> std::vector<unsigned> {
      return {
        std::begin(indexList),
        std::end(indexList)
      };
    }
  );

  decltype(convertedMappings) dynamicResultSet {
    dynamicMappings.indexMappings.begin(),
    dynamicMappings.indexMappings.end()
  };

  return temple::set_symmetric_difference(
    convertedMappings,
    dynamicResultSet
  ).empty();
}

using IndexAndMappingsPairType = std::pair<
  unsigned,
  Symmetry::constexprProperties::MappingsReturnType
>;

constexpr bool pairEqualityComparator(
  const IndexAndMappingsPairType& a,
  const IndexAndMappingsPairType& b
) {
  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    Symmetry::properties::floatingPointEqualityThreshold
  };

  return (
    comparator.isEqual(a.second.angularDistortion, b.second.angularDistortion)
    && comparator.isEqual(a.second.chiralDistortion, b.second.chiralDistortion)
  );
}

template<class SymmetryClassFrom, class SymmetryClassTo>
std::enable_if_t<
  SymmetryClassFrom::size == SymmetryClassTo::size + 1,
  bool
> doLigandGainTestIfAdjacent() {
  // Ligand loss situation

  /* Constexpr part */
  temple::Array<
    std::pair<
      unsigned,
      Symmetry::constexprProperties::MappingsReturnType
    >,
    SymmetryClassFrom::size
  > constexprMappings;

  for(unsigned i = 0; i < SymmetryClassFrom::size; ++i) {
    constexprMappings.at(i) = std::make_pair(
      i,
      Symmetry::constexprProperties::ligandLossMappings<
        SymmetryClassFrom,
        SymmetryClassTo
      >(i)
    );
  }

  // Group the results
  auto constexprGroups = temple::groupByEquality(
    constexprMappings,
    pairEqualityComparator // C++17 constexpr lambda
  );

  /* Dynamic part */
  std::vector<
    std::pair<
      unsigned,
      Symmetry::properties::SymmetryTransitionGroup
    >
  > dynamicMappings;

  for(unsigned i = 0; i < SymmetryClassFrom::size; ++i) {
    dynamicMappings.emplace_back(
      i,
      selectBestTransitionMappings(
        properties::ligandLossTransitionMappings(
          SymmetryClassFrom::name,
          SymmetryClassTo::name,
          i
        )
      )
    );
  }

  // Analyze all mappings - which indices have "identical" target mappings?
  auto dynamicGroups = temple::groupByEquality(
    dynamicMappings,
    [&](const auto& firstMappingPair, const auto& secondMappingPair) -> bool {
      return (
        temple::floating::isCloseRelative(
          firstMappingPair.second.angularDistortion,
          secondMappingPair.second.angularDistortion,
          Symmetry::properties::floatingPointEqualityThreshold
        ) && temple::floating::isCloseRelative(
          firstMappingPair.second.chiralDistortion,
          secondMappingPair.second.chiralDistortion,
          Symmetry::properties::floatingPointEqualityThreshold
        )
      );
    }
  );

  /* Comparison */
  // Quick check
  if(dynamicGroups.size() != constexprGroups.size()) {
    return false;
  }

  // Compare an unsigned set of sub-group sizes from each
  std::multiset<unsigned> dynamicGroupSizes, constexprGroupSizes;

  for(const auto& dynamicGroup : dynamicGroups) {
    dynamicGroupSizes.insert(dynamicGroup.size());
  }

  for(const auto& constexprGroup : constexprGroups) {
    constexprGroupSizes.insert(constexprGroup.size());
  }

  return dynamicGroupSizes == constexprGroupSizes;
}

// Base case in which source and target symmetries are non-adjacent
template<class SymmetryClassFrom, class SymmetryClassTo>
std::enable_if_t<
  (
    SymmetryClassFrom::size != SymmetryClassTo::size + 1
    && SymmetryClassFrom::size + 1 != SymmetryClassTo::size
    && SymmetryClassFrom::size != SymmetryClassTo::size
  ),
  bool
> doLigandGainTestIfAdjacent() {
  return true;
}

template<class SymmetryClassFrom, class SymmetryClassTo>
struct LigandGainTest {
  static bool value() {
    return doLigandGainTestIfAdjacent<SymmetryClassFrom, SymmetryClassTo>();
  }
};
#endif

template<class SymmetryClass>
struct RotationGenerationTest {
  static bool value() {

    // This is a DynamicSet of SymmetryClass-sized Arrays
    auto constexprRotations = constexprProperties::generateAllRotations<SymmetryClass>(
      constexprProperties::startingIndexSequence<SymmetryClass>()
    );

    // This is a std::set of SymmetryClass-sized std::vectors
    auto dynamicRotations = properties::generateAllRotations(
      SymmetryClass::name,
      temple::iota<unsigned>(SymmetryClass::size)
    );

    auto convertedRotations = temple::map_stl(
      temple::toSTL(constexprRotations),
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
        temple::set_symmetric_difference(
          convertedRotations,
          dynamicRotations
        ).empty()
      );
    }

    if(!pass) {
      std::cout << "Rotation generation differs for "
        << SymmetryClass::stringName
        << " symmetry: Sizes of generated sets are different. "
        << "constexpr - " << convertedRotations.size() << " != "
        << dynamicRotations.size() << " - dynamic" << std::endl;
      std::cout << " Maximum #rotations: " << constexprProperties::maxRotations<SymmetryClass>()
        << std::endl;

      std::cout << " Converted constexpr:" << std::endl;
      for(const auto& element : convertedRotations) {
        std::cout << " {" << temple::condense(element)
          << "}\n";
      }

      std::cout << " Dynamic:" << std::endl;
      for(const auto& element : dynamicRotations) {
        std::cout << " {" << temple::condense(element)
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

std::string getGraphvizNodeName(const Symmetry::Name& symmetryName) {
  auto stringName = Symmetry::name(symmetryName);

  stringName.erase(
    std::remove_if(
      stringName.begin(),
      stringName.end(),
      [](const char& singleChar) -> bool {
        return (
          singleChar == ' '
          || singleChar == '-'
        );
      }
    ),
    stringName.end()
  );
  return stringName;
}

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
template<typename SymmetryClassFrom, typename SymmetryClassTo>
std::enable_if_t<
  SymmetryClassFrom::size + 1 == SymmetryClassTo::size,
  void
> doWriteIfAdjacent() {
  auto constexprMapping = allMappings.at(
    static_cast<unsigned>(SymmetryClassFrom::name),
    static_cast<unsigned>(SymmetryClassTo::name)
  ).value();

  unsigned multiplicity = constexprMapping.mappings.size();

  std::cout << "  " << getGraphvizNodeName(SymmetryClassFrom::name)
    << " -> " << getGraphvizNodeName(SymmetryClassTo::name)
    << " [";

  if(multiplicity <= 3) {
    std::vector<std::string> repeatColor (
      multiplicity,
      "black"
    );

    std::cout << "color=\"" << temple::condense(repeatColor, ":invis:") << "\"";
  } else {
    std::cout << "color=\"" << "black" << "\"";
    std::cout << ", style=\"dashed\"";
  }

  std::cout << ", label=\""
    << temple::Math::round(constexprMapping.angularDistortion, 2);

  if(multiplicity > 3) {
    std::cout << " (" << multiplicity << ")";
  }

  std::cout << "\"];\n";
}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
std::enable_if_t<
  SymmetryClassFrom::size + 1 != SymmetryClassTo::size,
  void
> doWriteIfAdjacent() {}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
struct WriteLigandMapping {
  static bool value() {
    doWriteIfAdjacent<SymmetryClassFrom, SymmetryClassTo>();
    return true;
  }
};
#endif


BOOST_AUTO_TEST_CASE(constexprPropertiesTests) {
  // Full test of rotation algorithm equivalency for all symmetries
  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::map<
        Symmetry::data::allSymmetryDataTypes,
        RotationGenerationTest
      >(),
      temple::Identity {}
    ),
    "There is a discrepancy between constexpr and dynamic rotation generation"
  );

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
  // Write out all mappings
  temple::TupleType::mapAllPairs<
    Symmetry::data::allSymmetryDataTypes,
    WriteLigandMapping
  >();

  // Test transitions generation/evaluation algorithm equivalency for all
  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::mapAllPairs<
        Symmetry::data::allSymmetryDataTypes,
        LigandGainTest
      >(),
      temple::Identity {}
    ),
    "There is a discrepancy between constexpr and dynamic ligand gain mapping"
    << " generation!"
  );
#endif
}

#ifdef USE_CONSTEXPR_NUM_UNLINKED_ASSIGNMENTS

template<typename SymmetryClass>
struct NumUnlinkedTestFunctor {
  static bool value() {
    auto constexprResults = temple::toSTL(
      Symmetry::allNumUnlinkedAssignments.at(
        static_cast<unsigned>(SymmetryClass::name)
      )
    );

    /* Value for 0 is equal to value for 1, so calculate one less.
     */
    for(unsigned i = 0; i < constexprResults.size() - 1; ++i) {
      if(
        constexprResults.at(i)
        != Symmetry::properties::numUnlinkedAssignments(SymmetryClass::name, i + 1)
      ) {
        return false;
      }
    }

    return true;
  }
};

BOOST_AUTO_TEST_CASE(numUnlinkedTests) {
  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::map<
        Symmetry::data::allSymmetryDataTypes,
        NumUnlinkedTestFunctor
      >()
    ),
    "There is a discrepancy between constexpr and dynamic num unlinked "
    "assignments generation!"
  );

  for(const auto& symmetryName : Symmetry::allNames) {
    std::cout << Symmetry::name(symmetryName) << ": {";
    for(unsigned i = 0; i < Symmetry::size(symmetryName); ++i) {
      std::cout << Symmetry::properties::numUnlinkedAssignments(symmetryName, i);
      if(i != Symmetry::size(symmetryName) - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "}\n";
  }
}
#endif

static_assert(
  nSymmetries == std::tuple_size<data::allSymmetryDataTypes>::value,
  "nSymmetries does not equal number of symmetry data class types in "
  "allSymmetryDataTypes"
);

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
BOOST_AUTO_TEST_CASE(mappingsAreAvailable) {
  /* In every case where allMappings has a value, getMapping must also return
   * a some optional
   */
  bool pass = true;
  for(const auto& fromSymmetry : Symmetry::allNames) {
    auto i = static_cast<unsigned>(fromSymmetry);
    for(const auto& toSymmetry : Symmetry::allNames) {
      auto j = static_cast<unsigned>(toSymmetry);
      if(
        i < j
        && allMappings.at(i, j).hasValue()
          != static_cast<bool>(
            Symmetry::getMapping(fromSymmetry, toSymmetry)
          )
      ) {
        pass = false;
        break;
      }
    }
  }

  BOOST_CHECK_MESSAGE(
    pass,
    "Not all constexpr mappings from allMappings are available from getMapping!"
  );
}
#endif
