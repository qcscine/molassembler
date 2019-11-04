/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE SymmetryTests

#include <boost/test/unit_test.hpp>
#include <Eigen/Geometry>

#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Zip.h"
#include "temple/Functional.h"
#include "temple/Functor.h"
#include "temple/SetAlgorithms.h"
#include "temple/Stringify.h"
#include "temple/constexpr/Numeric.h"
#include "temple/constexpr/ToSTL.h"
#include "temple/constexpr/TupleTypePairs.h"

#include "shapes/PropertyCaching.h"
#include "shapes/Data.h"

#include <set>
#include <numeric>
#include <iostream>

using namespace Scine;
using namespace Symmetry;

template<typename EnumType>
constexpr inline auto underlying(const EnumType e) {
  return static_cast<std::underlying_type_t<EnumType>>(e);
}

std::vector<unsigned> rotate(
  const std::vector<unsigned>& toRotate,
  const std::vector<unsigned>& rotationVector
) {
  std::vector<unsigned> rotated (toRotate.size());

  for(unsigned i = 0; i < toRotate.size(); i++) {
    assert(rotationVector[i] < toRotate.size());
    rotated[i] = toRotate[
      rotationVector[i]
    ];
  }

  return rotated;
}

template<typename ShapeClass>
struct LockstepTest {
  static bool value() {
    return Symmetry::nameIndex(ShapeClass::shape) == underlying(ShapeClass::shape);
  }
};

BOOST_AUTO_TEST_CASE(SymmetryTypeAndPositionInEnumLockstep) {
  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::map<
        Symmetry::data::allShapeDataTypes,
        LockstepTest
      >()
    ),
    "Not all symmetries have the same order in Name and allShapeDataTypes"
  );
}

BOOST_AUTO_TEST_CASE(AngleFunctionInputSymmetry) {
  // every angle function must be symmetrical on input of valid unsigned indices
  for(const Shape shape: allShapes) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(shape) && passesAll; i++) {
      for(unsigned j = i + 1; j < size(shape); j++) {
        if(angleFunction(shape)(i, j) != angleFunction(shape)(j, i)) {
          passesAll = false;
          std::cout << name(shape)
            << " is not symmetrical w.r.t. input indices: falsified by ("
            << i << ", " << j <<") -> (" << angleFunction(shape)(i, j)
            << ", " << angleFunction(shape)(j, i) << ")." << std::endl;
          break;
        }
      }
    }

    BOOST_CHECK(passesAll);
  }

}

BOOST_AUTO_TEST_CASE(AngleFunctionZeroForIdenticalInput) {
  // every angle function must return 0 for identical indices
  for(const Shape shape: allShapes) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(shape); i++) {
      if(angleFunction(shape)(i, i) != 0) {
        passesAll = false;
        std::cout << name(shape)
          << "'s angle function does not return zero for identical indices ("
          << i << ", " << i << ")." << std::endl;
        break;
      }
    }

    BOOST_CHECK(passesAll);
  }
}

BOOST_AUTO_TEST_CASE(AnglesWithinRadiansBounds) {
  for(const Shape shape : allShapes) {
    bool passesAll = true;

    for(unsigned i = 0; i < size(shape); i++) {
      for(unsigned j = 0; j < size(shape); j++) {
        if(
          !(
            0 <= angleFunction(shape)(i, j)
          ) || !(
            angleFunction(shape)(i, j) <= M_PI
          )
        ) {
          passesAll = false;
          std::cout << name(shape)
            << "'s angle function is not within radians bounds for indices ("
            << i << ", " << j << ") -> " << angleFunction(shape)(i, j)
            << std::endl;
          break;
        }
      }
    }

    BOOST_CHECK(passesAll);
  }
}

BOOST_AUTO_TEST_CASE(AnglesMatchCoordinates) {

  /* The results of the angle functions ought to match the geometries specified
   * by the coordinates
   */

  for(const auto& shape: allShapes) {
    auto getCoordinates =  [&](const unsigned index) -> Eigen::Vector3d {
      return symmetryData().at(shape).coordinates.col(index);
    };

    bool all_pass = true;

    for(unsigned i = 0; i < size(shape); i++) {
      for(unsigned j = i + 1; j < size(shape); j++) {
        auto angleInCoordinates = std::acos(
          getCoordinates(i).dot(
            getCoordinates(j)
          ) / (
            getCoordinates(i).norm() * getCoordinates(j).norm()
          )
        );

        auto angleDifference = angleInCoordinates - angleFunction(shape)(i, j);

        // Tolerate only one degree difference
        if(std::fabs(angleDifference) > 1) {
          all_pass = false;

          std::cout << name(shape)
            << ": angleFunction != angles from coordinates ("
            << i << ", " << j << "): " << angleDifference
            << ", angleFunction = " << angleFunction(shape)(i, j)
            << ", angle from coordinates = " << angleInCoordinates << std::endl;
        }
      }
    }

    BOOST_CHECK(all_pass);
  }
}

BOOST_AUTO_TEST_CASE(AllTetrahedraPositive) {
  /* Checks if sequence that tetrahedra are defined in leads to a positive
   * volume when calculated via
   *
   *  (1 - 4) dot [ (2 - 4) x (3 - 4) ]
   *
   */
  for(const auto& shape: allShapes) {
    auto getCoordinates = [&](const boost::optional<unsigned>& indexOption) -> Eigen::Vector3d {
      if(indexOption) {
        return symmetryData().at(shape).coordinates.col(indexOption.value());
      }

      return {0, 0, 0};
    };

    bool all_pass = true;

    for(const auto& tetrahedron: tetrahedra(shape)) {

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
        std::cout << name(shape) << ": Tetrahedron {";

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

BOOST_AUTO_TEST_CASE(TetrahedraDefinitionIndicesUnique) {
  for(const auto& shape : allShapes) {
    for(const auto& tetrahedron : tetrahedra(shape)) {
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

BOOST_AUTO_TEST_CASE(SmallestAngleValueCorrect) {
  const double comparisonSmallestAngle = temple::min(
    temple::map(
      allShapes,
      [](const Shape& shape) -> double {
        double symmetrySmallestAngle = angleFunction(shape)(0, 1);

        for(unsigned i = 0; i < size(shape); i++) {
          for(unsigned j = i + 1; j < size(shape); j++) {
            double angle = angleFunction(shape)(i, j);
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
template<class ShapeClassFrom, class ShapeClassTo>
std::enable_if_t<
  (
    ShapeClassFrom::size + 1 == ShapeClassTo::size
    || ShapeClassFrom::size == ShapeClassTo::size
  ),
  bool
> doLigandGainTestIfAdjacent() {
  /* Struct:
   * .mappings - dynamic array of fixed-size index mappings
   * .angularDistortion, .chiralDistortion - doubles
   */
  auto constexprMappings = allMappings.at(
    static_cast<unsigned>(ShapeClassFrom::shape),
    static_cast<unsigned>(ShapeClassTo::shape)
  ).value();
  /* Vector of structs:
   * .indexMapping - vector containing the index mapping
   * .totalDistortion, .chiralDistortion - doubles
   */
  auto dynamicMappings = properties::selectBestTransitionMappings(
    properties::symmetryTransitionMappings(
      ShapeClassFrom::shape,
      ShapeClassTo::shape
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

template<class ShapeClassFrom, class ShapeClassTo>
std::enable_if_t<
  ShapeClassFrom::size == ShapeClassTo::size + 1,
  bool
> doLigandGainTestIfAdjacent() {
  // Ligand loss situation

  /* Constexpr part */
  temple::Array<
    std::pair<
      unsigned,
      Symmetry::constexprProperties::MappingsReturnType
    >,
    ShapeClassFrom::size
  > constexprMappings;

  for(unsigned i = 0; i < ShapeClassFrom::size; ++i) {
    constexprMappings.at(i) = std::make_pair(
      i,
      Symmetry::constexprProperties::ligandLossMappings<
        ShapeClassFrom,
        ShapeClassTo
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

  for(unsigned i = 0; i < ShapeClassFrom::size; ++i) {
    dynamicMappings.emplace_back(
      i,
      selectBestTransitionMappings(
        properties::ligandLossTransitionMappings(
          ShapeClassFrom::shape,
          ShapeClassTo::shape,
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
template<class ShapeClassFrom, class ShapeClassTo>
std::enable_if_t<
  (
    ShapeClassFrom::size != ShapeClassTo::size + 1
    && ShapeClassFrom::size + 1 != ShapeClassTo::size
    && ShapeClassFrom::size != ShapeClassTo::size
  ),
  bool
> doLigandGainTestIfAdjacent() {
  return true;
}

template<class ShapeClassFrom, class ShapeClassTo>
struct LigandGainTest {
  static bool value() {
    return doLigandGainTestIfAdjacent<ShapeClassFrom, ShapeClassTo>();
  }
};
#endif

template<class ShapeClass>
struct RotationGenerationTest {
  static bool value() {

    // This is a DynamicSet of ShapeClass-sized Arrays
    auto constexprRotations = constexprProperties::generateAllRotations<ShapeClass>(
      constexprProperties::startingIndexSequence<ShapeClass>()
    );

    // This is a std::set of ShapeClass-sized std::vectors
    auto dynamicRotations = properties::generateAllRotations(
      ShapeClass::shape,
      temple::iota<unsigned>(ShapeClass::size)
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
      std::cout << "In symmetry " << ShapeClass::stringName << ", "
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
        << ShapeClass::stringName
        << " shape: Sizes of generated sets are different. "
        << "constexpr - " << convertedRotations.size() << " != "
        << dynamicRotations.size() << " - dynamic" << std::endl;
      std::cout << " Maximum #rotations: " << constexprProperties::maxRotations<ShapeClass>()
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

std::string getGraphvizNodeName(const Symmetry::Shape shape) {
  auto stringName = Symmetry::name(shape);

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

BOOST_AUTO_TEST_CASE(constexprPropertiesTests) {
  // Full test of rotation algorithm equivalency for all symmetries
  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::map<
        Symmetry::data::allShapeDataTypes,
        RotationGenerationTest
      >()
    ),
    "There is a discrepancy between constexpr and dynamic rotation generation"
  );

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
  // Test transitions generation/evaluation algorithm equivalency for all
  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::mapAllPairs<
        Symmetry::data::allShapeDataTypes,
        LigandGainTest
      >()
    ),
    "There is a discrepancy between constexpr and dynamic ligand gain mapping"
    << " generation!"
  );
#endif
}

template<typename ShapeClass>
struct NumUnlinkedTestFunctor {
  static bool value() {
    for(unsigned i = 1; i < ShapeClass::size; ++i) {
      unsigned constexprResult = constexprProperties::numUnlinkedStereopermutations<ShapeClass>(i);

      unsigned dynamicResult = properties::numUnlinkedStereopermutations(ShapeClass::shape, i);

      if(constexprResult != dynamicResult) {
        std::cout << "Mismatch for " << Symmetry::name(ShapeClass::shape) << " and " << i << " identical ligands between constexpr and dynamic number of unlinked: " << constexprResult << " vs. " << dynamicResult << "\n";
        return false;
      }

      // Cross-check with constexpr hasMultiple
      bool constexprHasMultiple = constexprProperties::hasMultipleUnlinkedStereopermutations<ShapeClass>(i);
      if((constexprResult > 1) != constexprHasMultiple) {
        std::cout << "Mismatch between constexpr count and constexpr "
          << "hasMultiple unlinked ligands for "
          << Symmetry::name(ShapeClass::shape) << " and "
          << i << " identical ligands: " << constexprResult << " and "
          << std::boolalpha << constexprHasMultiple << "\n";
        return false;
      }

      // Cross-check with dynamic hasMultiple
      bool dynamicHasMultiple = properties::hasMultipleUnlinkedStereopermutations(ShapeClass::shape, i);
      if((constexprResult > 1) != dynamicHasMultiple) {
        std::cout << "Mismatch between constexpr count and dynamic "
          << "hasMultiple unlinked ligands for "
          << Symmetry::name(ShapeClass::shape) << " and "
          << i << " identical ligands: " << constexprResult << " and "
          << std::boolalpha << dynamicHasMultiple << "\n";
        return false;
      }
    }

    return true;
  }
};

BOOST_AUTO_TEST_CASE(numUnlinkedAlgorithms) {
  using TestTypes = std::tuple<
    data::Line,
    data::Bent,
    data::EquilateralTriangle, // 3
    data::VacantTetrahedron,
    data::T,
    data::Tetrahedron, // 4
    data::Square,
    data::Seesaw,
    data::TrigonalPyramid,
    data::SquarePyramid, // 5
    data::TrigonalBipyramid,
    data::Pentagon,
    data::Octahedron, // 6
    data::TrigonalPrism,
    data::PentagonalPyramid,
    data::Hexagon
  >;

  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::TupleType::map<TestTypes, NumUnlinkedTestFunctor>()
    ),
    "Not all numbers of unlinked stereopermutations match across constexpr "
    " and dynamic algorithms"
  );

  BOOST_CHECK(properties::numUnlinkedStereopermutations(Symmetry::Shape::Line, 0) == 1);
  BOOST_CHECK(properties::numUnlinkedStereopermutations(Symmetry::Shape::Bent, 0) == 1);
  BOOST_CHECK(properties::numUnlinkedStereopermutations(Symmetry::Shape::EquilateralTriangle, 0) == 1);
  BOOST_CHECK(properties::numUnlinkedStereopermutations(Symmetry::Shape::Tetrahedron, 0) == 2);
  BOOST_CHECK(properties::numUnlinkedStereopermutations(Symmetry::Shape::Octahedron, 0) == 30);
}

static_assert(
  nShapes == std::tuple_size<data::allShapeDataTypes>::value,
  "nShapes does not equal number of shape data class types in "
  "allShapeDataTypes"
);

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
BOOST_AUTO_TEST_CASE(mappingsAreAvailable) {
  /* In every case where allMappings has a value, getMapping must also return
   * a some optional
   */
  bool pass = true;
  for(const auto& fromShape : Symmetry::allShapes) {
    auto i = static_cast<unsigned>(fromShape);
    for(const auto& toShape : Symmetry::allShapes) {
      auto j = static_cast<unsigned>(toShape);
      if(
        i < j
        && allMappings.at(i, j).hasValue()
          != static_cast<bool>(
            Symmetry::getMapping(fromShape, toShape)
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

BOOST_AUTO_TEST_CASE(angleBoundsTests) {
  BOOST_CHECK(Symmetry::minimumAngle(Symmetry::Shape::T) == M_PI / 2);
  BOOST_CHECK(Symmetry::maximumAngle(Symmetry::Shape::T) == M_PI);

  BOOST_CHECK(Symmetry::minimumAngle(Symmetry::Shape::Octahedron) == M_PI / 2);
  BOOST_CHECK(Symmetry::maximumAngle(Symmetry::Shape::Octahedron) == M_PI);

  BOOST_CHECK(Symmetry::minimumAngle(Symmetry::Shape::TrigonalBipyramid) == M_PI / 2);
  BOOST_CHECK(Symmetry::maximumAngle(Symmetry::Shape::TrigonalBipyramid) == M_PI);

  BOOST_CHECK(
    Symmetry::minimumAngle(Symmetry::Shape::Tetrahedron) == Symmetry::maximumAngle(Symmetry::Shape::Tetrahedron)
  );
}

std::ostream& operator << (std::ostream& os, const std::vector<char>& chars) {
  std::cout << "{";
  for(const char x : chars) {
    std::cout << x;
  }
  std::cout << "}";
  return os;
}

BOOST_AUTO_TEST_CASE(PositionGroups) {
  BOOST_CHECK(properties::positionGroups(Shape::Line) == std::vector<char> (2, 'A'));
  BOOST_CHECK(properties::positionGroups(Shape::Tetrahedron) == std::vector<char> (4, 'A'));
  BOOST_CHECK(properties::positionGroups(Shape::TrigonalBipyramid) == (std::vector<char> {'A','A','A','B','B'}));
  BOOST_CHECK(properties::positionGroups(Shape::CappedTrigonalPrism) == (std::vector<char> {'A','B','C','B','C','D','D'}));
  BOOST_CHECK(properties::positionGroups(Shape::Octahedron) == std::vector<char> (6, 'A'));
  BOOST_CHECK(properties::positionGroups(Shape::Cube) == std::vector<char> (8, 'A'));
  BOOST_CHECK(properties::positionGroups(Shape::Icosahedron) == std::vector<char> (12, 'A'));
}
