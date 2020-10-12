/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>
#include <Eigen/Geometry>

#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Functor.h"
#include "Molassembler/Temple/SetAlgorithms.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/Stl17.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/constexpr/ToStl.h"
#include "Molassembler/Temple/constexpr/TupleTypePairs.h"

#include "Molassembler/Shapes/PropertyCaching.h"
#include "Molassembler/Shapes/Data.h"

#include <set>
#include <numeric>
#include <iostream>
#include <iomanip>

using namespace Scine::Molassembler;
using namespace Shapes;

template<typename EnumType>
constexpr inline auto underlying(const EnumType e) {
  return static_cast<std::underlying_type_t<EnumType>>(e);
}

template<typename ShapeClass>
struct LockstepTest {
  static bool value() {
    return Shapes::nameIndex(ShapeClass::shape) == underlying(ShapeClass::shape);
  }
};

BOOST_AUTO_TEST_CASE(SymmetryTypeAndPositionInEnumLockstep, *boost::unit_test::label("Shapes")) {
  BOOST_CHECK_MESSAGE(
    Temple::all_of(
      Temple::Tuples::map<
        Shapes::Data::allShapeDataTypes,
        LockstepTest
      >()
    ),
    "Not all symmetries have the same order in Name and allShapeDataTypes"
  );
}

BOOST_AUTO_TEST_CASE(AngleFunctionInputSymmetry, *boost::unit_test::label("Shapes")) {
  // every angle function must be symmetrical on input of valid unsigned indices
  for(const Shape shape: allShapes) {
    bool passesAll = true;

    for(Vertex i(0); i < size(shape) && passesAll; ++i) {
      for(Vertex j(i + 1); j < size(shape); ++j) {
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

BOOST_AUTO_TEST_CASE(AngleFunctionZeroForIdenticalInput, *boost::unit_test::label("Shapes")) {
  // every angle function must return 0 for identical indices
  for(const Shape shape: allShapes) {
    bool passesAll = true;

    for(Vertex i(0); i < size(shape); i++) {
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

BOOST_AUTO_TEST_CASE(AnglesWithinRadiansBounds, *boost::unit_test::label("Shapes")) {
  for(const Shape shape : allShapes) {
    bool passesAll = true;

    for(Vertex i(0); i < size(shape); i++) {
      for(Vertex j(0); j < size(shape); j++) {
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

BOOST_AUTO_TEST_CASE(AnglesMatchCoordinates, *boost::unit_test::label("Shapes")) {

  /* The results of the angle functions ought to match the geometries specified
   * by the coordinates
   */

  for(const auto& shape: allShapes) {
    auto getCoordinates =  [shape](const Vertex index) -> Eigen::Vector3d {
      return coordinates(shape).col(index);
    };

    bool all_pass = true;

    for(Vertex i {0}; i < size(shape); i++) {
      for(Vertex j {i + 1}; j < size(shape); j++) {
        const double angleInCoordinates = std::acos(
          Temple::Stl17::clamp(
            getCoordinates(i).dot(getCoordinates(j)) / (
              getCoordinates(i).norm() * getCoordinates(j).norm()
            ),
            -1.0,
            1.0
          )
        );

        const double angleDifference = angleInCoordinates - angleFunction(shape)(i, j);

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

BOOST_AUTO_TEST_CASE(AllTetrahedraPositive, *boost::unit_test::label("Shapes")) {
  /* Checks if sequence that tetrahedra are defined in leads to a positive
   * volume when calculated via
   *
   *  (1 - 4) dot [ (2 - 4) x (3 - 4) ]
   *
   */
  for(const auto& shape: allShapes) {
    auto getCoordinates = [shape](const boost::optional<Vertex>& indexOption) -> Eigen::Vector3d {
      if(indexOption) {
        return coordinates(shape).col(indexOption.value());
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

BOOST_AUTO_TEST_CASE(TetrahedraDefinitionIndicesUnique, *boost::unit_test::label("Shapes")) {
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

BOOST_AUTO_TEST_CASE(SmallestAngleValueCorrect, *boost::unit_test::label("Shapes")) {
  auto shapeSmallestAngle = [](const Shape shape) -> double {
    return Temple::accumulate(
      Temple::Adaptors::allPairs(
        Temple::iota<Vertex>(Vertex(size(shape)))
      ),
      M_PI,
      [shape](const double currentSmallest, const auto& vertexPair) {
        return std::min(
          currentSmallest,
          Temple::invoke(angleFunction(shape), vertexPair)
        );
      }
    );
  };

  const double comparisonSmallestAngle = Temple::min(
    Temple::map(allShapes, shapeSmallestAngle)
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
  auto dynamicMappings = Properties::selectBestTransitionMappings(
    Properties::shapeTransitionMappings(
      ShapeClassFrom::shape,
      ShapeClassTo::shape
    )
  );

  Temple::Floating::ExpandedRelativeEqualityComparator<double> comparator {
    Properties::floatingPointEqualityThreshold
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
  auto convertedMappings = Temple::map_stl(
    Temple::toSTL(constexprMappings.mappings),
    [&](const auto& indexList) -> std::vector<Vertex> {
      std::vector<Vertex> v;
      for(unsigned i : indexList) {
        v.emplace_back(i);
      }
      return v;
    }
  );

  decltype(convertedMappings) dynamicResultSet {
    std::begin(dynamicMappings.indexMappings),
    std::end(dynamicMappings.indexMappings)
  };

  return Temple::set_symmetric_difference(
    convertedMappings,
    dynamicResultSet
  ).empty();
}

using IndexAndMappingsPairType = std::pair<
  unsigned,
  Shapes::ConstexprProperties::MappingsReturnType
>;

constexpr bool pairEqualityComparator(
  const IndexAndMappingsPairType& a,
  const IndexAndMappingsPairType& b
) {
  Temple::Floating::ExpandedRelativeEqualityComparator<double> comparator {
    Shapes::Properties::floatingPointEqualityThreshold
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
  Temple::Array<
    std::pair<
      unsigned,
      Shapes::ConstexprProperties::MappingsReturnType
    >,
    ShapeClassFrom::size
  > constexprMappings;

  for(unsigned i = 0; i < ShapeClassFrom::size; ++i) {
    constexprMappings.at(i) = std::make_pair(
      i,
      Shapes::ConstexprProperties::ligandLossMappings<
        ShapeClassFrom,
        ShapeClassTo
      >(i)
    );
  }

  // Group the results
  auto constexprGroups = Temple::groupByEquality(
    constexprMappings,
    pairEqualityComparator // C++17 constexpr lambda
  );

  /* Dynamic part */
  std::vector<
    std::pair<unsigned, Shapes::Properties::ShapeTransitionGroup>
  > dynamicMappings;

  for(unsigned i = 0; i < ShapeClassFrom::size; ++i) {
    dynamicMappings.emplace_back(
      i,
      selectBestTransitionMappings(
        Properties::ligandLossTransitionMappings(
          ShapeClassFrom::shape,
          ShapeClassTo::shape,
          i
        )
      )
    );
  }

  // Analyze all mappings - which indices have "identical" target mappings?
  auto dynamicGroups = Temple::groupByEquality(
    dynamicMappings,
    [&](const auto& firstMappingPair, const auto& secondMappingPair) -> bool {
      return (
        Temple::Floating::isCloseRelative(
          firstMappingPair.second.angularDistortion,
          secondMappingPair.second.angularDistortion,
          Shapes::Properties::floatingPointEqualityThreshold
        ) && Temple::Floating::isCloseRelative(
          firstMappingPair.second.chiralDistortion,
          secondMappingPair.second.chiralDistortion,
          Shapes::Properties::floatingPointEqualityThreshold
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
    auto constexprRotations = ConstexprProperties::generateAllRotations<ShapeClass>(
      ConstexprProperties::startingIndexSequence<ShapeClass>()
    );

    // This is a std::set of ShapeClass-sized std::vectors
    auto dynamicRotations = Properties::generateAllRotations(
      ShapeClass::shape,
      Temple::iota<Vertex>(Vertex(ShapeClass::size))
    );

    auto convertedRotations = Temple::map_stl(
      Temple::toSTL(constexprRotations),
      [&](const auto& indexList) -> std::vector<Vertex> {
        std::vector<Vertex> v;
        for(unsigned i : indexList) {
          v.emplace_back(i);
        }
        return v;
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
        Temple::set_symmetric_difference(
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
      std::cout << " Maximum #rotations: " << ConstexprProperties::maxRotations<ShapeClass>()
        << std::endl;

      std::cout << " Converted constexpr:" << std::endl;
      for(const auto& element : convertedRotations) {
        std::cout << " {" << Temple::condense(element)
          << "}\n";
      }

      std::cout << " Dynamic:" << std::endl;
      for(const auto& element : dynamicRotations) {
        std::cout << " {" << Temple::condense(element)
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
 * value (generateAllRotations -> shapeData), the value-initialize variant
 * leads to a static initialization fiasco, where it is unclear whether the
 * value data member or shapeData is initialized first.
 *
 * Only in the case of static constexpr is the value-initialize variant
 * semantically equivalent to the data member variant.
 */
};

std::string getGraphvizNodeName(const Shapes::Shape shape) {
  auto stringName = Shapes::name(shape);

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

BOOST_AUTO_TEST_CASE(ConstexprRotationsMatch, *boost::unit_test::label("Shapes")) {
  // Full test of rotation algorithm equivalency for all symmetries
  BOOST_CHECK_MESSAGE(
    Temple::all_of(
      Temple::Tuples::map<
        Shapes::Data::allShapeDataTypes,
        RotationGenerationTest
      >()
    ),
    "There is a discrepancy between constexpr and dynamic rotation generation"
  );

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
  // Test transitions generation/evaluation algorithm equivalency for all
  BOOST_CHECK_MESSAGE(
    Temple::all_of(
      Temple::Tuples::mapAllPairs<
        Shapes::Data::allShapeDataTypes,
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
      unsigned constexprResult = ConstexprProperties::numUnlinkedStereopermutations<ShapeClass>(i);

      unsigned dynamicResult = Properties::numUnlinkedStereopermutations(ShapeClass::shape, i);

      if(constexprResult != dynamicResult) {
        std::cout << "Mismatch for " << Shapes::name(ShapeClass::shape) << " and " << i << " identical ligands between constexpr and dynamic number of unlinked: " << constexprResult << " vs. " << dynamicResult << "\n";
        return false;
      }

      // Cross-check with constexpr hasMultiple
      bool constexprHasMultiple = ConstexprProperties::hasMultipleUnlinkedStereopermutations<ShapeClass>(i);
      if((constexprResult > 1) != constexprHasMultiple) {
        std::cout << "Mismatch between constexpr count and constexpr "
          << "hasMultiple unlinked ligands for "
          << Shapes::name(ShapeClass::shape) << " and "
          << i << " identical ligands: " << constexprResult << " and "
          << std::boolalpha << constexprHasMultiple << "\n";
        return false;
      }

      // Cross-check with dynamic hasMultiple
      bool dynamicHasMultiple = Properties::hasMultipleUnlinkedStereopermutations(ShapeClass::shape, i);
      if((constexprResult > 1) != dynamicHasMultiple) {
        std::cout << "Mismatch between constexpr count and dynamic "
          << "hasMultiple unlinked ligands for "
          << Shapes::name(ShapeClass::shape) << " and "
          << i << " identical ligands: " << constexprResult << " and "
          << std::boolalpha << dynamicHasMultiple << "\n";
        return false;
      }
    }

    return true;
  }
};

BOOST_AUTO_TEST_CASE(NumUnlinkedStereopermsMatch, *boost::unit_test::label("Shapes")) {
  using TestTypes = std::tuple<
    Data::Line,
    Data::Bent,
    Data::EquilateralTriangle, // 3
    Data::VacantTetrahedron,
    Data::T,
    Data::Tetrahedron, // 4
    Data::Square,
    Data::Seesaw,
    Data::TrigonalPyramid,
    Data::SquarePyramid, // 5
    Data::TrigonalBipyramid,
    Data::Pentagon,
    Data::Octahedron, // 6
    Data::TrigonalPrism,
    Data::PentagonalPyramid,
    Data::Hexagon
  >;

  BOOST_CHECK_MESSAGE(
    Temple::all_of(
      Temple::Tuples::map<TestTypes, NumUnlinkedTestFunctor>()
    ),
    "Not all numbers of unlinked stereopermutations match across constexpr "
    " and dynamic algorithms"
  );

  BOOST_CHECK(Properties::numUnlinkedStereopermutations(Shapes::Shape::Line, 0) == 1);
  BOOST_CHECK(Properties::numUnlinkedStereopermutations(Shapes::Shape::Bent, 0) == 1);
  BOOST_CHECK(Properties::numUnlinkedStereopermutations(Shapes::Shape::EquilateralTriangle, 0) == 1);
  BOOST_CHECK(Properties::numUnlinkedStereopermutations(Shapes::Shape::Tetrahedron, 0) == 2);
  BOOST_CHECK(Properties::numUnlinkedStereopermutations(Shapes::Shape::Octahedron, 0) == 30);
}

static_assert(
  nShapes == std::tuple_size<Data::allShapeDataTypes>::value,
  "nShapes does not equal number of shape data class types in "
  "allShapeDataTypes"
);

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
BOOST_AUTO_TEST_CASE(ConstexprMappingsExist, *boost::unit_test::label("Shapes")) {
  /* In every case where allMappings has a value, getMapping must also return
   * a some optional
   */
  bool pass = true;
  for(const auto& fromShape : Shapes::allShapes) {
    auto i = static_cast<unsigned>(fromShape);
    for(const auto& toShape : Shapes::allShapes) {
      auto j = static_cast<unsigned>(toShape);
      if(
        i < j
        && allMappings.at(i, j).hasValue()
          != static_cast<bool>(
            Shapes::getMapping(fromShape, toShape)
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

BOOST_AUTO_TEST_CASE(AngleBoundsCorrect, *boost::unit_test::label("Shapes")) {
  BOOST_CHECK(Shapes::minimumAngle(Shapes::Shape::T) == M_PI / 2);
  BOOST_CHECK(Shapes::maximumAngle(Shapes::Shape::T) == M_PI);

  BOOST_CHECK(Shapes::minimumAngle(Shapes::Shape::Octahedron) == M_PI / 2);
  BOOST_CHECK(Shapes::maximumAngle(Shapes::Shape::Octahedron) == M_PI);

  BOOST_CHECK(Shapes::minimumAngle(Shapes::Shape::TrigonalBipyramid) == M_PI / 2);
  BOOST_CHECK(Shapes::maximumAngle(Shapes::Shape::TrigonalBipyramid) == M_PI);

  BOOST_CHECK(
    Shapes::minimumAngle(Shapes::Shape::Tetrahedron) == Shapes::maximumAngle(Shapes::Shape::Tetrahedron)
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

BOOST_AUTO_TEST_CASE(PositionGroups, *boost::unit_test::label("Shapes")) {
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::Line) == std::vector<char> (2, 'A'));
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::Tetrahedron) == std::vector<char> (4, 'A'));
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::TrigonalBipyramid) == (std::vector<char> {'A','A','A','B','B'}));
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::CappedTrigonalPrism) == (std::vector<char> {'A','B','C','B','C','D','D'}));
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::Octahedron) == std::vector<char> (6, 'A'));
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::Cube) == std::vector<char> (8, 'A'));
  BOOST_CHECK(Properties::positionGroupCharacters(Shape::Icosahedron) == std::vector<char> (12, 'A'));
}

BOOST_AUTO_TEST_CASE(DimensionalityProperty, *boost::unit_test::label("Shapes")) {
  BOOST_CHECK(threeDimensional(Shape::Tetrahedron));
  BOOST_CHECK(!threeDimensional(Shape::Line));
  BOOST_CHECK(!threeDimensional(Shape::Bent));
  BOOST_CHECK(!threeDimensional(Shape::Pentagon));
  BOOST_CHECK(!threeDimensional(Shape::Hexagon));
  BOOST_CHECK(threeDimensional(Shape::Icosahedron));
}
