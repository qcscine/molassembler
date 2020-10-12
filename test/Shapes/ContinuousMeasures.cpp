/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details. for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Shapes/CoordinateSystemTransformation.h"
#include "Molassembler/Shapes/Shapes.h"
#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Shapes/InertialMoments.h"
#include "Molassembler/Shapes/ContinuousMeasures.h"

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Adaptors/Iota.h"

#include <iostream>
#include "Molassembler/Temple/Stringify.h"

using namespace Scine::Molassembler;
using namespace Shapes;

// From InertialMoments.cpp
extern const std::string& topName(Top top);
extern void randomlyRotate(Eigen::Ref<Continuous::PositionCollection> vs);

inline Vertex operator "" _v (unsigned long long v) { return Vertex(v); }

template<typename EnumType>
constexpr inline auto underlying(const EnumType e) {
  return static_cast<std::underlying_type_t<EnumType>>(e);
}

Continuous::PositionCollection addOrigin(const Continuous::PositionCollection& vs) {
  const unsigned N = vs.cols();
  Continuous::PositionCollection positions(3, N + 1);
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) = vs.col(i);
  }

  // Add origin point explicitly to consideration
  positions.col(N) = Eigen::Vector3d::Zero(3);
  return positions;
}

void distort(Eigen::Ref<Continuous::PositionCollection> positions, const double distortionNorm = 0.01) {
  const unsigned N = positions.cols();
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) += distortionNorm * Eigen::Vector3d::Random().normalized();
  }
}

const std::string& pointGroupString(PointGroup group) {
  static const std::vector<std::string> strings {
    "C1", "Ci", "Cs",
    "C2", "C3", "C4", "C5", "C6", "C7", "C8",
    "C2h", "C3h", "C4h", "C5h", "C6h", "C7h", "C8h",
    "C2v", "C3v", "C4v", "C5v", "C6v", "C7v", "C8v",
    "S4", "S6", "S8",
    "D2", "D3", "D4", "D5", "D6", "D7", "D8",
    "D2h", "D3h", "D4h", "D5h", "D6h", "D7h", "D8h",
    "D2d", "D3d", "D4d", "D5d", "D6d", "D7d", "D8d",
    "T", "Td", "Th",
    "O", "Oh",
    "I", "Ih",
    "Cinfv", "Dinfh"
  };

  return strings.at(
    static_cast<std::underlying_type<PointGroup>::type>(group)
  );
}

BOOST_AUTO_TEST_CASE(PointGroupMeasures, *boost::unit_test::label("Shapes")) {
  unsigned shapesTestedCount = 0;
  unsigned shapesPassCount = 0;

  for(const Shape shape : allShapes) {
#ifndef NDEBUG
    // Skip sizes greater 4 in debug builds
    if(size(shape) > 4) {
      continue;
    }
#endif

    // Even release builds bow out here
    if(size(shape) > 8) {
      continue;
    }

    auto positions = addOrigin(coordinates(shape));
    // distort(positions);
    auto normalized = Continuous::normalize(positions);

    // Add a random coordinate transformation
    normalized = Eigen::Quaterniond::UnitRandom().toRotationMatrix() * normalized;

    // Standardize the top
    const Top top = standardizeTop(normalized);
    if(top == Top::Asymmetric) {
      reorientAsymmetricTop(normalized);
    }

    const double pgCSM = Continuous::pointGroup(
      normalized,
      pointGroup(shape)
    );

    ++shapesTestedCount;

    if(pgCSM >= 0.01) {
      BOOST_TEST_MESSAGE(
        "Expected CSM(" << pointGroupString(pointGroup(shape))
        << ") < 0.01 for " << Shapes::name(shape)
        << "shape, got " << pgCSM << " (top is " << topName(top) << ")"
      );
    } else {
      ++shapesPassCount;
    }
  }

  BOOST_CHECK_MESSAGE(
    shapesTestedCount - shapesPassCount <= 1,
    "A single shape continuous symmetry measure is allowed to fail in the tests, but not multiple."
  );
}

std::ostream& operator << (std::ostream& os, const PointGroup group) {
  const auto elements = Elements::symmetryElements(group);
  const auto groupings = Elements::npGroupings(elements);
  os << pointGroupString(group) << ": {";
  for(const auto& element : elements) {
    os << element -> name() << ", ";
  }
  os << "}\n";
  for(const auto& iterPair : groupings) {
    for(const auto& grouping : iterPair.second) {
      os << "  np = " << iterPair.first << " along " << grouping.probePoint.transpose() << " -> " << Temple::stringify(grouping.groups) << "\n";
      os << "  ";
      os << Temple::stringifyContainer(grouping.groups,
        [&](const auto& grp) -> std::string {
          return Temple::stringifyContainer(grp,
            [&elements](const unsigned elementIdx) -> std::string {
              return elements.at(elementIdx)->name();
            }
          );
        }
      ) << "\n";
    }
  }

  return os;
}

BOOST_AUTO_TEST_CASE(PointGroupElementGroupings, *boost::unit_test::label("Shapes")) {
  const PointGroup limit = PointGroup::Ih;
  for(unsigned g = 0; g <= underlying(limit); ++g) {
    const PointGroup pointGroup = static_cast<PointGroup>(g);
    const auto elements = Elements::symmetryElements(pointGroup);
    const auto groupings = Elements::npGroupings(elements);

    BOOST_CHECK_MESSAGE(
      groupings.count(1) > 0,
      "There is no single-point element group for point group " << pointGroupString(pointGroup)
    );

    bool anyGroupSizeMismatches = false;
    for(const auto& sizeGroupingsPair : groupings) {
      BOOST_CHECK_MESSAGE(
        elements.size() % sizeGroupingsPair.first == 0,
        "Grouping does not evenly divide the group "
        << pointGroupString(pointGroup) << ", G = " << elements.size()
        << ", group size = " << sizeGroupingsPair.first
      );

      bool pass = Temple::all_of(
        sizeGroupingsPair.second,
        [&](const Elements::ElementGrouping& grouping) -> bool {
          const auto size = grouping.groups.front().size();
          return Temple::all_of(
            grouping.groups,
            [&](const auto& groupSubVector) -> bool {
              return groupSubVector.size() == size;
            }
          );
        }
      );

      if(!pass) {
        anyGroupSizeMismatches = true;
      }

      BOOST_CHECK_MESSAGE(
        pass,
        "Not all subgroups of " << pointGroupString(pointGroup) << " have the same size!"
      );
    }
    if(anyGroupSizeMismatches) {
      std::cout << static_cast<PointGroup>(g);
    }
  }
}

// BOOST_AUTO_TEST_CASE(PaperCSMExamples, *boost::unit_test::label("Shapes")) {
//   /* From Continuous Symmetry Measures. 2. Symmetry Groups and the Tetrahedron,
//    * Zabrodsky, Peleg, Avnir. J. Am. Chem. Soc. 1993, 115, 8278-8289
//    *
//    * Example on p. 8283 from their ref 12
//    */
//   PositionCollection tetrahedralPositions(3, 4);
//   tetrahedralPositions.col(0) = Eigen::Vector3d {0.0, 0.0, 1.645};
//   tetrahedralPositions.col(1) = Eigen::Vector3d {0.0, 1.518860, -0.347028};
//   tetrahedralPositions.col(2) = Eigen::Vector3d {-1.286385, -0.700083, -0.391603};
//   tetrahedralPositions.col(3) = Eigen::Vector3d {1.179085, -0.755461, -0.372341};
//
//   const double expectedTdCSM = 0.17;
//   const double calculatedTd = csm::pointGroup(
//     Detail::normalize(tetrahedralPositions),
//     PointGroup::Td
//   ).value_or(1000);
//
//   std::cout << "Calculated Td CSM: " << calculatedTd << " (expect " << expectedTdCSM << ")\n";
//   BOOST_CHECK_CLOSE(expectedTdCSM, calculatedTd, 1);
//
//   /* Example from page 8284 */
//   PositionCollection c3vPositions(3, 4);
//   c3vPositions.col(0) = Eigen::Vector3d {0.0, 0.0, 1.645};
//   c3vPositions.col(1) = Eigen::Vector3d {0.0, 0.87461971, -0.48460962};
//   c3vPositions.col(2) = Eigen::Vector3d {0.75128605, -0.39338892, -0.52991926};
//
//   /* NOTE: Sign inversion on x is speculative, paper reprints P3 = P4 positions
//    * although this does not make sense to me.
//    */
//   c3vPositions.col(3) = Eigen::Vector3d {-0.75128605, -0.39338892, -0.52991926};
//
//   const double expectedC3vCSM = 1.16;
//   const double calculatedC3v = csm::pointGroup(
//     Detail::normalize(c3vPositions),
//     PointGroup::C3v
//   ).value_or(1000);
//
//   std::cout << "Calculated C3v CSM: " << calculatedC3v << " (expect " << expectedC3vCSM << ")\n";
//   BOOST_CHECK_CLOSE(expectedC3vCSM, calculatedC3v, 1);
// }

BOOST_AUTO_TEST_CASE(SquareC4D4PointGroups, *boost::unit_test::label("Shapes")) {
  const double pointGroupCSM = Continuous::pointGroup(
    Continuous::normalize(coordinates(Shape::Square)),
    PointGroup::C4
  );
  BOOST_CHECK_MESSAGE(
    std::fabs(pointGroupCSM) < 1e-10,
    "C4 point group CSM on square planar coordinates is not zero, but " << pointGroupCSM
  );

  const double D4CSM = Continuous::pointGroup(
    Continuous::normalize(coordinates(Shape::Square)),
    PointGroup::D4
  );

  BOOST_CHECK_MESSAGE(
    0 < D4CSM && D4CSM < 1e-10,
    "D4 CSM on square planar coordinates is not zero, but " << D4CSM
  );
}

BOOST_AUTO_TEST_CASE(FixedCnAxis, *boost::unit_test::label("Shapes")) {
  const std::vector<std::pair<Shape, unsigned>> highestOrderAxis {
    {Shape::Bent, 2},
    {Shape::EquilateralTriangle, 3},
    {Shape::VacantTetrahedron, 3},
    {Shape::T, 2},
    {Shape::Tetrahedron, 3},
    {Shape::Square, 4},
    {Shape::Seesaw, 2},
    {Shape::TrigonalPyramid, 3},
    {Shape::SquarePyramid, 4},
    {Shape::TrigonalBipyramid, 3},
    {Shape::Pentagon, 5},
    {Shape::Octahedron, 4},
    {Shape::TrigonalPrism, 3},
    {Shape::PentagonalPyramid, 5},
    {Shape::PentagonalBipyramid, 5},
    {Shape::SquareAntiprism, 4}
  };

  const Eigen::Matrix3d axes = Eigen::Matrix3d::Identity();

  constexpr double acceptanceThreshold = 0.3;

  for(const auto& nameOrderPair : highestOrderAxis) {
    auto positions = addOrigin(coordinates(nameOrderPair.first));
    distort(positions);

    auto normalizedPositions = Continuous::normalize(positions);
    standardizeTop(normalizedPositions);

    boost::optional<unsigned> highestFoundOrder;
    for(unsigned i = 0; i < 3; ++i) {
      Eigen::Vector3d axis = axes.col(i);
      for(unsigned n = 2; n < 6; ++n) {
        const double cn = Continuous::Fixed::element(normalizedPositions, Elements::Rotation::Cn(axis, n));

        if(cn < acceptanceThreshold) {
          highestFoundOrder = n;
        }
      }
    }

    BOOST_CHECK_MESSAGE(
      highestFoundOrder,
      "No Cn axis found along any principal moment axis for "
      << Shapes::name(nameOrderPair.first)
    );
    if(highestFoundOrder) {
      BOOST_CHECK_MESSAGE(
        highestFoundOrder.value() == nameOrderPair.second,
        "Expected to find Cn of order " << nameOrderPair.second << ", but found " << highestFoundOrder.value() << " instead."
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(AlleneS4, *boost::unit_test::label("Shapes")) {
  Continuous::PositionCollection allenePositions(3, 7);
  allenePositions << 0.0, 0.0, 0.0, 0.0, 0.0, 0.928334, -0.928334,
                     0.0, 0.0, 0.0, 0.928334, -0.928334, 0.0, 0.0,
                     0.0, 1.3102, -1.3102, 1.866201, 1.866201, -1.866201, -1.866201;

  auto normalizedPositions = Continuous::normalize(allenePositions);
  const Top top = standardizeTop(normalizedPositions);
  BOOST_CHECK_MESSAGE(
    top == Top::Prolate,
    "Top isn't prolate, but " << topName(top)
  );

  const double S4CSM = Continuous::Fixed::element(normalizedPositions, Elements::Rotation::Sn(Eigen::Vector3d::UnitZ(), 4));
  BOOST_CHECK_MESSAGE(
    S4CSM < 0.1,
    "CSM(S4) = " << S4CSM << " of allene is over recognition threshold (0.1)"
  );

  const double D2dCSM = Continuous::pointGroup(normalizedPositions, PointGroup::D2d);
  BOOST_CHECK_MESSAGE(
    D2dCSM < 0.1,
    "CSM(D2d) = " << D2dCSM << " of allene is over recognition threshold (0.1)"
  );

  const auto optimizedS4Result = Continuous::element(
    normalizedPositions,
    Elements::Rotation::Sn(
      Eigen::Vector3d::UnitZ() + 0.1 * Eigen::Vector3d::Random().normalized(),
      4
    )
  );

  BOOST_CHECK_MESSAGE(
    optimizedS4Result.second.axis.isApprox(Eigen::Vector3d::UnitZ(), 1e-2),
    "Axis of optimized S4 is not +z, but " << optimizedS4Result.second.axis.transpose()
  );

  BOOST_CHECK_LT(optimizedS4Result.first, 0.1);
}

BOOST_AUTO_TEST_CASE(ReflectionPlaneOptimization, *boost::unit_test::label("Shapes")) {
  // Generate 8 points in the xy plane
  Continuous::PositionCollection planarPositions(3, 8);
  for(unsigned i = 0; i < 8; ++i) {
    Eigen::Vector3d v = 3 * Eigen::Vector3d::Random();
    v.z() = 0;
    planarPositions.col(i) = v;
  }

  auto normalized = Continuous::normalize(planarPositions);

  const double zPlaneCSM = Continuous::Fixed::element(
    normalized,
    Elements::Reflection {Eigen::Vector3d::UnitZ()}
  );

  BOOST_CHECK_LT(zPlaneCSM, 0.1);

  const auto optimizedSigma = Continuous::element(
    normalized,
    Elements::Reflection {
      Eigen::Vector3d::UnitZ() + 0.1 * Eigen::Vector3d::Random().normalized()
    }
  );

  BOOST_CHECK_MESSAGE(
    optimizedSigma.second.normal.isApprox(Eigen::Vector3d::UnitZ(), 1e-2),
    "Optimized sigma plane's normal is not +z, but " << optimizedSigma.second.normal.transpose()
  );

  BOOST_CHECK_LT(optimizedSigma.first, 0.1);
}

BOOST_AUTO_TEST_CASE(AsymmetricTopStandardization, *boost::unit_test::label("Shapes")) {
  std::vector<Shape> asymmetricTopsWithC2 {
    Shape::Bent,
    Shape::T,
    Shape::Seesaw
  };

  for(const Shape shape : asymmetricTopsWithC2) {
    auto shapeCoordinates = addOrigin(coordinates(shape));
    auto normalizedPositions = Continuous::normalize(shapeCoordinates);
    const Top top = standardizeTop(normalizedPositions);
    BOOST_CHECK_MESSAGE(
      top == Top::Asymmetric,
      "Expected asymmetric top for " << name(shape) << ", got "
      << topName(top) << " instead"
    );

    const unsigned highestAxisOrder = reorientAsymmetricTop(normalizedPositions);
    BOOST_CHECK_EQUAL(highestAxisOrder, 2);

    // Ensure rotation of highest order axis to z worked
    const double CnCSM = Continuous::Fixed::element(
      normalizedPositions,
      Elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), 2)
    );
    BOOST_CHECK_MESSAGE(
      CnCSM < 1e-10,
      "Expected Cn of order 2 along z < 1e-10, got " << CnCSM << " instead."
    );
  }
}

BOOST_AUTO_TEST_CASE(ShapeMeasuresYieldForwardPermutation, *boost::unit_test::label("Shapes")) {
  const Shape testShape = Shape::Tetrahedron;

  auto shapeCoordinates = Continuous::normalize(
    addOrigin(coordinates(testShape))
  );

  Continuous::PositionCollection shuffled (3, shapeCoordinates.cols());
  BOOST_REQUIRE_MESSAGE(size(testShape) == 4, "Test setup is no longer valid");
  const std::vector<Vertex> shufflePermutation {{3_v, 4_v, 0_v, 2_v, 1_v}};
  for(unsigned i = 0; i < 5; ++i) {
    shuffled.col(shufflePermutation.at(i)) = shapeCoordinates.col(i);
  }

  auto faithful = Continuous::shapeFaithfulPaperImplementation(shuffled, testShape);
  auto alternate = Continuous::shapeAlternateImplementation(shuffled, testShape);
  auto heuristics = Continuous::shapeHeuristics(shuffled, testShape);

  BOOST_CHECK_CLOSE(faithful.measure, alternate.measure, 1);
  BOOST_CHECK_CLOSE(faithful.measure, heuristics.measure, 1);

  auto isBackwardsMapping = [](const std::vector<Vertex>& permutation) -> bool {
    /* A forward permutation maps from positions to shape indices.
     * A backwards permutation maps from shape indices to positions.
     *
     * Only backwards permutations have the (added) centroid at the back, and
     * hence their [0, S] = [0, N - 1] = [0, N) ranges are sub-permutations
     * (if ordered, identical to iota(N)) that can be rotated with shape
     * rotations.
     */
    return permutation.back() == permutation.size() - 1;
  };

  auto invert = [](const std::vector<Vertex>& p) -> std::vector<Vertex> {
    const unsigned P = p.size();
    std::vector<Vertex> inverse(P);
    for(Vertex i {0}; i < P; ++i) {
      inverse.at(p.at(i)) = i;
    }
    return inverse;
  };

  BOOST_CHECK(shufflePermutation == invert(faithful.mapping));
  BOOST_CHECK(shufflePermutation == invert(alternate.mapping));
  BOOST_CHECK(shufflePermutation == invert(heuristics.mapping));

  BOOST_CHECK_MESSAGE(
    !isBackwardsMapping(faithful.mapping),
    "Faithful mapping is not a forward mapping: " << Temple::condense(faithful.mapping)
  );

  BOOST_CHECK_MESSAGE(
    !isBackwardsMapping(alternate.mapping),
    "Alternate mapping is not a forward mapping: " << Temple::condense(alternate.mapping)
  );

  BOOST_CHECK_MESSAGE(
    !isBackwardsMapping(heuristics.mapping),
    "Heuristics mapping is not a forward mapping: " << Temple::condense(heuristics.mapping)
  );

  BOOST_CHECK_MESSAGE(
    faithful.mapping == alternate.mapping,
    "Faithful mapping does not match alternate mapping.\n"
    << "Faithful: " << Temple::condense(faithful.mapping) << "\n"
    << "Alternate: " << Temple::condense(alternate.mapping)
  );

  BOOST_CHECK_MESSAGE(
    faithful.mapping == heuristics.mapping,
    "Faithful mapping does not match heuristics mapping.\n"
    << "Faithful: " << Temple::condense(faithful.mapping) << "\n"
    << "Heuristics: " << Temple::condense(heuristics.mapping)
  );
}

BOOST_AUTO_TEST_CASE(ShapeMeasuresAlternateAlgorithm, *boost::unit_test::label("Shapes")) {
#ifdef NDEBUG
  constexpr unsigned testingShapeSizeLimit = 7;
#else
  constexpr unsigned testingShapeSizeLimit = 5;
#endif

  for(const Shape shape : allShapes) {
    if(size(shape) > testingShapeSizeLimit) {
      continue;
    }

    auto shapeCoordinates = Continuous::normalize(
      addOrigin(coordinates(shape))
    );
    const double unrotated = Continuous::shapeAlternateImplementation(shapeCoordinates, shape).measure;
    BOOST_CHECK_MESSAGE(
      unrotated < 1e-10,
      "Expected CShM < 1e-10 for unrotated coordinates of " << name(shape) << ", but got " << unrotated
    );
    randomlyRotate(shapeCoordinates);
    const double rotated = Continuous::shapeAlternateImplementation(shapeCoordinates, shape).measure;
    BOOST_CHECK_MESSAGE(
      rotated < 0.1,
      "Expected CShM < 1e-2 for rotated coordinates of " << name(shape) << ", but got " << rotated
    );

    for(unsigned i = 1; i < 5; ++i) {
      const double distortionNorm = 0.1 * i;
      auto distorted = shapeCoordinates;
      distort(distorted, distortionNorm);
      distorted = Continuous::normalize(distorted);

      const double faithful = Continuous::shapeFaithfulPaperImplementation(
        distorted,
        shape
      ).measure;
      const double alternate = Continuous::shapeAlternateImplementation(distorted, shape).measure;
      BOOST_CHECK_CLOSE(faithful, alternate, 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(ShapeMeasuresHeuristics, *boost::unit_test::label("Shapes")) {
#ifdef NDEBUG
  constexpr unsigned testingShapeSizeLimit = 7;
#else
  constexpr unsigned testingShapeSizeLimit = 5;
#endif

  constexpr unsigned repeats = 5;

  for(const Shape shape : allShapes) {
    // Cannot use heuristics on fewer than five vertices
    if(size(shape) < 4) {
      continue;
    }

    if(size(shape) > testingShapeSizeLimit) {
      continue;
    }

    auto shapeCoordinates = Continuous::normalize(
      addOrigin(coordinates(shape))
    );

    for(unsigned i = 1; i < 6; ++i) {
      const double distortionNorm = 0.1 * i;

      std::vector<double> referenceValues;
      std::vector<double> errors;
      for(unsigned j = 0; j < repeats; ++j) {
        auto distorted = shapeCoordinates;
        distort(distorted, distortionNorm);
        distorted = Continuous::normalize(distorted);

        const double alternate = Continuous::shapeAlternateImplementation(distorted, shape).measure;
        referenceValues.push_back(alternate);
        const double heuristic = Continuous::shapeHeuristics(distorted, shape).measure;
        const double error = std::fabs(alternate - heuristic);

        if(error > 1e-2) {
          std::cout << "Error " << error << " for shape " << name(shape) << " and distortion " << distortionNorm << "\n";
        }
        errors.push_back(std::fabs(alternate - heuristic));
      }

      const double referenceAverage = Temple::accumulate(referenceValues, 0.0, std::plus<>()) / repeats;
      const double errorAverage = Temple::accumulate(errors, 0.0, std::plus<>()) / repeats;
      BOOST_CHECK_MESSAGE(
        errorAverage < 0.01 * referenceAverage,
        "Expected error average below 1% of reference value average, but mu(error) = "
        << errorAverage << " >= 0.01 * " << referenceAverage << " = "
        << (0.01 * referenceAverage) << " for shape " << name(shape) << " and distortion norm " << distortionNorm
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(MinimumDistortionConstants, *boost::unit_test::label("Shapes")) {
  /* NOTES
   * - These constants are from https://pubs.acs.org/doi/10.1021/ja036479n
   *
   *   k_XY = sqrt( CShM_X(Y) ) = sqrt( CShM_Y(X) )
   *
   *   where CShM_X(Y) is the shape measure regarding shape X of coordinates of
   *   the ideal shape Y
   *
   * - VOC-5 is our square pyramid, their square pyramid has 105° base-apex angles
   */

  const std::vector<
    std::tuple<Shape, Shape, double>
  > minimumDistortionConstants {
    {Shape::Tetrahedron, Shape::Square, 5.774}, // T-4, SP-4

    /* NOTE: Seesaw / Sawhorse angle of "equatorial" vertices can vary,
     * unclear here from the paper and close, but mismatching with our
     * definition (120°, essentially a single-vacant trigonal bipyramid)
     */
//    {Shape::Tetrahedron, Shape::Seesaw, 3.129}, // T-4, SW-4 [sic, is SS-4]
//    {Shape::Seesaw, Shape::Square, 4.365}, // SW-4, SP-4

    {Shape::SquarePyramid, Shape::TrigonalBipyramid, 2.710}, // VOC-5, TBPY-5
    {Shape::SquarePyramid, Shape::Pentagon, 5.677}, // VOC-5, PP-5
    {Shape::TrigonalBipyramid, Shape::Pentagon, 6.088}, // TBPY-5, PP-5

    // NOTE: All of the below that are commented out are missing the required shapes

    {Shape::Octahedron, Shape::TrigonalPrism, 4.091}, // OC-6, TPR-6
    {Shape::Octahedron, Shape::PentagonalPyramid, 5.517}, // OC-6, PPY-6
//    {Shape::Octahedron, Shape::Hexagon, 5.774}, // OC-6, HP-6
    {Shape::TrigonalPrism, Shape::PentagonalPyramid, 4.125}, // TPR-6, PPY-6
//    {Shape::TrigonalPrism, Shape::Hexagon, 5.803}, // TPR-6, HP-6
//    {Shape::PentagonalPyramid, Shape::Hexagon, 5.352}, // PPY-6, HP-6


//    {Shape::CappedOctahedron, Shape::CappedTrigonalPrism, 1.236},
//    {Shape::CappedOctahedron, Shape::PentagonalBipyramid, 2.899},
//    {Shape::CappedOctahedron, Shape::HexagonalPyramid, 4.130},
//    {Shape::CappedOctahedron, Shape::Heptagon, 6.146},
//    {Shape::CappedTrigonalPrism, Shape::PentagonalBipyramid, 2.577},
//    {Shape::CappedTrigonalPrism, Shape::HexagonalPyramid, 4.467},
//    {Shape::CappedTrigonalPrism, Shape::Heptagon, 5.989},
//    {Shape::PentagonalBipyramid, Shape::HexagonalPyramid, 5.166},
//    {Shape::PentagonalBipyramid, Shape::Heptagon, 5.934},
//    {Shape::HexagonalPyramid, Shape::Heptagon, 5.047},

//    {Shape::Cube, Shape::TrigonalDodecahedron, 2.820},
//    {Shape::Cube, Shape::SquareAntiprism, 3.315},
//    {Shape::Cube, Shape::HexagonalBipyramid, 2.897},
//    {Shape::Cube, Shape::HeptagonalPyramid, 5.533},
//    {Shape::TrigonalDodecahedron, Shape::SquareAntiprism, 1.688},
//    {Shape::TrigonalDodecahedron, Shape::HexagonalBipyramid, 3.960},
//    {Shape::TrigonalDodecahedron, Shape::HeptagonalPyramid, 4.989},
//    {Shape::SquareAntiprism, Shape::HexagonalBipyramid, 4.296},
//    {Shape::SquareAntiprism, Shape::HeptagonalPyramid, 4.953},
//    {Shape::HexagonalBipyramid, Shape::HeptagonalPyramid, 4.865},
  };

  auto testF = [](const Shape a, const Shape b, const double expectedMinimumDistortion) {
    auto makeShapeCoordinates = [](const Shape shape) {
      return Continuous::normalize(
        addOrigin(coordinates(shape))
      );
    };

    const double ab = Continuous::shapeAlternateImplementation(makeShapeCoordinates(a), b).measure;
    const double ba = Continuous::shapeAlternateImplementation(makeShapeCoordinates(b), a).measure;

    BOOST_CHECK_CLOSE(ab, ba, 0.1);

    const double calculatedMinimumDistortion = std::sqrt((ab + ba) / 2);
    BOOST_CHECK_CLOSE(calculatedMinimumDistortion, expectedMinimumDistortion, 2);
  };

  Temple::forEach(minimumDistortionConstants, testF);
}
