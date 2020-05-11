/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details. for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Shapes/CoordinateSystemTransformation.h"
#include "molassembler/Shapes/InertialMoments.h"
#include "molassembler/Shapes/ContinuousMeasures.h"
#include "molassembler/Shapes/Shapes.h"
#include "molassembler/Shapes/Data.h"

using namespace Scine;
using namespace Shapes;

// From ContinuousMeasures.cpp
extern Continuous::PositionCollection addOrigin(const Continuous::PositionCollection& vs);

void randomlyRotate(Eigen::Ref<Continuous::PositionCollection> vs) {
  vs = rotationMatrix(CoordinateSystem {}, CoordinateSystem::random()) * vs;
}

const std::vector<std::string>& topNames() {
  static const std::vector<std::string> strings {
    "Line",
    "Asymmetric",
    "Prolate",
    "Oblate",
    "Spherical"
  };

  return strings;
}

const std::string& topName(Top top) {
  return topNames().at(
    static_cast<std::underlying_type<Top>::type>(top)
  );
}

BOOST_AUTO_TEST_CASE(InertialStandardization) {
  const std::vector<
    std::pair<Shape, Top>
  > tops {
    {Shape::Line, Top::Line},
    {Shape::Bent, Top::Asymmetric},
    {Shape::Seesaw, Top::Asymmetric},
    {Shape::TrigonalBipyramid, Top::Prolate},
    {Shape::PentagonalBipyramid, Top::Oblate},
    {Shape::Square, Top::Oblate},
    {Shape::Octahedron, Top::Spherical},
    {Shape::Tetrahedron, Top::Spherical}
  };

  for(const auto& nameTopPair : tops) {
    auto positions = addOrigin(coordinates(nameTopPair.first));

    // Apply a random coordinate transformation

    // Analyze it
    auto normalizedPositions = Continuous::normalize(positions);
    Top top = standardizeTop(normalizedPositions);
    BOOST_CHECK_MESSAGE(
      top == nameTopPair.second,
      "Top standardization failed. Expected "
      << topName(nameTopPair.second)
      << " for shape " << Shapes::name(nameTopPair.first)
      << ", got " << topName(top) << " instead."
    );
  }
}
