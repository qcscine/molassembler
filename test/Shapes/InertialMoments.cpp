/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Shapes/CoordinateSystemTransformation.h"
#include "Molassembler/Shapes/InertialMoments.h"
#include "Molassembler/Shapes/ContinuousMeasures.h"
#include "Molassembler/Shapes/Shapes.h"
#include "Molassembler/Shapes/Data.h"

using namespace Scine::Molassembler;
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

BOOST_AUTO_TEST_CASE(InertialStandardization, *boost::unit_test::label("Shapes")) {
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
