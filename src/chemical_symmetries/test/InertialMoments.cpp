/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "chemical_symmetries/CoordinateSystemTransformation.h"
#include "chemical_symmetries/InertialMoments.h"
#include "chemical_symmetries/ContinuousMeasures.h"
#include "chemical_symmetries/Shapes.h"
#include "chemical_symmetries/Symmetries.h"

using namespace Scine;
using namespace Symmetry;

// From ContinuousMeasures.cpp
extern continuous::PositionCollection addOrigin(const continuous::PositionCollection& vs);

void randomlyRotate(Eigen::Ref<continuous::PositionCollection> vs) {
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
    {Shape::Disphenoid, Top::Asymmetric},
    {Shape::TrigonalBipyramid, Top::Prolate},
    {Shape::PentagonalBipyramid, Top::Oblate},
    {Shape::Square, Top::Oblate},
    {Shape::Octahedron, Top::Spherical},
    {Shape::Tetrahedron, Top::Spherical}
  };

  for(const auto nameTopPair : tops) {
    auto positions = addOrigin(symmetryData().at(nameTopPair.first).coordinates);

    // Apply a random coordinate transformation

    // Analyze it
    auto normalizedPositions = continuous::normalize(positions);
    Top top = standardizeTop(normalizedPositions);
    BOOST_CHECK_MESSAGE(
      top == nameTopPair.second,
      "Top standardization failed. Expected "
      << topName(nameTopPair.second)
      << " for symmetry " << Symmetry::name(nameTopPair.first)
      << ", got " << topName(top) << " instead."
    );
  }
}
