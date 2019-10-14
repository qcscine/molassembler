/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "chemical_symmetries/Flowchart.h"
#include "chemical_symmetries/InertialMoments.h"
#include "chemical_symmetries/Symmetries.h"

#include <iostream>

using namespace Scine;
using namespace Symmetry;

extern const std::string& pointGroupString(PointGroup group);
extern continuous::PositionCollection addOrigin(const continuous::PositionCollection& vs);

BOOST_AUTO_TEST_CASE(Flowcharting) {
  for(const Shape shape : allShapes) {
    if(size(shape) >= 4) {
      break;
    }

    auto normalizedCoordinates = continuous::normalize(
      addOrigin(symmetryData().at(shape).coordinates)
    );
    standardizeTop(normalizedCoordinates);

    std::cout << "Flowcharting " << name(shape) << "\n";
    const PointGroup expectedPointGroup = pointGroup(shape);
    const auto flowchartResult = flowchart(normalizedCoordinates);
    std::cout << "-> " << pointGroupString(flowchartResult.first) << " with certainty " << flowchartResult.second << "\n";

    BOOST_CHECK_MESSAGE(
      expectedPointGroup == flowchartResult.first,
      "Expected " << pointGroupString(expectedPointGroup)
      << " for symmetry " << name(shape)
      << ", got " << pointGroupString(flowchartResult.first)
      << " instead."
    );
  }
}
