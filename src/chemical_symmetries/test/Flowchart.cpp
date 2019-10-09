/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "chemical_symmetries/Flowchart.h"
#include "chemical_symmetries/Symmetries.h"

#include <iostream>

using namespace Scine;
using namespace Symmetry;

extern const std::string& pointGroupString(PointGroup group);
extern PositionCollection addOrigin(const PositionCollection& vs);

BOOST_AUTO_TEST_CASE(Flowcharting) {
  for(const Name symmetryName : allNames) {
    if(size(symmetryName) >= 4) {
      break;
    }

    auto normalizedCoordinates = detail::normalize(
      addOrigin(symmetryData().at(symmetryName).coordinates)
    );
    standardizeTop(normalizedCoordinates);

    std::cout << "Flowcharting " << name(symmetryName) << "\n";
    const PointGroup expectedPointGroup = pointGroup(symmetryName);
    const auto flowchartResult = flowchart(normalizedCoordinates);
    std::cout << "-> " << pointGroupString(flowchartResult.first) << " with certainty " << flowchartResult.second << "\n";

    BOOST_CHECK_MESSAGE(
      expectedPointGroup == flowchartResult.first,
      "Expected " << pointGroupString(expectedPointGroup)
      << " for symmetry " << name(symmetryName)
      << ", got " << pointGroupString(flowchartResult.first)
      << " instead."
    );
  }
}
