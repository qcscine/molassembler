/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "shapes/CoordinateSystemTransformation.h"

using namespace Scine;
using namespace Symmetry;

BOOST_AUTO_TEST_CASE(CoordinateSystemTransformation) {
  for(unsigned i = 0; i < 20; ++i) {
    // Create two random coordinate systems
    const CoordinateSystem a = CoordinateSystem::random();
    const CoordinateSystem b = CoordinateSystem::random();

    const auto rot = rotationMatrix(a, b);

    BOOST_CHECK((rot * a.x).isApprox(b.x, 1e-10));
    BOOST_CHECK((rot * a.y).isApprox(b.y, 1e-10));
    BOOST_CHECK((rot * a.z).isApprox(b.z, 1e-10));
  }
}
