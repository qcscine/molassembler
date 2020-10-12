/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details. for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Shapes/CoordinateSystemTransformation.h"

BOOST_AUTO_TEST_CASE(CoordinateSystemTransformation, *boost::unit_test::label("Shapes")) {
  using namespace Scine::Molassembler::Shapes;
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
