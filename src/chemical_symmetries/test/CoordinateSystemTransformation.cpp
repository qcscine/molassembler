/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "chemical_symmetries/CoordinateSystemTransformation.h"

using namespace Scine;
using namespace Symmetry;

BOOST_AUTO_TEST_CASE(CoordinateSystemTransformation) {
  for(unsigned i = 0; i < 20; ++i) {
    // Create two random coordinate systems
    const Eigen::Vector3d c1 = Eigen::Vector3d::Random().normalized();
    const Eigen::Vector3d c2 = Eigen::Vector3d::Random().cross(c1).normalized();
    const CoordinateSystem a {c1, c2};

    // Default xyz system
    const Eigen::Vector3d c3 = Eigen::Vector3d::Random().normalized();
    const Eigen::Vector3d c4 = Eigen::Vector3d::Random().cross(c3).normalized();
    const CoordinateSystem b {c3, c4};

    const auto rot = rotationMatrix(a, b);

    BOOST_CHECK((rot * a.x).isApprox(b.x, 1e-10));
    BOOST_CHECK((rot * a.y).isApprox(b.y, 1e-10));
    BOOST_CHECK((rot * a.z).isApprox(b.z, 1e-10));
  }
}
