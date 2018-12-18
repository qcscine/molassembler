/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Detail/DelibHelpers.h"
#include "temple/Random.h"
#include "molassembler/Options.h"

#include <Eigen/Geometry>

#include <cmath>

BOOST_AUTO_TEST_CASE(dihedralTests) {
  using namespace Scine::molassembler::DelibHelpers;

  Scine::Utils::PositionCollection positions(4, 3);
  positions.row(0) = Scine::Utils::Position {1, 0, 0};
  positions.row(1) = Scine::Utils::Position {0, 0, 0};
  positions.row(2) = Scine::Utils::Position {0, 1, 0};

  const Eigen::Vector3d lastPosition {1, 1, 0};

  positions.row(3) = lastPosition;

  for(
    const double randomAngle :
    temple::random::getN<double>(
      -M_PI + 0.01,
      M_PI - 0.01,
      100,
      Scine::molassembler::randomnessEngine()
    )
  ) {

    positions.row(3) = Eigen::AngleAxisd(
      randomAngle,
      Eigen::Vector3d::UnitY()
    ) * lastPosition;

    const double dihedral = getDihedral(positions, 0, 1, 2, 3);

    BOOST_CHECK_MESSAGE(
       std::fabs(dihedral - randomAngle) < 1e-10,
      "Twist angle: " << randomAngle << ", reported angle: " << dihedral
    );

    const double reverseDihedral = getDihedral(positions, 3, 2, 1, 0);

    BOOST_CHECK_MESSAGE(
      std::fabs(dihedral - reverseDihedral) < 1e-10,
      "Dihedarl of reverse sequence is not identical: " << dihedral <<", " << reverseDihedral
    );
  }
}
