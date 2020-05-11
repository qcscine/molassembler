/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Detail/Cartesian.h"
#include "molassembler/Temple/Random.h"
#include "molassembler/Options.h"

#include <Eigen/Geometry>

#include <cmath>

BOOST_AUTO_TEST_CASE(dihedralTests) {
  using namespace Scine;
  using namespace Molassembler;

  Scine::Utils::PositionCollection positions(4, 3);
  positions.row(0) = Scine::Utils::Position {1, 0, 0};
  positions.row(1) = Scine::Utils::Position {0, 0, 0};
  positions.row(2) = Scine::Utils::Position {0, 1, 0};

  const Eigen::Vector3d lastPosition {1, 1, 0};

  positions.row(3) = lastPosition;

  for(
    const double randomAngle :
    Temple::Random::getN<double>(
      -M_PI + 0.01,
      M_PI - 0.01,
      100,
      Scine::Molassembler::randomnessEngine()
    )
  ) {

    positions.row(3) = Eigen::AngleAxisd(
      randomAngle,
      Eigen::Vector3d::UnitY()
    ) * lastPosition;

    const double forwardDihedral = cartesian::dihedral(
      positions.row(0),
      positions.row(1),
      positions.row(2),
      positions.row(3)
    );

    BOOST_CHECK_MESSAGE(
       std::fabs(forwardDihedral - randomAngle) < 1e-10,
      "Twist angle: " << randomAngle << ", reported angle: " << forwardDihedral
    );

    const double reverseDihedral = cartesian::dihedral(
      positions.row(3),
      positions.row(2),
      positions.row(1),
      positions.row(0)
    );

    BOOST_CHECK_MESSAGE(
      std::fabs(forwardDihedral - reverseDihedral) < 1e-10,
      "Dihedral of reverse sequence is not identical: " << forwardDihedral
        << ", " << reverseDihedral
    );
  }
}
