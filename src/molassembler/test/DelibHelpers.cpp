#define BOOST_TEST_MODULE DelibHelperTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "detail/DelibHelpers.h"
#include "temple/Random.h"

#include <Eigen/Geometry>

BOOST_AUTO_TEST_CASE(dihedralTests) {
  using namespace molassembler::DelibHelpers;

  temple::Generator rng;

  Delib::PositionCollection positions;
  positions.push_back(
    Delib::Position {
      Eigen::Vector3d {
        1, 0, 0
      }
    }
  );
  positions.push_back(
    Delib::Position {
      Eigen::Vector3d {
        0, 0, 0
      }
    }
  );
  positions.push_back(
    Delib::Position {
      Eigen::Vector3d {
        0, 1, 0
      }
    }
  );

  const Eigen::Vector3d lastPosition {1, 1, 0};

  positions.push_back(
    Delib::Position {lastPosition}
  );

  for(const double& randomAngle : rng.getN<double>(-M_PI + 0.01, M_PI - 0.01, 100)) {

    positions[3] = Delib::Position {
      Eigen::AngleAxisd(
        randomAngle,
        Eigen::Vector3d::UnitY()
      ) * lastPosition
    };

    const auto dihedral = getDihedral(positions, 0, 1, 2, 3);

    BOOST_CHECK_MESSAGE(
       std::fabs(dihedral - randomAngle) < 1e-10,
      "Twist angle: " << randomAngle << ", reported angle: " << dihedral
    );
  }
}
