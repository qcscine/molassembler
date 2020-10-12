/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "boost/test/unit_test.hpp"

#include "Molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "Molassembler/DistanceGeometry/TetrangleSmoothing.h"

using namespace Scine::Molassembler;
using namespace DistanceGeometry;



BOOST_AUTO_TEST_CASE(TriangleSmoothingFloydExplicit, *boost::unit_test::label("DG")) {
  Eigen::Matrix4d input;
  input <<   0.0,   1.0, 100.0,   1.0,
             1.0,   0.0,   1.0, 100.0,
             0.5,   1.0,   0.0,   1.0,
             1.0,   0.5,   1.0,   0.0;

  Eigen::Matrix4d expected;
  expected << 0.0, 1.0, 2.0, 1.0,
              1.0, 0.0, 1.0, 2.0,
              0.5, 1.0, 0.0, 1.0,
              1.0, 0.5, 1.0, 0.0;

  Eigen::MatrixXd a = input;
  DistanceBoundsMatrix::smooth(a);
  BOOST_CHECK_MESSAGE(
    a.isApprox(expected, 1e-10),
    "Triangle smoothing of example matrix does not give triangle inequality limits: Expected\n"
    << expected << "\nGot:\n" << a << "\n"
  );
}

BOOST_AUTO_TEST_CASE(TetrangleSmoothingExplicit, *boost::unit_test::label("DG")) {
  Eigen::Matrix4d input;
  input <<   0.0,   1.0, 100.0,   1.0,
             1.0,   0.0,   1.0, 100.0,
             0.5,   1.0,   0.0,   1.0,
             1.0,   0.5,   1.0,   0.0;

  Eigen::Matrix4d expected;
  expected << 0.0, 1.0, 2.0, 1.0,
              1.0, 0.0, 1.0, 2.0,
              0.5, 1.0, 0.0, 1.0,
              1.0, 0.5, 1.0, 0.0;

  Eigen::Matrix4d a = input;
  const unsigned iterations = tetrangleSmooth(a);

  BOOST_CHECK_MESSAGE(
    (a.array() <= expected.array()).all(),
    "Tetrangle smoothing (" << iterations << " iterations) of example matrix "
    << "does not give values smaller than the inequality limits: Expected at most\n"
    << expected << "\nGot:\n" << a << "\n"
  );
}

BOOST_AUTO_TEST_CASE(TriangleSmoothingDetectsViolations, *boost::unit_test::label("DG")) {
  Eigen::Matrix3d impossibleBounds;
  impossibleBounds << 0.0, 1.0, 4.0,
                      1.0, 0.0, 2.0,
                      4.0, 2.0, 0.0;

  BOOST_CHECK_THROW(DistanceBoundsMatrix::smooth(impossibleBounds), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TetrangleSmoothingDetectsViolations, *boost::unit_test::label("DG")) {
  Eigen::Matrix3d impossibleBounds;
  impossibleBounds << 0.0, 1.0, 4.0,
                      1.0, 0.0, 2.0,
                      4.0, 2.0, 0.0;

  BOOST_CHECK_THROW(tetrangleSmooth(impossibleBounds), std::runtime_error);
}
