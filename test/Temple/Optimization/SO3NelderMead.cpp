/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Temple/Optimization/SO3NelderMead.h"

#include <iostream>

using namespace Scine;

BOOST_AUTO_TEST_CASE(SO3NelderMead) {
  struct EigenValueDecomposition {
    static inline double square(double x) noexcept {
      return x * x;
    }

    double operator() (const Eigen::Matrix3d& rotation) const {
      Eigen::Matrix3d X;
      X <<  5,  2,  1,
            2,  7,  3,
            1,  3, 10;
      Eigen::Matrix3d partiallyDiagonal = rotation * X * rotation.transpose();
      double offDiagonalSquares = 0;
      for(unsigned i = 0; i < 2; ++i) {
        for(unsigned j = i + 1; j < 3; ++j) {
          offDiagonalSquares += square(partiallyDiagonal(i, j)) + square(partiallyDiagonal(j, i));
        }
      }
      return offDiagonalSquares;
    }
  };

  struct NelderMeadChecker {
    bool shouldContinue(unsigned iteration, double lowestValue, double stddev) {
      //std::cout << "Iteration " << iteration << ": " << lowestValue << " +- " << stddev << "\n";
      return iteration < 1000 && lowestValue >= 1e-5 && stddev > 1e-7;
    }
  };

  using OptimizerType = Temple::SO3NelderMead<>;

  { /* Squared distance commutative */
    const auto A = OptimizerType::Manifold::randomRotation();
    const auto B = OptimizerType::Manifold::randomRotation();
    BOOST_CHECK_MESSAGE(
      std::fabs(OptimizerType::Manifold::distanceSquared(A, B) - OptimizerType::Manifold::distanceSquared(B, A)) < 1e-10,
      "Squared distance calculation is not commutative"
    );
  }

  { /* Geodesic interpolation works */
    const Eigen::Matrix3d zPi = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const Eigen::Matrix3d zPiHalf = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    BOOST_CHECK_MESSAGE(
      OptimizerType::Manifold::geodesic(zPi, zPiHalf, 0).isApprox(zPiHalf, 1e-4),
      "geodesicExtrapolation(a, b, 0) != b"
    );
    BOOST_CHECK_MESSAGE(
      OptimizerType::Manifold::geodesic(zPi, zPiHalf, 1).isApprox(zPi, 1e-4),
      "geodesicExtrapolation(a, b, 1) != a"
    );

    const Eigen::AngleAxisd interpolation {OptimizerType::Manifold::geodesic(zPiHalf, zPi, 0.5)};
    BOOST_CHECK_MESSAGE(
      interpolation.axis().isApprox(Eigen::Vector3d::UnitZ(), 1e-4),
      "Axis of interpolated rotation is no longer z, but: " << interpolation.axis().transpose()
    );
    BOOST_CHECK_MESSAGE(
      std::fabs(interpolation.angle() - 3 * M_PI / 4) < 1e-4,
      "Angle of interpolated rotation is not 3 pi / 4 (" << (3 * M_PI / 4) << "), but " << interpolation.angle()
    );
  }

  bool pass = false;
  OptimizerType::OptimizationReturnType result;
  const unsigned maxAttempts = 10;
  for(unsigned i = 0; i < maxAttempts; ++i) {
    auto simplexVertices = OptimizerType::randomParameters();
    result = Temple::SO3NelderMead<>::minimize(
      simplexVertices,
      EigenValueDecomposition {},
      NelderMeadChecker {}
    );
    if(std::fabs(result.value) < 1e-2) {
      pass = true;
      break;
    } else {
      std::cout << "Attempt " << i << " lowest SO(3) minim. value = " << result.value << "\n";
    }
  }

  BOOST_CHECK_MESSAGE(
    pass,
    "SO(3) Nelder-Mead does not find minimization of EigenValueDecomposition problem in three attempts, value is "
    << result.value << " after " << result.iterations
    << " iterations."
  );
}
