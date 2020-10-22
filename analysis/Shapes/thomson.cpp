/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Shapes/ContinuousMeasures.h"

#include "Molassembler/Temple/Optimization/Lbfgs.h"
#include "Molassembler/Temple/Adaptors/Iota.h"
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Functional.h"

#include <Eigen/Core>
#include <iostream>

using namespace Scine;
using namespace Molassembler;
using namespace Shapes;

using M = Eigen::Matrix<double, 3, Eigen::Dynamic>;

M cappedOctahedron() {
  M m(3, 7);
  m.col(0) = Eigen::Vector3d { 0.000000,  0.000000,  1.000000};
  m.col(1) = Eigen::Vector3d { 0.956305,  0.000000,  0.292372};
  m.col(2) = Eigen::Vector3d {-0.478152,  0.828184,  0.292372};
  m.col(3) = Eigen::Vector3d {-0.478152, -0.828184,  0.292372};
  m.col(4) = Eigen::Vector3d { 0.400888,  0.694358, -0.597625};
  m.col(5) = Eigen::Vector3d {-0.801776,  0.000000, -0.597625};
  m.col(6) = Eigen::Vector3d { 0.400888, -0.694358, -0.597625};
  return m;
}

M cappedTrigonalPrism() {
  M m(3, 7);
  m.col(0) = Eigen::Vector3d { 0.000000,  0.000000,  1.000000};
  m.col(1) = Eigen::Vector3d { 0.990268,  0.000000,  0.139173};
  m.col(2) = Eigen::Vector3d { 0.000000,  0.990268,  0.139173};
  m.col(3) = Eigen::Vector3d {-0.990268,  0.000000,  0.139173};
  m.col(4) = Eigen::Vector3d {-0.000000, -0.990268,  0.139173};
  m.col(5) = Eigen::Vector3d { 0.414628,  0.414628, -0.810042};
  m.col(6) = Eigen::Vector3d {-0.414628, -0.414628, -0.810042};
  return m;
}

M bicappedTrigonalPrism() {
  M m(3, 8);

  const double C0 = 0.288675134594812882254574390251;
  const double C1 = 0.497890957890680203327709376178;
  const double C2 = 0.577350269189625764509148780502;
  const double C3 = 0.862372435695794524549321018676;

  m.col(0) = Eigen::Vector3d { 0.5,  0.5,  C0};
  m.col(1) = Eigen::Vector3d { 0.5, -0.5,  C0};
  m.col(2) = Eigen::Vector3d {-0.5,  0.5,  C0};
  m.col(3) = Eigen::Vector3d {-0.5, -0.5,  C0};
  m.col(4) = Eigen::Vector3d { 0.0,  0.5, -C2};
  m.col(5) = Eigen::Vector3d { 0.0, -0.5, -C2};
  m.col(6) = Eigen::Vector3d {  C3,  0.0, -C1};
  m.col(7) = Eigen::Vector3d { -C3,  0.0, -C1};

  return m;
}

M tricappedTrigonalPrism() {
  M m(3, 9);

  const double C0 = 0.288675134594812882254574390251;
  const double C1 = 0.497890957890680203327709376178;
  const double C2 = 0.577350269189625764509148780502;
  const double C3 = 0.862372435695794524549321018676;
  const double C4 = 0.995781915781360406655418752356;

  m.col(0) = Eigen::Vector3d { 0.5, -C0,  0.5};
  m.col(1) = Eigen::Vector3d { 0.5, -C0, -0.5};
  m.col(2) = Eigen::Vector3d {-0.5, -C0,  0.5};
  m.col(3) = Eigen::Vector3d {-0.5, -C0, -0.5};
  m.col(4) = Eigen::Vector3d { 0.0,  C2,  0.5};
  m.col(5) = Eigen::Vector3d { 0.0,  C2, -0.5};
  m.col(6) = Eigen::Vector3d {  C3,  C1,  0.0};
  m.col(7) = Eigen::Vector3d { -C3,  C1,  0.0};
  m.col(8) = Eigen::Vector3d { 0.0, -C4,  0.0};

  return m;
}

M cappedSquareAntiprism() {
  M m(3, 9);

  const double C0 = 0.42044820762685727151556273811;
  const double C1 = 0.70710678118654752440084436210;
  const double C2 = 1.12755498881340479591640710022;

  m.col(0) = Eigen::Vector3d { 0.0,   C1,  C0};
  m.col(1) = Eigen::Vector3d { 0.0,  -C1,  C0};
  m.col(2) = Eigen::Vector3d {  C1,  0.0,  C0};
  m.col(3) = Eigen::Vector3d { -C1,  0.0,  C0};
  m.col(4) = Eigen::Vector3d { 0.5,  0.5, -C0};
  m.col(5) = Eigen::Vector3d { 0.5, -0.5, -C0};
  m.col(6) = Eigen::Vector3d {-0.5,  0.5, -C0};
  m.col(7) = Eigen::Vector3d {-0.5, -0.5, -C0};
  m.col(8) = Eigen::Vector3d { 0.0,  0.0,  C2};

  return m;
}

M bicappedSquareAntiprism() {
  M m(3, 10);
  const double C0 = 0.42044820762685727151556273811;
  const double C1 = 0.70710678118654752440084436210;
  const double C2 = 1.12755498881340479591640710022;

  m.col(0) = Eigen::Vector3d {  C1,  0.0,  C0};
  m.col(1) = Eigen::Vector3d { -C1,  0.0,  C0};
  m.col(2) = Eigen::Vector3d { 0.0,   C1,  C0};
  m.col(3) = Eigen::Vector3d { 0.0,  -C1,  C0};
  m.col(4) = Eigen::Vector3d { 0.5,  0.5, -C0};
  m.col(5) = Eigen::Vector3d { 0.5, -0.5, -C0};
  m.col(6) = Eigen::Vector3d {-0.5,  0.5, -C0};
  m.col(7) = Eigen::Vector3d {-0.5, -0.5, -C0};
  m.col(8) = Eigen::Vector3d { 0.0,  0.0,  C2};
  m.col(9) = Eigen::Vector3d { 0.0,  0.0, -C2};

  return m;
}

M trigonalDodecahedron() {
  M m(3, 8);

  const double C0 = 0.205561565853259538077851199951;
  const double C1 = 0.644584273224154984541338729084;
  const double C2 = 0.783930924232563649920229273668;

  m.col(0) = Eigen::Vector3d { 0.5,  0.0, -C2};
  m.col(1) = Eigen::Vector3d {-0.5,  0.0, -C2};
  m.col(2) = Eigen::Vector3d { 0.0,  0.5,  C2};
  m.col(3) = Eigen::Vector3d { 0.0, -0.5,  C2};
  m.col(4) = Eigen::Vector3d {  C1,  0.0,  C0};
  m.col(5) = Eigen::Vector3d { -C1,  0.0,  C0};
  m.col(6) = Eigen::Vector3d { 0.0,   C1, -C0};
  m.col(7) = Eigen::Vector3d { 0.0,  -C1, -C0};

  return m;
}

using F = std::function<M()>;

Eigen::VectorXd transformToAngles(M m) {
  const unsigned N = m.cols();
  Eigen::VectorXd angles(2 * N);
  for(unsigned i = 0; i < N; ++i) {
    assert(std::fabs(m.col(i).norm() - 1) < 1e-10);
    angles(2 * i) = std::acos(m.col(i).z());
    angles(2 * i + 1) = std::atan2(m.col(i).y(), m.col(i).x());
  }

  assert(
    Temple::all_of(
      Temple::Adaptors::allPairs(Temple::Adaptors::range(N)),
      [&](const unsigned i, const unsigned j) -> bool {
        return !angles.template segment<2>(2 * i).isApprox(
          angles.template segment<2>(2 * j),
          1e-10
        );
      }
    )
  );

  return angles;
}

M transformToPositions(const Eigen::VectorXd& angles) {
  const unsigned N = angles.size() / 2;
  M m(3, N);
  for(unsigned i = 0; i < N; ++i) {
    m.col(i) = Eigen::Vector3d {
      std::sin(angles(2 * i)) * std::cos(angles(2 * i + 1)),
      std::sin(angles(2 * i)) * std::sin(angles(2 * i + 1)),
      std::cos(angles(2 * i))
    };
  }
  return m;
}

struct AnglePotential {
  static Eigen::Vector3d makeVector(unsigned i, const Eigen::VectorXd& angles) {
    const double thetaI = angles(2 * i);
    const double phiI = angles(2 * i + 1);
    return {
      std::sin(thetaI) * std::cos(phiI),
      std::sin(thetaI) * std::sin(phiI),
      std::cos(thetaI)
    };
  }

  static double inverseDistanceSum(const Eigen::VectorXd& angles) {
    const unsigned N = angles.size() / 2;
    double sum = 0.0;
    for(unsigned i = 0; i < N; ++i) {
      const Eigen::Vector3d ri = makeVector(i, angles);
      for(unsigned j = i + 1; j < N; ++j) {
        const Eigen::Vector3d rj = makeVector(j, angles);
        sum += 1.0 / (ri - rj).norm();
      }
    }
    return sum;
  }

  void operator() (
    const Eigen::VectorXd& angles,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient
  ) {
    value = inverseDistanceSum(angles);
    gradient = Temple::Optimization::numericalGradient(&inverseDistanceSum, angles);
  }
};

void writeXYZ(const std::string& name, const M& m) {
  const unsigned N = m.cols();
  std::ofstream xyz(name);
  xyz << (N + 1) << "\n\n";

  xyz << std::fixed;
  for(unsigned i = 0; i < N; ++i) {
    xyz << "H  " << m.col(i).x() << " " << m.col(i).y() << " " << m.col(i).z() << "\n";
  }

  xyz << "Fe " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
}

void spherize(Shape shape, const M& p) {
  M positions = p;
  positions.colwise().normalize();

  using Optimizer = Temple::Lbfgs<>;
  auto parameters = transformToAngles(positions);

  std::cout << "positions:\n" << positions.transpose() << "\n";
  std::cout << "parameters: " << parameters.transpose() << "\n";
  std::cout << name(shape) << " f = " << AnglePotential::inverseDistanceSum(parameters) << "\n";

  struct Checker {
    bool shouldContinue(const unsigned iterations, const Optimizer::StepValues& step) {
      static unsigned counter = 0;
      writeXYZ("iter_" + std::to_string(counter) + ".xyz", transformToPositions(step.parameters.current));
      ++counter;
      const double gradNorm = step.gradients.current.norm();
      std::cout << "norm grad f = " << gradNorm << "\n";
      return iterations < 1000 && gradNorm > 1e-4;
    }
  };

  Optimizer optimizer;
  auto result = optimizer.minimize(
    parameters,
    AnglePotential {},
    Checker {}
  );

  auto post = transformToPositions(parameters);

  std::cout << " f = " << result.value << ", norm grad f = " << result.gradient.norm() << "\n";
  std::cout << post.transpose() << "\n";
}

int main() {
  spherize(Shape::CappedOctahedron, cappedOctahedron());
  // spherize(Shape::CappedTrigonalPrism, cappedTrigonalPrism());
  // spherize(Shape::BicappedTrigonalPrism, bicappedTrigonalPrism()},
  // spherize(Shape::CappedSquareAntiprism, cappedSquareAntiprism());
  // spherize(Shape::TrigonalDodecahedron, trigonalDodecahedron());
}
