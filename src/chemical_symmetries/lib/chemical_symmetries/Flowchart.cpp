#include "chemical_symmetries/Flowchart.h"

#include "temple/constexpr/JSF.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Adaptors/Zip.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Stringify.h"

#include <iostream>

#include <Eigen/Geometry>

/* TODO
 * - Use statistics better (replace use of average and certainty calculation
 *   with a statistical test over the probability distribution function). The
 *   best-fitting distribution is the Weibull. Using random point clouds and
 *   fitting the distribution, you can use the cumulative density function to
 *   argue that a particular csm value is not the result of a random point
 *   cloud, but more structured.
 * - Need to open up recognition functions to allow restrictions of the search
 *   space or passing explicit rotation matrices to test when setting up the
 *   simplex.
 */

namespace Scine{
namespace Symmetry {
namespace distributions {

struct Weibull {
  double shape;
  double scale;

  double pdf(const double x) const {
    if(x < 0) {
      return 0.0;
    }
    return (shape / scale) * std::pow(x / scale, shape - 1) * std::exp(-std::pow(x / scale, shape));
  }

  double cdf(const double x) const {
    if(x < 0) {
      return 0.0;
    }

    return 1.0 - std::exp(-std::pow(x / scale, shape));
  }
};

namespace parameters {

/* Weibull distribution parameters fitted to 1000 inversion CSM samples for
 * N = 2 to N = 8 points (N does not count the origin). Points have uniform
 * direction and radius ~ norm(mu = 1.0, sigma = 0.2).
 */
static constexpr std::array<Weibull, 7> inversion {{
  {1.37266101068785, 19.6787332218856},
  {1.58376472686282, 9.02421976122357},
  {2.42079642881566, 10.9152762243839},
  {2.75217499779139, 8.02096277130682},
  {3.41495571561817, 8.91666122941434},
  {3.30943108907646, 7.38795307826982},
  {4.01346300739247, 7.88065600872248}
}};

/* Weibull distribution parameters fitted to 1000 Cinf CSM samples for
 * N = 2 to N = 8 points (N does not count the origin). Points have uniform
 * direction and radius ~ norm(mu = 1.0, sigma = 0.2).
 */
static constexpr std::array<Weibull, 7> Cinf {{
  {1.16820081756571, 16.5667847270495},
  {1.96060538456455, 19.6203168014297},
  {2.65552778791327, 21.3725178846392},
  {3.21615985344926, 22.9719471115388},
  {3.74123238890814, 24.4653029244794},
  {4.05645938407287, 25.9935790116030},
  {4.36929088095962, 27.1728431327878}
}};

/* Weibull distribution parameters fitted to 1000 sigma CSM samples for
 * N = 3 to N = 8 points (N does not count the origin). Points have uniform
 * direction and radius ~ norm(mu = 1.0, sigma = 0.2).
 *
 * N = 2 is omitted since any two points and the origin can be encompassed in
 * a symmetry plane with zero CSM.
 */
static constexpr std::array<Weibull, 6> sigma {{
  {0.749932612735404, 0.340017576942346},
  {1.197579311240620, 0.894300802240169},
  {1.601645462693400, 1.193257820437570},
  {1.873029164899980, 1.520395558392160},
  {2.313458407275550, 1.851023412032510},
  {2.580016853674300, 2.110094243561410}
}};

} // namespace parameters

Weibull inversion(const unsigned N) {
  assert(N >= 2);
  return parameters::inversion.at(N - 2);
}

Weibull Cinf(const unsigned N) {
  assert(N >= 2);
  return parameters::Cinf.at(N - 2);
}

Weibull sigma(const unsigned N) {
  assert(N != 2 && "Do not try a distribution for N = 2, CSM is always zero!");
  return parameters::sigma.at(N - 3);
}

} // namespace distributions

template<typename PRNG>
Eigen::Vector3d normallyDistributedVectorInSphere(const double radius, PRNG& prng) {
  std::normal_distribution<double> radiusDistribution {1.0, 0.2};
  std::normal_distribution<double> normal {};
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << normal(prng),
         normal(prng),
         normal(prng);
  }

  return radius * (1 + radiusDistribution(prng)) * v / v.norm();
}

template<typename PRNG>
Eigen::Matrix<double, 3, Eigen::Dynamic> generateRandomCoordinates(unsigned P, PRNG& prng) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> positions(3, P + 1);
  positions.col(0) = Eigen::Vector3d::Zero();
  for(unsigned i = 1; i <= P; ++i) {
    positions.col(i) = normallyDistributedVectorInSphere(1.0, prng);
  }
  return positions;
}

template<typename PRNG, typename F>
double averageRandomCsm(const unsigned N, PRNG& prng, F&& f) {
  assert(N > 1);
  const unsigned nExperiments = 100;
  return temple::average(
    temple::map(
      temple::adaptors::range(nExperiments),
      [&](unsigned /* i */) -> double {
        auto normalized = detail::normalize(generateRandomCoordinates(N - 1, prng));
        Top top = standardizeTop(normalized);
        if(top == Top::Asymmetric) {
          reorientAsymmetricTop(normalized);
        }
        return f(normalized);
      }
    )
  );
}

double normalizedAngle(
  const Eigen::Vector3d& a,
  const Eigen::Vector3d& b
) {
  assert(std::fabs(a.norm() - 1) < 1e-5 && "a is not normalized");
  assert(std::fabs(b.norm() - 1) < 1e-5 && "b is not normalized");

  return std::acos(a.dot(b));
}

bool angleWithinBounds(
  const Eigen::Vector3d& a,
  const Eigen::Vector3d& b,
  const double expected,
  const double tolerance
) {
  const double angle = normalizedAngle(a, b);
  return (expected - tolerance <= angle && angle <= expected + tolerance);
}

Eigen::Vector3d perpendicular(const Eigen::Vector3d& v) {
  return v.cross(Eigen::Vector3d::Random()).normalized();
}

namespace predicates {

double inversion(const PositionCollection& positions, temple::jsf::JSF64& prng) {
  const unsigned N = positions.cols();
  const double inversionCsm = csm::element(positions, elements::Inversion {});
  const double randomInversionCsm = averageRandomCsm(N, prng,
    [](const PositionCollection& positions) -> double {
      return csm::element(positions, elements::Inversion {});
    }
  );
  const double iCertainty = std::max(randomInversionCsm - inversionCsm, 0.0) / randomInversionCsm;
  std::cout << "CSM(i) = " << inversionCsm << ", random CSM for " << N << " points: " << randomInversionCsm << ". Certainty = " << iCertainty << "\n";
  return iCertainty;
}

struct SigmaResult {
  double certainty;
  Eigen::Vector3d normal;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
SigmaResult sigma(
  const PositionCollection& positions,
  const Eigen::Vector3d& normal,
  temple::jsf::JSF64& prng
) {
  const unsigned N = positions.cols();
  const auto sigmaHResult = csm::optimize(
    positions,
    elements::Reflection {normal}
  );

  const double randomSigmaCsm = averageRandomCsm(N, prng,
    [](const PositionCollection& positions) -> double {
      return csm::optimize(
        positions,
        elements::Reflection {Eigen::Vector3d::UnitZ()}
      ).first;
    }
  );
  const double sigmaHCertainty = std::max(randomSigmaCsm - sigmaHResult.first, 0.0) / randomSigmaCsm;

  return {
    sigmaHCertainty,
    sigmaHResult.second.normal
  };
}

} // namespace predicates

PointGroup computePointGroup(
  const PointGroup base,
  const unsigned mainAxisOrder
) {
  assert(mainAxisOrder >= 2);
  return static_cast<PointGroup>(
    static_cast<std::underlying_type<PointGroup>::type>(base) + (mainAxisOrder - 2)
  );
}

std::pair<PointGroup, double> linearPointGroups(
  const PositionCollection& normalizedPositions,
  temple::jsf::JSF64& prng,
  const double cumulativeCertainty
) {
  // Yes branch for Cinf -> Linear molecules group
  // Test for inversion
  const double iCertainty = predicates::inversion(normalizedPositions, prng);
  if(iCertainty >= 0.5) {
    return {PointGroup::Dinfh, cumulativeCertainty * iCertainty};
  }

  return {PointGroup::Cinfv, cumulativeCertainty * (1 - iCertainty)};
}

std::pair<PointGroup, double> smallPointGroups(
  const PositionCollection& normalizedPositions,
  temple::jsf::JSF64& prng,
  double cumulativeCertainty
) {
  // No Cn axis found -> Small symmetries group
  // Test for reflection plane
  auto sigma = predicates::sigma(normalizedPositions, Eigen::Vector3d::UnitZ(), prng);

  if(sigma.certainty >= 0.5) {
    return {PointGroup::Cs, cumulativeCertainty * sigma.certainty};
  }

  /* No branch for reflection */
  cumulativeCertainty *= (1 - sigma.certainty);

  // Test for inversion center
  const double iCertainty = predicates::inversion(normalizedPositions, prng);
  if(iCertainty >= 0.5) {
    return {PointGroup::Ci, cumulativeCertainty * iCertainty};
  }

  return {PointGroup::C1, cumulativeCertainty * (1 - iCertainty)};
}

std::pair<PointGroup, double> mediumPointGroups(
  const PositionCollection& normalizedPositions,
  const elements::Rotation& mainAxis,
  const unsigned N,
  temple::jsf::JSF64& prng,
  double cumulativeCertainty
) {
  // Test for n C2s orthogonal to the main axis
  const Eigen::Vector3d axisPerpendicular = perpendicular(mainAxis.axis);
  const auto firstC2Result = csm::optimize(
    normalizedPositions,
    elements::Rotation::Cn(axisPerpendicular, 2)
  );
  const double randomC2Csm = averageRandomCsm(N, prng,
    [](const PositionCollection& positions) -> double {
      return csm::optimize(
        positions,
        elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), 2)
      ).first;
    }
  );
  const double firstC2Certainty = std::max(randomC2Csm - firstC2Result.first, 0.0) / randomC2Csm;

  if(firstC2Certainty >= 0.5 && angleWithinBounds(firstC2Result.second.axis, mainAxis.axis, M_PI / 2, M_PI / 6)) {
    // Yes branch for n C2 perp to Cn
    cumulativeCertainty *= firstC2Certainty;
    const auto sigmaH = predicates::sigma(normalizedPositions, mainAxis.axis, prng);

    if(sigmaH.certainty >= 0.5 && normalizedAngle(sigmaH.normal, mainAxis.axis) <= M_PI / 6) {
      return {
        computePointGroup(PointGroup::D2h, mainAxis.n),
        cumulativeCertainty * sigmaH.certainty
      };
    }

    // No branch for sigma_h
    cumulativeCertainty *= (1 - sigmaH.certainty);

    // Look for sigma_v between C2s
    const Eigen::Vector3d sigmaVNormal = elements::Rotation::Cn(mainAxis.axis, 2 * mainAxis.n).matrix() * firstC2Result.second.axis;
    const auto sigmaV = predicates::sigma(normalizedPositions, sigmaVNormal, prng);

    if(sigmaV.certainty >= 0.5 && angleWithinBounds(sigmaV.normal, mainAxis.axis, M_PI / 2, M_PI / 6)) {
      // Yes branch for n sigma_vs
      return {
        computePointGroup(PointGroup::D2d, mainAxis.n),
        cumulativeCertainty * sigmaV.certainty
      };
    }

    // No branch for n sigma_v
    return {
      computePointGroup(PointGroup::D2, mainAxis.n),
      cumulativeCertainty * (1 - sigmaV.certainty)
    };
  }

  // No branch for n C2 orthogonal to Cn
  cumulativeCertainty *= (1 - firstC2Certainty);
  // Look for S2n
  if(2 <= mainAxis.n && mainAxis.n <= 4) {
    const auto S2nResult = csm::optimize(
      normalizedPositions,
      elements::Rotation::Sn(mainAxis.axis, 2 * mainAxis.n)
    );
    const double randomS2nCsm = averageRandomCsm(N, prng,
      [&mainAxis](const PositionCollection& positions) -> double {
        return csm::optimize(
          positions,
          elements::Rotation::Sn(Eigen::Vector3d::UnitZ(), 2 * mainAxis.n)
        ).first;
      }
    );
    const double S2nCertainty = std::max(randomS2nCsm - S2nResult.first, 0.0) / randomS2nCsm;

    if(S2nCertainty >= 0.5 && normalizedAngle(S2nResult.second.axis, mainAxis.axis) <= M_PI / 6) {
      // Yes branch for S2n
      return {
        static_cast<PointGroup>(
          static_cast<std::underlying_type<PointGroup>::type>(PointGroup::S4) + (mainAxis.n - 2)
        ),
        cumulativeCertainty * S2nCertainty
      };
    }

    // No branch for S2n
    cumulativeCertainty *= (1 - S2nCertainty);
  }

  // sigmaH?
  auto sigmaH = predicates::sigma(normalizedPositions, mainAxis.axis, prng);
  if(sigmaH.certainty >= 0.5 && normalizedAngle(sigmaH.normal, mainAxis.axis) <= M_PI / 6) {
    return {
      computePointGroup(PointGroup::C2h, N),
      cumulativeCertainty * sigmaH.certainty
    };
  }

  // No branch, look for n sigma_v
  cumulativeCertainty *= (1 - sigmaH.certainty);
  auto sigmaV = predicates::sigma(normalizedPositions, perpendicular(mainAxis.axis), prng);
  if(sigmaV.certainty >= 0.5 && angleWithinBounds(sigmaV.normal, mainAxis.axis, M_PI / 2, M_PI / 6)) {
    return {
      computePointGroup(PointGroup::C2v, N),
      cumulativeCertainty * sigmaV.certainty
    };
  }

  // No branch for n sigma_v
  return {
    computePointGroup(PointGroup::C2, N),
    cumulativeCertainty * (1 - sigmaV.certainty)
  };
}

std::pair<PointGroup, double> flowchart(const PositionCollection& normalizedPositions) {
  assert(
    temple::accumulate(
      temple::adaptors::range(normalizedPositions.cols()),
      0u,
      [&](unsigned carry, unsigned col) -> unsigned {
        if(normalizedPositions.col(col).isApprox(Eigen::Vector3d::Zero(), 1e-5)) {
          return carry + 1;
        }

        return carry;
      }
    ) < 2u
  );

  static temple::jsf::JSF64 prng;
  const unsigned N = normalizedPositions.cols();

  // For reproducible behavior
  prng.seed(498123941);

  /* At each point of the decision tree, try to quantify whether the presence
   * of a symmetry element is random or not between 0 and 1. If the measure is
   * > 0.5, continue with yes propagating that measure along.
   */
  const double cinfCsm = csm::optimizeCinf(normalizedPositions);
  const double randomCinfCsm = averageRandomCsm(N, prng, csm::optimizeCinf);
  const double cinfCertainty = std::max(randomCinfCsm - cinfCsm, 0.0) / randomCinfCsm;

  std::cout << "CSM(Cinf) = " << cinfCsm << ", random CSM for " << N << " points: " << randomCinfCsm << ". Certainty = " << cinfCertainty << "\n";

  double cumulativeCertainty = 1.0;

  if(cinfCertainty >= 0.5) {
    cumulativeCertainty *= cinfCertainty;
    return linearPointGroups(normalizedPositions, prng, cumulativeCertainty);
  }

  // No Cinf branch
  cumulativeCertainty *= (1 - cinfCertainty);

  // Find highest-order proper rotation axis
  double highestMainAxisCertainty = 0.0;
  boost::optional<
    std::pair<double, elements::Rotation>
  > mainAxisPairOption;
  for(unsigned order = 8; order >= 2; --order) {
    const double randomOrderCn = averageRandomCsm(N, prng,
      [order](const PositionCollection& positions) -> double {
        return csm::optimize(
          positions,
          elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), order)
        ).first;
      }
    );

    const auto optimizedAxes = temple::map(
      temple::adaptors::range(3),
      [&](const unsigned col) -> std::pair<double, elements::Rotation> {
        return csm::optimize(
          normalizedPositions,
          elements::Rotation::Cn(
            Eigen::Matrix3d::Identity().col(col),
            order
          )
        );
      }
    );

    const auto bestAxis = *std::min_element(
      std::begin(optimizedAxes),
      std::end(optimizedAxes),
      [](const auto& a, const auto& b) -> bool {
        return a.first < b.first;
      }
    );

    const double axisCertainty = std::max(randomOrderCn - bestAxis.first, 0.0) / randomOrderCn;

    std::cout << "Best C" << order << " has CSM = " << bestAxis.first << ", random is " << randomOrderCn << ", certainty = " << axisCertainty << "\n";

    highestMainAxisCertainty = std::max(highestMainAxisCertainty, axisCertainty);

    if(axisCertainty >= 0.5) {
      mainAxisPairOption = bestAxis;
      break;
    }
  }

  if(!mainAxisPairOption) {
    cumulativeCertainty *= (1 - highestMainAxisCertainty);
    return smallPointGroups(normalizedPositions, prng, cumulativeCertainty);
  }

  /* Yes branch for main axis */
  cumulativeCertainty *= highestMainAxisCertainty;
  const elements::Rotation& mainAxis = mainAxisPairOption->second;

  /* Testing for 6 C5 only makes sense if the order of the main axis is exactly 6 */
  if(mainAxis.n == 6) {
    // TODO reasoning about probabilities here should include the prior (there is already a C6 present)
    // TODO reasoning here maybe should cover the certainties of all C5 axes
    // Generate a C5 axis offset some 60° from the C6 axis and check
    const Eigen::Vector3d axisPerpendicular = perpendicular(mainAxis.axis);
    const auto hypotheticalC5 = elements::Rotation::Cn(axisPerpendicular, 5);
    const double randomC5Csm = averageRandomCsm(N, prng,
      [hypotheticalC5](const PositionCollection& positions) -> double {
        return csm::optimize(positions, hypotheticalC5).first;
      }
    );
    const auto C5Result = csm::optimize(normalizedPositions, hypotheticalC5);
    const double C5Certainty = std::max(randomC5Csm - C5Result.first, 0.0) / randomC5Csm;

    if(C5Certainty >= 0.5) {
      // Yes branch for 6 C5
      cumulativeCertainty *= C5Certainty;
      // Test for inversion
      const double iCertainty = predicates::inversion(normalizedPositions, prng);
      if(iCertainty >= 0.5) {
        return {PointGroup::Ih, cumulativeCertainty * iCertainty};
      }

      return {PointGroup::I, cumulativeCertainty * (1 - iCertainty)};
    } else {
      // No branch does influence the cumulative certainty
      cumulativeCertainty *= (1 - C5Certainty);
    }
  }

  /* Testing for 3 C4 only makes sense if the order of the main axis is exactly 4 */
  if(mainAxis.n == 4) {
    /* TODO Super unsure about this block. The reasoning is weird, the cutoffs are
     * weird, everything is weird
     */

    auto approximatelyOrthogonal = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b) -> bool {
      const double angle = std::acos(a.dot(b));
      constexpr double tolerance = M_PI / 6;
      return (M_PI - tolerance <= angle && angle <= M_PI + tolerance);
    };

    const Eigen::Vector3d axisPerpendicular = perpendicular(mainAxis.axis);
    const auto firstHypotheticalC4 = elements::Rotation::Cn(axisPerpendicular, 4);
    const double randomC4Csm = averageRandomCsm(N, prng,
      [firstHypotheticalC4](const PositionCollection& positions) -> double {
        return csm::optimize(positions, firstHypotheticalC4).first;
      }
    );
    const auto firstC4Result = csm::optimize(normalizedPositions, firstHypotheticalC4);
    const double firstC4Certainty = std::max(randomC4Csm - firstC4Result.first, 0.0) / randomC4Csm;

    if(firstC4Certainty >= 0.5 && approximatelyOrthogonal(firstC4Result.second.axis, mainAxis.axis)) {
      const auto secondHypotheticalC4 = elements::Rotation::Cn(
        Eigen::AngleAxisd(M_PI / 2, mainAxis.axis) * firstC4Result.second.axis,
        4
      );
      const auto secondC4Result = csm::optimize(normalizedPositions, secondHypotheticalC4);
      const double secondC4Certainty = std::max(randomC4Csm - secondC4Result.first, 0.0) / randomC4Csm;

      if(
        secondC4Certainty >= 0.5
        && approximatelyOrthogonal(secondC4Result.second.axis, mainAxis.axis)
        && approximatelyOrthogonal(secondC4Result.second.axis, firstC4Result.second.axis)
      ) {
        // Yes branch for 3 C4
        cumulativeCertainty *= firstC4Certainty * secondC4Certainty;

        const double iCertainty = predicates::inversion(normalizedPositions, prng);
        if(iCertainty >= 0.5) {
          return {PointGroup::Oh, cumulativeCertainty * iCertainty};
        }

        return {PointGroup::O, cumulativeCertainty * (1 - iCertainty)};
      }
    }

    // TODO What about secondC4Certainty?
    cumulativeCertainty *= (1 - firstC4Certainty);
  }

  /* Testing for 4 C3 only makes sense if the order of the main axis is exactly 3 */
  if(mainAxis.n == 3) {
    auto angleIsRoughlyRight = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b) -> bool {
      const double angle = std::acos(a.dot(b));
      constexpr double tolerance = M_PI / 6;
      constexpr double expected = 2 * M_PI / 3;
      return (expected - tolerance <= angle && angle <= expected + tolerance);
    };

    const Eigen::Vector3d axisPerpendicular = perpendicular(mainAxis.axis);
    auto firstHypotheticalC3 = elements::Rotation::Cn(
      Eigen::AngleAxisd(2 * M_PI / 3, axisPerpendicular).toRotationMatrix() * mainAxis.axis,
      3
    );
    const auto firstC3Result = csm::optimize(normalizedPositions, firstHypotheticalC3);
    const double randomC3Csm = averageRandomCsm(N, prng,
      [](const PositionCollection& positions) -> double {
        return csm::optimize(
          positions,
          elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), 3)
        ).first;
      }
    );

    const double firstC3Certainty = std::max(randomC3Csm - firstC3Result.first, 0.0) / randomC3Csm;
    if(firstC3Certainty >= 0.5 && angleIsRoughlyRight(firstC3Result.second.axis, mainAxis.axis)) {
      /* TODO keep looking for the other C3 axes by optimizing the result of
       * rotating the found secondary axis by the main axis
       */
      // Yes branch for 4 C3. Look for i
      cumulativeCertainty *= firstC3Certainty;
      const double iCertainty = predicates::inversion(normalizedPositions, prng);
      if(iCertainty >= 0.5) {
        // Yes branch for i
        return {PointGroup::Th, cumulativeCertainty * iCertainty};
      }

      // No branch for i, look for 6 sigma
      Eigen::Matrix<double, 3, 4> C3Axes;
      C3Axes.col(0) = mainAxis.axis;
      C3Axes.col(1) = firstC3Result.second.axis;
      C3Axes.col(2) = mainAxis.matrix() * C3Axes.col(1);
      C3Axes.col(3) = mainAxis.matrix() * C3Axes.col(2);
      const double randomSigmaCsm = averageRandomCsm(N, prng,
        [](const PositionCollection& positions) -> double {
          return csm::optimize(
            positions,
            elements::Reflection {Eigen::Vector3d::UnitZ()}
          ).first;
        }
      );
      const double sixSigmaCertainty = temple::accumulate(
        temple::adaptors::transform(
          temple::adaptors::allPairs(
            temple::adaptors::range(4)
          ),
          [&](const unsigned i, const unsigned j) -> double {
            const Eigen::Vector3d planeNormal = C3Axes.col(i).cross(C3Axes.col(j));
            // TODO element instead of optimize here to avoid the spatial checking mess after optimization?
            const double sigmaCsm = csm::element(
              normalizedPositions,
              elements::Reflection {planeNormal}
            );
            return std::max(randomSigmaCsm - sigmaCsm, 0.0) / randomSigmaCsm;
          }
        ),
        1.0,
        std::multiplies<>()
      );

      if(sixSigmaCertainty >= 0.5) {
        return {PointGroup::Td, cumulativeCertainty * sixSigmaCertainty};
      }

      return {PointGroup::T, cumulativeCertainty * (1 - sixSigmaCertainty)};
    }

    // Checking and not finding a C3 axis does influence the branching
    cumulativeCertainty *= (1 - firstC3Certainty);
  }

  return mediumPointGroups(normalizedPositions, mainAxis, N, prng, cumulativeCertainty);
}

} // namespace Symmetry
} // namespace Scine
