#include "chemical_symmetries/Flowchart.h"

#include "temple/constexpr/JSF.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Adaptors/Zip.h"
#include "temple/Functional.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Stringify.h"

#include <iostream>

namespace Scine{
namespace Symmetry {

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

std::pair<PointGroup, double> flowchart(const PositionCollection& normalizedPositions) {
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
    // Yes branch for Cinf, test i
    const double inversionCsm = csm::element(normalizedPositions, elements::Inversion {});
    const double randomInversionCsm = averageRandomCsm(N, prng,
      [](const PositionCollection& positions) -> double {
        return csm::element(positions, elements::Inversion {});
      }
    );
    const double iCertainty = std::max(randomInversionCsm - inversionCsm, 0.0) / randomInversionCsm;
    std::cout << "CSM(i) = " << inversionCsm << ", random CSM for " << N << " points: " << randomInversionCsm << ". Certainty = " << iCertainty << "\n";

    if(iCertainty >= 0.5) {
      return {PointGroup::Dinfh, cumulativeCertainty * iCertainty};
    }

    return {PointGroup::Cinfv, cumulativeCertainty * (1 - iCertainty)};
  }

  cumulativeCertainty *= (1 - cinfCertainty);

  /* Test for Cn axes up to n = 8 */
  const auto randomOrderCn = temple::map(
    temple::adaptors::range(2, 9),
    [&](const unsigned order) -> double {
      return averageRandomCsm(N, prng,
        [order](const PositionCollection& positions) -> double {
          return csm::element(
            positions,
            elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), order)
          );
        }
      );
    }
  );

  const auto bestOrderCn = temple::map(
    temple::adaptors::range(2, 9),
    [&](const unsigned order) -> std::pair<double, elements::Rotation> {
      // Try x, y, and z for the specified order
      const auto bestAxes = temple::map(
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

      // Return the best axis out of those
      return *std::min_element(
        std::begin(bestAxes),
        std::end(bestAxes),
        [](const auto& a, const auto& b) -> bool {
          return a.first < b.first;
        }
      );
    }
  );

  const auto orderCertainties = temple::map(
    temple::adaptors::zip(bestOrderCn, randomOrderCn),
    [](const auto& csmRotationPair, const double randomCsm) -> double {
      return std::max(randomCsm - csmRotationPair.first, 0.0) / randomCsm;
    }
  );

  std::cout << temple::stringify(randomOrderCn) << "\n";
  std::cout << temple::stringify(temple::map(bestOrderCn, temple::functor::first)) << "\n";
  std::cout << temple::stringify(orderCertainties) << "\n";

  // No branch for Cinf, find major Cn up to C8
  return {PointGroup::C1, 0.0};
}

} // namespace Symmetry
} // namespace Scine
