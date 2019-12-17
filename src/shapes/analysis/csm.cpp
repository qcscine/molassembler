/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "shapes/Data.h"
#include "shapes/ContinuousMeasures.h"
#include "shapes/InertialMoments.h"
#include "shapes/Diophantine.h"

#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"
#include "temple/Stringify.h"
#include "temple/constexpr/Jsf.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Core>
#include <random>
#include <iostream>
#include <iomanip>

using namespace Scine;
using namespace Shapes;

template<typename PRNG>
Eigen::Vector3d randomVectorOnSphere(const double radius, PRNG& prng) {
  std::normal_distribution<double> normal {};
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << normal(prng),
         normal(prng),
         normal(prng);
  }

  return radius * v / v.norm();
}

template<typename PRNG>
Eigen::Vector3d randomVectorInSphere(const double radius, PRNG& prng) {
  std::uniform_real_distribution<double> uniform {};
  std::normal_distribution<double> normal {};
  const double u = std::cbrt(uniform(prng));
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << normal(prng),
         normal(prng),
         normal(prng);
  }

  return radius * u * v / v.norm();
}

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
Eigen::Matrix<double, 3, Eigen::Dynamic> generateCoordinates(unsigned P, PRNG& prng) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> positions(3, P + 1);
  positions.col(0) = Eigen::Vector3d::Zero();
  for(unsigned i = 1; i <= P; ++i) {
    positions.col(i) = normallyDistributedVectorInSphere(1.0, prng);
  }
  return positions;
}

constexpr unsigned nExperiments = 1000;

template<typename PRNG, typename F>
std::vector<double> averageRandomCsm(const unsigned N, PRNG& prng, F&& f) {
  assert(N >= 2);
  return temple::map(
    temple::adaptors::range(nExperiments),
    [&](unsigned /* i */) -> double {
      auto normalized = continuous::normalize(generateCoordinates(N, prng));
      Top top = standardizeTop(normalized);
      if(top == Top::Asymmetric) {
        reorientAsymmetricTop(normalized);
      }
      return f(normalized);
    }
  );
}

std::ostream& operator << (std::ostream& os, const std::vector<double>& values) {
  const auto end = std::end(values);
  for(auto it = std::begin(values); it != end; ++it) {
    os << *it;
    if(it != end - 1) {
      os << ", ";
    }
  }

  return os;
}

/* Does it make more sense to just maximize the CSM using two parameters at
 * each point and use LBFGS with numerical derivatives? More principled than
 * many random points, no guarantee maximum CSM is actually found.
 *
 * Procedure for finding upper bounds on the CSM for which clear distinctions
 * can be made between unrelated point groups (unrelated meaning neither is a
 * subset of the other).
 * - Generate points P for which CSM(P, D3h) = 0, and points Q for which CSM(Q, Oh) = 0.
 * - Add uniform direction fixed magnitude (fuzz) distortions to P and Q, making P' and Q'
 * - Find the maximum CSM(P', D3h) and maximum CSM(Q', Oh)
 * - Find the minimum CSM(Q', D3h) and minimum CSM(P', Oh)
 * - If the maximum CSM(P', D3h) < minimum CSM(Q', D3h)
 *   and the maximum CSM(Q', Oh) < minimum CSM(P', Oh),
 *   increase the fuzz magnitude
 *
 * To find the maximum and minimum CSM for a particular fuzziness of a set of
 * group symmetric points by numerical optimization:
 * - Generate A-symmetric points: Go through each solution of the diophantine
 *   for the usable group sizes and the required number of points (exclude
 *   groups whose probe point is the origin):
 *   - If only a single group size has multiplier > 0, generate group-symmetric
 *     points by applying a symmetry element from each group to the normalized
 *     probe point.
 *   - If there are multiple group sizes with multiplier > 0, the norms of each
 *     probe point are meta-parameters that need to be optimized too.
 */

struct RScriptWriter {
  std::ofstream file;

  RScriptWriter(std::string str) : file(str) {
    writeHeader();
  }

  void writeHeader() {
    file << "symmetryNames <- c(\"" << temple::condense(
      temple::map(Shapes::allShapes, [](auto name) { return Shapes::name(name); }),
      "\",\""
    ) << "\")\n";
    file << "symmetrySizes <- c(" << temple::condense(
      temple::map(Shapes::allShapes, [](auto name) { return Shapes::size(name); })
    ) << ")\n";
    file << "results <- array(numeric(), c(" << Shapes::allShapes.size() << ", " << nExperiments << "))\n";
  }

  void writeSeed(int seed) {
    file << "seed <- " << seed << "\n";
  }

  void addResults(const Shapes::Shape name, const std::vector<double>& results) {
    const unsigned symmetryIndex = nameIndex(name) + 1;
    file << "results[" << symmetryIndex << ",] <- c(" << results << ")\n";
  }

  template<typename F, typename PRNG>
  void addElementArray(const std::string& nameBase, F&& f, PRNG&& prng) {
    file << std::scientific;
    file << nameBase << "Array <- array(numeric(), c(7, " << nExperiments << "))\n";

    std::array<int, 7> seeds;
    for(unsigned i = 0; i < 7; ++i) {
      seeds[i] = prng();
    }

#pragma omp parallel for
    for(unsigned N = 2; N <= 8; ++N) {
      temple::jsf::JSF64 localPrng {seeds.at(N - 2)};
      const auto values = averageRandomCsm(N, localPrng, std::forward<F>(f));

#pragma omp critical
      {
        std::cout << "CSM(" << nameBase << ", " << N << ") = " << temple::average(values) << " +- " << temple::stddev(values) << "\n";
        file << nameBase << "Array[" << (N - 1) << ",] <- c(" << values << ")\n";
      }
    }

    file << "elementArrays[[\"" << nameBase << "\"]] <- " << nameBase << "Array)\n";
  }
};

int main(int argc, char* argv[]) {
  bool showElements = false;
  /* Set up program options */
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help,h", "Produce help message")
    (
      "seed,s",
      boost::program_options::value<int>(),
      "Seed to initialize PRNG with."
    )
    (
      "elements,e",
      boost::program_options::bool_switch(&showElements),
      "Show element CSM statistics instead of point groups"
    )
  ;

  /* Parse */
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv).
    options(options_description).
    style(
      boost::program_options::command_line_style::unix_style
      | boost::program_options::command_line_style::allow_long_disguise
    ).run(),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  if(options_variables_map.count("help") > 0) {
    std::cout << options_description << "\n";
    return 0;
  }

  RScriptWriter writer {
    showElements ? "elements.R" : "point_groups_data.R"
  };
  temple::jsf::JSF64 prng;
  if(options_variables_map.count("seed")) {
    const int seed = options_variables_map["seed"].as<int>();
    prng.seed(seed);
    writer.writeSeed(seed);
    std::cout << "PRNG seeded from parameters: " << seed << ".\n";
  } else {
    std::random_device randomDevice;
    const int seed = std::random_device {}();
    std::cout << "PRNG seeded from random_device: " << seed << ".\n";
    prng.seed(seed);
    writer.writeSeed(seed);
  }

  writer.file << std::scientific;

  if(showElements) {
    writer.file << "elementArrays <- list()\n";

    /* Inversion */
    writer.addElementArray(
      "inversion",
      [](const continuous::PositionCollection& positions) -> double {
        return continuous::element(positions, elements::Inversion {});
      },
      prng
    );

    /* Cinf */
    writer.addElementArray(
      "Cinf",
      [](const continuous::PositionCollection& positions) -> double {
        return continuous::Cinf(positions);
      },
      prng
    );

    /* Sigma */
    writer.addElementArray(
      "sigma",
      [](const continuous::PositionCollection& positions) -> double {
        return continuous::element(positions, elements::Reflection {Eigen::Vector3d::UnitZ()}).first;
      },
      prng
    );

    /* Cn axes */
    for(unsigned order = 2; order <= 8; ++order) {
      writer.addElementArray(
        "C" + std::to_string(order),
        [order](const continuous::PositionCollection& positions) -> double {
          return continuous::element(positions,
            elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), order)
          ).first;
        },
        prng
      );
    }

    /* Sn axes */
    for(unsigned order = 4; order <= 8; order += 2) {
      writer.addElementArray(
        "S" + std::to_string(order),
        [order](const continuous::PositionCollection& positions) -> double {
          return continuous::element(positions,
            elements::Rotation::Sn(Eigen::Vector3d::UnitZ(), order)
          ).first;
        },
        prng
      );
    }
  }

  if(!showElements) {
    std::cout << "Average CSM for uniform coordinates in sphere:\n";
    for(const Shape shape : allShapes) {
      const PointGroup group = pointGroup(shape);
      for(unsigned N = 2; N < 8; ++N) {
        /* Generate 100 random coordinates within a uniform sphere for each
         * symmetry and evaluate the CSM
         */
        const auto values = temple::map(
          temple::adaptors::range(nExperiments),
          [&](unsigned /* i */) -> double {
            auto normalized = continuous::normalize(generateCoordinates(N, prng));
            Top top = standardizeTop(normalized);
            if(top == Top::Asymmetric) {
              reorientAsymmetricTop(normalized);
            }
            return continuous::pointGroup(normalized, group);
          }
        );
        const double csmAverage = temple::average(values);
        const double csmStddev = temple::stddev(values);

        writer.addResults(shape, values);
        std::cout << name(shape) << " - " << N << ": " << csmAverage << " +- " << csmStddev << "\n";
      }
    }
  }

  return 0;
}
