/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "chemical_symmetries/Symmetries.h"
#include "chemical_symmetries/Recognition.h"
#include "chemical_symmetries/Diophantine.h"

#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"
#include "temple/Stringify.h"
#include "temple/constexpr/JSF.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Core>
#include <random>
#include <iostream>
#include <iomanip>

using namespace Scine;
using namespace Symmetry;

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
      auto normalized = detail::normalize(generateCoordinates(N, prng));
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

std::vector<PositionCollection> groupPoints(
  const elements::ElementsList elements,
  const elements::NPGroupingsMapType& npGroupings,
  const unsigned N
) {
  std::vector<unsigned> usableGroupSizes;
  for(const auto& iterPair : npGroupings) {
    if(
      iterPair.first <= N
      && !temple::all_of(
        iterPair.second,
        [](const auto& elementGrouping) {
          return elementGrouping.probePoint.isApprox(Eigen::Vector3d::Zero(), 1e-5);
        }
      )
    ) {
      usableGroupSizes.push_back(iterPair.first);
    }
  }

  std::vector<unsigned> groupSizeMultipliers;
  if(!diophantine::has_solution(usableGroupSizes, N)) {
    throw std::logic_error("No P = N groupings for this point group");
  }

  if(!diophantine::first_solution(groupSizeMultipliers, usableGroupSizes, N)) {
    throw std::logic_error("Diophantine failed");
  }

  do {

  } while(diophantine::next_solution(groupSizeMultipliers, usableGroupSizes, N));
}

struct RScriptWriter {
  std::ofstream file;

  RScriptWriter(std::string str) : file(str) {
    writeHeader();
  }

  void writeHeader() {
    file << "symmetryNames <- c(\"" << temple::condense(
      temple::map(Symmetry::allNames, [](auto name) { return Symmetry::name(name); }),
      "\",\""
    ) << "\")\n";
    file << "symmetrySizes <- c(" << temple::condense(
      temple::map(Symmetry::allNames, [](auto name) { return Symmetry::size(name); })
    ) << ")\n";
    file << "results <- array(numeric(), c(" << Symmetry::allNames.size() << ", " << nExperiments << "))\n";
  }

  void writeSeed(int seed) {
    file << "seed <- " << seed << "\n";
  }

  void addResults(const Symmetry::Name name, const std::vector<double>& results) {
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
  bool fuzzPointGroups = false;
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
    (
      "fuzz,f",
      boost::program_options::bool_switch(&fuzzPointGroups),
      "Fuzz point groups"
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
      [](const PositionCollection& positions) -> double {
        return csm::element(positions, elements::Inversion {});
      },
      prng
    );

    /* Cinf */
    writer.addElementArray(
      "Cinf",
      [](const PositionCollection& positions) -> double {
        return csm::optimizeCinf(positions);
      },
      prng
    );

    /* Sigma */
    writer.addElementArray(
      "sigma",
      [](const PositionCollection& positions) -> double {
        return csm::optimize(positions, elements::Reflection {Eigen::Vector3d::UnitZ()}).first;
      },
      prng
    );

    /* Cn axes */
    for(unsigned order = 2; order <= 8; ++order) {
      writer.addElementArray(
        "C" + std::to_string(order),
        [order](const PositionCollection& positions) -> double {
          return csm::optimize(positions,
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
        [order](const PositionCollection& positions) -> double {
          return csm::optimize(positions,
            elements::Rotation::Sn(Eigen::Vector3d::UnitZ(), order)
          ).first;
        },
        prng
      );
    }
  }

  if(fuzzPointGroups) {
    const unsigned N = 6;
    const unsigned nStructures = 100;
    const double fuzzDelta = 0.05;
    const double maxFuzz = 1.00;
    const std::vector<PointGroup> groups {
      PointGroup::Oh,
      PointGroup::D3h
    };

    for(const PointGroup group : groups) {
      /* Need to use diophantine of all groupings smaller or equal in size to N
       * to generate all types of structure compatible with the point group
       */
      auto elements = elements::symmetryElements(group);
      auto npGroupings = elements::npGroupings(elements);


      // Generate random point group symmetric structures
      auto structures = temple::map(
        temple::adaptors::range(100),
        [&](unsigned /* i */) -> PositionCollection {

        }
      );


      for(double fuzz = 0; fuzz <= maxFuzz; fuzz += fuzzDelta) {

      }
    }
  }

  if(!fuzzPointGroups && !showElements) {
    std::cout << "Average CSM for uniform coordinates in sphere:\n";
    for(const Name& symmetryName : allNames) {
      const PointGroup group = pointGroup(symmetryName);
      for(unsigned N = 2; N < 8; ++N) {
        /* Generate 100 random coordinates within a uniform sphere for each
         * symmetry and evaluate the CSM
         */
        const auto values = temple::map(
          temple::adaptors::range(nExperiments),
          [&](unsigned /* i */) -> double {
            auto normalized = detail::normalize(generateCoordinates(N, prng));
            Top top = standardizeTop(normalized);
            if(top == Top::Asymmetric) {
              reorientAsymmetricTop(normalized);
            }
            return csm::pointGroup(normalized, group);
          }
        );
        const double csmAverage = temple::average(values);
        const double csmStddev = temple::stddev(values);

        writer.addResults(symmetryName, values);
        std::cout << name(symmetryName) << " - " << N << ": " << csmAverage << " +- " << csmStddev << "\n";
      }
    }
  }

  return 0;
}
