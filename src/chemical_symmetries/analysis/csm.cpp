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

#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"
#include "temple/Stringify.h"
#include "temple/constexpr/JSF.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Core>
#include <random>
#include <iostream>

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

constexpr unsigned nExperiments = 100;

struct RScriptWriter {
  std::ofstream file;

  RScriptWriter() : file("uniform_csm_data.R") {
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
    file << "results[" << symmetryIndex << ",] <- c(" << temple::condense(results) << ")\n";
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

  RScriptWriter writer;
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

  if(showElements) {
    for(unsigned N = 2; N < 8; ++N) {
      const auto values = temple::map(
        temple::adaptors::range(2 * nExperiments),
        [&](unsigned /* i */) -> double {
          auto normalized = detail::normalize(generateCoordinates(N, prng));
          Top top = standardizeTop(normalized);
          if(top == Top::Asymmetric) {
            reorientAsymmetricTop(normalized);
          }
          return csm::element(normalized, elements::Inversion {});
        }
      );
      std::cout << "CSM(i, " << N << ") = " << temple::average(values) << " +- " << temple::stddev(values) << "\n";
    }

    for(unsigned N = 2; N < 8; ++N) {
      const auto values = temple::map(
        temple::adaptors::range(2 * nExperiments),
        [&](unsigned /* i */) -> double {
          auto normalized = detail::normalize(generateCoordinates(N, prng));
          Top top = standardizeTop(normalized);
          if(top == Top::Asymmetric) {
            reorientAsymmetricTop(normalized);
          }
          return csm::optimizeCinf(normalized);
        }
      );
      std::cout << "CSM(Cinf, " << N << ") = " << temple::average(values) << " +- " << temple::stddev(values) << "\n";
    }
  } else {
    const std::map<Name, PointGroup> expected {
      {Name::Linear, PointGroup::Cinfv},
      {Name::Bent, PointGroup::C2v},
      {Name::TrigonalPlanar, PointGroup::D3h},
      {Name::CutTetrahedral, PointGroup::C3v},
      {Name::TShaped, PointGroup::C2v},
      {Name::Tetrahedral, PointGroup::Td},
      {Name::SquarePlanar, PointGroup::D4h},
      {Name::Seesaw, PointGroup::C2v},
      {Name::TrigonalPyramidal, PointGroup::C3v},
      {Name::SquarePyramidal, PointGroup::C4v},
      {Name::TrigonalBiPyramidal, PointGroup::D3h},
      {Name::PentagonalPlanar, PointGroup::D5h},
      {Name::Octahedral, PointGroup::Oh},
      {Name::TrigonalPrismatic, PointGroup::D3h},
      {Name::PentagonalPyramidal, PointGroup::C5v},
      {Name::PentagonalBiPyramidal, PointGroup::D5h}
    };

    std::cout << "Average CSM for uniform coordinates in sphere:\n";

    for(const auto& pair : expected) {
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
            return csm::pointGroup(normalized, pair.second);
          }
        );
        const double csmAverage = temple::average(values);
        const double csmStddev = temple::stddev(values);

        writer.addResults(pair.first, values);
        std::cout << name(pair.first) << " - " << N << ": " << csmAverage << " +- " << csmStddev << "\n";
      }
    }
  }

  return 0;
}
