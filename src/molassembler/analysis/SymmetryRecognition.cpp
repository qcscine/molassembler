/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "boost/mpl/list.hpp"
#include "boost/mpl/for_each.hpp"
#include "boost/mpl/size.hpp"

#include "molassembler/Detail/Cartesian.h"

#include "chemical_symmetries/Symmetries.h"
#include "chemical_symmetries/Recognition.h"
#include "chemical_symmetries/DynamicProperties.h"
#include "chemical_symmetries/TauCriteria.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/Stringify.h"
#include "temple/constexpr/JSF.h"
#include "temple/constexpr/Numeric.h"

#include <iostream>
#include <fstream>
#include <random>

/* TODO
 * - CSM based does not understand subsets of symmetries and has no threshold
 *   for it. Any marginally bent three-particle arrangement will not be
 *   classified as linear.
 *
 * - Ideas
 *   - Find maximum value of CSM for each point group and rescale CSMs accordingly
 *   - Find average value of CSM for random point cloud and rescale CSMs accordingly
 *   - Add the geometry index hybrid with angular deviation
 *   - Jan: Flowchart point groups instead of using CSM value comparisons,
 *     assign probabilities to decisions instead of using thresholds to follow
 *     single paths. Evaluate trees for all relevant point groups and then
 *     compare cumulative (multiplicative) probabilities.
 */

using namespace std::string_literals;
using namespace Scine;

using Positions = Eigen::Matrix<double, 3, Eigen::Dynamic>;
using PRNG = temple::jsf::JSF64;

struct Recognizer {
  virtual ~Recognizer() = default;

  /*! @brief Figure out which symmetry is present
   *
   * @param positions Positions of all particles of the symmetry. The first
   * column is the central particle of the symmetry.
   */
  virtual Symmetry::Name identify(const Positions& positions) const = 0;
  virtual std::string name() const = 0;
};

struct AngularDeviation {
  template<typename F>
  static Symmetry::Name identify(const Positions& positions, F&& f) {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Name>;

    return temple::accumulate(
      Symmetry::allNames,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Name::Linear},
      [&](const CarryType& carry, const Symmetry::Name name) -> CarryType {
        if(Symmetry::size(name) != S) {
          return carry;
        }
        const auto angleFunction = Symmetry::angleFunction(name);

        /* Minimize angular deviations over all rotations of maximally
         * asymmetric symmetry case
         */
        const double penalty = temple::accumulate(
          Symmetry::properties::generateAllRotations(name, temple::iota<unsigned>(S)),
          std::numeric_limits<double>::max(),
          [&](const double minAngularDeviation, const auto& rotation) -> double {
            const double angleDeviation = temple::sum(
              temple::adaptors::transform(
                temple::adaptors::allPairs(
                  temple::adaptors::range(S)
                ),
                [&](const unsigned siteI, const unsigned siteJ) -> double {
                  return f(
                    molassembler::cartesian::angle(
                      positions.col(1 + siteI),
                      positions.col(0),
                      positions.col(1 + siteJ)
                    ) - angleFunction(rotation.at(siteI), rotation.at(siteJ))
                  );
                }
              )
            );

            return std::min(minAngularDeviation, angleDeviation);
          }
        );

        if(penalty < carry.first) {
          return {penalty, name};
        }

        return carry;
      }
    ).second;
  }
};

struct PureAngularDeviationAbs final : public Recognizer {
  struct FabsFunctor {
    template<typename T>
    T operator() (T x) const {
      return std::fabs(x);
    }
  };

  Symmetry::Name identify(const Positions& positions) const final {
    return AngularDeviation::identify(positions, FabsFunctor {});
  }

  std::string name() const final {
    return "Pure angular deviation abs. value";
  }
};

struct PureAngularDeviationSquare final : public Recognizer {
  struct SquareFunctor {
    template<typename T>
    T operator() (T x) const {
      return x * x;
    }
  };

  Symmetry::Name identify(const Positions& positions) const final {
    return AngularDeviation::identify(positions, SquareFunctor {});
  }

  std::string name() const final {
    return "Pure angular deviation square";
  }
};

struct AngularDeviationGeometryIndexHybrid final : public Recognizer {
  Symmetry::Name identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Name>;

    /* Exclude symmetries using geometry indices if a relevant size */
    std::vector<Symmetry::Name> excludedSymmetries;
    if(S == 4 || S == 5) {
      const double tau = Symmetry::tau(
        temple::sort(
          temple::map(
            temple::adaptors::allPairs(temple::adaptors::range(S)),
            [&](const unsigned i, const unsigned j) -> double {
              return molassembler::cartesian::angle(
                positions.col(1 + i),
                positions.col(0),
                positions.col(1 + j)
              );
            }
          )
        )
      );

      if(S == 4) {
        /* Thresholds
         * - τ₄' = 0 -> Symmetry is square planar
         * - τ₄' = 0.24 -> Symmetry is seesaw
         * - τ₄' = 1 -> Symmetry is tetrahedral
         */
        if(tau < 0.12) {
          // Symmetry is square planar
          excludedSymmetries.push_back(Symmetry::Name::Seesaw);
          excludedSymmetries.push_back(Symmetry::Name::Tetrahedral);
        } else if(0.12 <= tau && tau < 0.62) {
          excludedSymmetries.push_back(Symmetry::Name::SquarePlanar);
          // Symmetry is seesaw
          excludedSymmetries.push_back(Symmetry::Name::Tetrahedral);
        } else if(0.62 <= tau) {
          excludedSymmetries.push_back(Symmetry::Name::SquarePlanar);
          excludedSymmetries.push_back(Symmetry::Name::Seesaw);
          // Symmetry is tetrahedral
        }
      } else if(S == 5) {
        /* Thresholds:
         * - τ₅ = 0 -> Symmetry is square pyramidal
         * - τ₅ = 1 -> Symmetry is trigonal bipyramidal
         */

        if(tau < 0.5) {
          excludedSymmetries.push_back(Symmetry::Name::TrigonalBiPyramidal);
        } else if(tau > 0.5) {
          excludedSymmetries.push_back(Symmetry::Name::SquarePyramidal);
        }
      }
    }

    return temple::accumulate(
      Symmetry::allNames,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Name::Linear},
      [&](const CarryType& carry, const Symmetry::Name name) -> CarryType {
        if(Symmetry::size(name) != S || temple::makeContainsPredicate(excludedSymmetries)(name)) {
          return carry;
        }
        const auto angleFunction = Symmetry::angleFunction(name);

        /* Minimize angular deviations over all rotations of maximally
         * asymmetric symmetry case
         */
        const double penalty = temple::accumulate(
          Symmetry::properties::generateAllRotations(name, temple::iota<unsigned>(S)),
          std::numeric_limits<double>::max(),
          [&](const double minAngularDeviation, const auto& rotation) -> double {
            const double angleDeviation = temple::sum(
              temple::adaptors::transform(
                temple::adaptors::allPairs(
                  temple::adaptors::range(S)
                ),
                [&](const unsigned siteI, const unsigned siteJ) -> double {
                  const double deviation = (
                    molassembler::cartesian::angle(
                      positions.col(1 + siteI),
                      positions.col(0),
                      positions.col(1 + siteJ)
                    ) - angleFunction(rotation.at(siteI), rotation.at(siteJ))
                  );
                  return deviation * deviation;
                }
              )
            );

            return std::min(minAngularDeviation, angleDeviation);
          }
        );

        if(penalty < carry.first) {
          return {penalty, name};
        }

        return carry;
      }
    ).second;
  }

  std::string name() const final {
    return "Angular deviation (squared) & geometry index hybrid";
  }
};

const std::map<Symmetry::Name, Symmetry::PointGroup> pointGroupMapping {
  {Symmetry::Name::Linear, Symmetry::PointGroup::Dinfh},
  {Symmetry::Name::Bent, Symmetry::PointGroup::C2v},
  {Symmetry::Name::TrigonalPlanar, Symmetry::PointGroup::D3h},
  {Symmetry::Name::CutTetrahedral, Symmetry::PointGroup::C3v},
  {Symmetry::Name::TShaped, Symmetry::PointGroup::C2v},
  {Symmetry::Name::Tetrahedral, Symmetry::PointGroup::Td},
  {Symmetry::Name::SquarePlanar, Symmetry::PointGroup::D4h},
  {Symmetry::Name::Seesaw, Symmetry::PointGroup::C2v},
  {Symmetry::Name::TrigonalPyramidal, Symmetry::PointGroup::C3v},
  {Symmetry::Name::SquarePyramidal, Symmetry::PointGroup::C4v},
  {Symmetry::Name::TrigonalBiPyramidal, Symmetry::PointGroup::D3h},
  {Symmetry::Name::PentagonalPlanar, Symmetry::PointGroup::D5h},
  {Symmetry::Name::Octahedral, Symmetry::PointGroup::Oh},
  {Symmetry::Name::TrigonalPrismatic, Symmetry::PointGroup::D3h},
  {Symmetry::Name::PentagonalPyramidal, Symmetry::PointGroup::C5v},
  {Symmetry::Name::PentagonalBiPyramidal, Symmetry::PointGroup::D5h},
  {Symmetry::Name::SquareAntiPrismatic, Symmetry::PointGroup::D4d}
};

struct PureCSM final : public Recognizer {
  Symmetry::Name identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Name>;

    Positions normalized = Symmetry::detail::normalize(positions);
    const Symmetry::Top top = Symmetry::standardizeTop(normalized);
    if(top == Symmetry::Top::Asymmetric) {
      Symmetry::reorientAsymmetricTop(normalized);
    }

    return temple::accumulate(
      Symmetry::allNames,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Name::Linear},
      [&](const CarryType& bestPair, const Symmetry::Name name) -> CarryType {
        if(Symmetry::size(name) != S) {
          return bestPair;
        }

        const double csm = Symmetry::csm::pointGroup(
          normalized,
          pointGroupMapping.at(name)
        );

        if(csm < bestPair.first) {
          return {csm, name};
        }

        return bestPair;
      }
    ).second;
  }

  std::string name() const final {
    return "Pure CSM";
  }
};

struct OrderRespectingCSM final : public Recognizer {
  Symmetry::Name identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Name>;

    Positions normalized = Symmetry::detail::normalize(positions);
    const Symmetry::Top top = Symmetry::standardizeTop(normalized);
    if(top == Symmetry::Top::Asymmetric) {
      Symmetry::reorientAsymmetricTop(normalized);
    }

    return temple::accumulate(
      Symmetry::allNames,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Name::Linear},
      [&](const CarryType& bestPair, const Symmetry::Name name) -> CarryType {
        if(Symmetry::size(name) != S) {
          return bestPair;
        }

        const double csm = Symmetry::csm::pointGroup(
          normalized,
          pointGroupMapping.at(name)
        );

        const unsigned thisPointGroupOrder = Symmetry::elements::order(pointGroupMapping.at(name));
        const unsigned carryPointGroupOrder = Symmetry::elements::order(pointGroupMapping.at(bestPair.second));

        if(csm / thisPointGroupOrder < bestPair.first / carryPointGroupOrder) {
          return {csm, name};
        }

        return bestPair;
      }
    ).second;
  }

  std::string name() const final {
    return "Order respecting CSM";
  }
};

struct PositionManipulatingCSM final : public Recognizer {
  Symmetry::Name identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Name>;

    Positions manipulated = positions;
    // Normalize all positions to norm 1 from the central atom
    for(unsigned i = 1; i < positions.cols(); ++i) {
      manipulated.col(i) = manipulated.col(0) + (manipulated.col(i) - manipulated.col(0)).normalized();
    }

    Positions normalized = Symmetry::detail::normalize(manipulated);
    const Symmetry::Top top = Symmetry::standardizeTop(normalized);
    if(top == Symmetry::Top::Asymmetric) {
      Symmetry::reorientAsymmetricTop(normalized);
    }

    return temple::accumulate(
      Symmetry::allNames,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Name::Linear},
      [&](const CarryType& bestPair, const Symmetry::Name name) -> CarryType {
        if(Symmetry::size(name) != S) {
          return bestPair;
        }

        const double csm = Symmetry::csm::pointGroup(
          normalized,
          pointGroupMapping.at(name)
        );

        if(csm < bestPair.first) {
          return {csm, name};
        }

        return bestPair;
      }
    ).second;
  }

  std::string name() const final {
    return "Position manipulating CSM";
  }
};

struct Random final : public Recognizer {
  std::reference_wrapper<PRNG> prngRef;

  Random(PRNG& globalPRNG) : prngRef(globalPRNG) {}

  Symmetry::Name identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    std::vector<Symmetry::Name> viableSymmetries;
    for(const auto& name : Symmetry::allNames) {
      if(Symmetry::size(name) == S) {
        viableSymmetries.push_back(name);
      }
    }
    return viableSymmetries.at(
      temple::random::getSingle<unsigned>(0, viableSymmetries.size() - 1, prngRef.get())
    );
  }

  std::string name() const final {
    return "Random";
  }
};


Eigen::Vector3d randomVectorOnSphere(const double radius, PRNG& prng) {
  std::normal_distribution<double> distribution {};
  Eigen::Vector3d v = Eigen::Vector3d::Zero();

  while(v.norm() < 0.01) {
    v << distribution(prng),
         distribution(prng),
         distribution(prng);
  }

  return radius * v / v.norm();
}

void distort(Eigen::Ref<Positions> positions, const double distortionNorm, PRNG& prng) {
  const unsigned N = positions.cols();
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) += randomVectorOnSphere(distortionNorm, prng);
  }
}

using RecognizersTuple = boost::mpl::list<PureAngularDeviationSquare, AngularDeviationGeometryIndexHybrid, PureCSM, OrderRespectingCSM>;
constexpr std::size_t nRecognizers = boost::mpl::size<RecognizersTuple>::value;

std::vector<std::string> recognizerNames() {
  std::vector<std::string> names;
  boost::mpl::for_each<RecognizersTuple>([&](auto x) { names.push_back(x.name()); });
  return names;
}

constexpr unsigned nDistortionValues = 11;
constexpr unsigned nRepeats = 100;
constexpr double maxDistortion = 1.0;

struct RScriptWriter {
  std::ofstream file;

  RScriptWriter() : file("recognition_data.R") {
    writeHeader();
  }

  void writeHeader() {
    file << "nDistortionValues <- " << nDistortionValues << "\n";
    file << "maxDistortion <- " << maxDistortion << "\n";
    file << "nRepeats <- " << nRepeats << "\n";
    file << "symmetryNames <- c(\"" << temple::condense(
      temple::map(Symmetry::allNames, [](auto name) { return Symmetry::name(name); }),
      "\",\""
    ) << "\")\n";
    file << "symmetrySizes <- c(" << temple::condense(
      temple::map(Symmetry::allNames, [](auto name) { return Symmetry::size(name); })
    ) << ")\n";
    file << "recognizers <- c(\"" << temple::condense(recognizerNames(), "\",\"") << "\")\n";
    file << "x.values <- c(0:" << (nDistortionValues - 1) << ") * " << maxDistortion << " / " << (nDistortionValues - 1) << "\n";
    file << "results <- array(numeric(), c(" << nRepeats << ", " << nDistortionValues << ", " << nRecognizers << ", " << Symmetry::allNames.size() << "))\n";
  }

  void writeSeed(int seed) {
    file << "seed <- " << seed << "\n";
  }

  void addResults(const Symmetry::Name name, const std::vector<std::vector<std::vector<Symmetry::Name>>>& results) {
    const unsigned symmetryIndex = Symmetry::nameIndex(name);
    for(unsigned recognizerIndex = 0; recognizerIndex < nRecognizers; ++recognizerIndex) {
      const std::string array = (
        "array(c("
        + temple::condense(
          temple::map(
            results.at(recognizerIndex),
            [](const auto& distortionResult) {
              return temple::condense(
                temple::map(
                  distortionResult,
                  [](Symmetry::Name n) -> unsigned {
                    return static_cast<std::underlying_type<Symmetry::Name>::type>(n);
                  }
                )
              );
            }
          )
        )
        + "), dim=c("
        + std::to_string(nRepeats)
        + ", "
        + std::to_string(nDistortionValues)
        + "))"
      );

      file << "results[,," << (recognizerIndex + 1) << ", " << (symmetryIndex + 1) << "] <- " << array << "\n";
    }
  }
};

int main(int argc, char* argv[]) {
  /* Set up program options */
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help,h", "Produce help message")
    (
      "seed,s",
      boost::program_options::value<int>(),
      "Seed to initialize PRNG with."
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

  /* Manage parse results */
  RScriptWriter scriptFile {};

  temple::jsf::JSF64 prng;
  if(options_variables_map.count("seed")) {
    const int seed = options_variables_map["seed"].as<int>();
    prng.seed(seed);
    std::cout << "PRNG seeded from parameters: " << seed << ".\n";
    scriptFile.writeSeed(seed);
  } else {
    std::random_device randomDevice;
    const int seed = std::random_device {}();
    std::cout << "PRNG seeded from random_device: " << seed << ".\n";
    prng.seed(seed);
    scriptFile.writeSeed(seed);
  }

  std::vector<std::unique_ptr<Recognizer>> recognizerPtrs;
  boost::mpl::for_each<RecognizersTuple>(
    [&recognizerPtrs](auto x) {
      using T = decltype(x);
      recognizerPtrs.push_back(std::make_unique<T>());
    }
  );

  const unsigned R = recognizerPtrs.size();

  for(const Symmetry::Name name : Symmetry::allNames) {
    const unsigned S = Symmetry::size(name);
    unsigned symmetriesOfSameSize = temple::accumulate(
      Symmetry::allNames,
      0u,
      [&S](const unsigned carry, Symmetry::Name n) -> unsigned {
        if(Symmetry::size(n) == S) {
          return carry + 1;
        }

        return carry;
      }
    );

    if(symmetriesOfSameSize == 1) {
      continue;
    }

    Positions basePositions (3, S + 1);
    basePositions.col(0) = Eigen::Vector3d::Zero();
    basePositions.rightCols(S) = Symmetry::symmetryData().at(name).coordinates;

    std::vector< // each recognizer has its own vector
      std::vector< // each distortion value
        std::vector<Symmetry::Name>
      >
    > recognizerCounts(R);

    for(unsigned i = 0; i < nDistortionValues; ++i) {
      for(unsigned r = 0; r < R; ++r) {
        recognizerCounts.at(r).push_back(std::vector<Symmetry::Name> {});
      }

      const double distortion = i * maxDistortion / (nDistortionValues - 1);
      for(unsigned j = 0; j < nRepeats; ++j) {
        Positions positions = basePositions;
        distort(positions, distortion, prng);

        for(unsigned r = 0; r < R; ++r) {
          recognizerCounts.at(r).back().push_back(
            recognizerPtrs.at(r)->identify(positions)
          );
        }
      }
    }

    scriptFile.addResults(name, recognizerCounts);

    std::cout << "Symmetry: " << Symmetry::name(name) << "\n";
  }
}
