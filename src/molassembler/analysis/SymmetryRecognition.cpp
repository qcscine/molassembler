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
#include "chemical_symmetries/ContinuousMeasures.h"
#include "chemical_symmetries/InertialMoments.h"
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
 * - Ideas
 *   - CShM classification
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
  virtual Symmetry::Shape identify(const Positions& positions) const = 0;
  virtual std::string name() const = 0;
};

struct AngularDeviation {
  template<typename F>
  static Symmetry::Shape identify(const Positions& positions, F&& f) {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Shape>;

    return temple::accumulate(
      Symmetry::allShapes,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Shape::Line},
      [&](const CarryType& carry, const Symmetry::Shape name) -> CarryType {
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

  Symmetry::Shape identify(const Positions& positions) const final {
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

  Symmetry::Shape identify(const Positions& positions) const final {
    return AngularDeviation::identify(positions, SquareFunctor {});
  }

  std::string name() const final {
    return "Pure angular deviation square";
  }
};

struct AngularDeviationGeometryIndexHybrid final : public Recognizer {
  Symmetry::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Shape>;

    /* Exclude symmetries using geometry indices if a relevant size */
    std::vector<Symmetry::Shape> excludedSymmetries;
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
          excludedSymmetries.push_back(Symmetry::Shape::Seesaw);
          excludedSymmetries.push_back(Symmetry::Shape::Tetrahedron);
        } else if(0.12 <= tau && tau < 0.62) {
          excludedSymmetries.push_back(Symmetry::Shape::Square);
          // Symmetry is seesaw
          excludedSymmetries.push_back(Symmetry::Shape::Tetrahedron);
        } else if(0.62 <= tau) {
          excludedSymmetries.push_back(Symmetry::Shape::Square);
          excludedSymmetries.push_back(Symmetry::Shape::Seesaw);
          // Symmetry is tetrahedral
        }
      } else if(S == 5) {
        /* Thresholds:
         * - τ₅ = 0 -> Symmetry is square pyramidal
         * - τ₅ = 1 -> Symmetry is trigonal bipyramidal
         */

        if(tau < 0.5) {
          excludedSymmetries.push_back(Symmetry::Shape::TrigonalBipyramid);
        } else if(tau > 0.5) {
          excludedSymmetries.push_back(Symmetry::Shape::SquarePyramid);
        }
      }
    }

    return temple::accumulate(
      Symmetry::allShapes,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Shape::Line},
      [&](const CarryType& carry, const Symmetry::Shape name) -> CarryType {
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

const std::map<Symmetry::Shape, Symmetry::PointGroup> pointGroupMapping {
  {Symmetry::Shape::Line, Symmetry::PointGroup::Dinfh},
  {Symmetry::Shape::Bent, Symmetry::PointGroup::C2v},
  {Symmetry::Shape::EquilateralTriangle, Symmetry::PointGroup::D3h},
  {Symmetry::Shape::VacantTetrahedron, Symmetry::PointGroup::C3v},
  {Symmetry::Shape::T, Symmetry::PointGroup::C2v},
  {Symmetry::Shape::Tetrahedron, Symmetry::PointGroup::Td},
  {Symmetry::Shape::Square, Symmetry::PointGroup::D4h},
  {Symmetry::Shape::Seesaw, Symmetry::PointGroup::C2v},
  {Symmetry::Shape::TrigonalPyramid, Symmetry::PointGroup::C3v},
  {Symmetry::Shape::SquarePyramid, Symmetry::PointGroup::C4v},
  {Symmetry::Shape::TrigonalBipyramid, Symmetry::PointGroup::D3h},
  {Symmetry::Shape::Pentagon, Symmetry::PointGroup::D5h},
  {Symmetry::Shape::Octahedron, Symmetry::PointGroup::Oh},
  {Symmetry::Shape::TrigonalPrism, Symmetry::PointGroup::D3h},
  {Symmetry::Shape::PentagonalPyramid, Symmetry::PointGroup::C5v},
  {Symmetry::Shape::PentagonalBipyramid, Symmetry::PointGroup::D5h},
  {Symmetry::Shape::SquareAntiprism, Symmetry::PointGroup::D4d}
};

struct PureCSM final : public Recognizer {
  Symmetry::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Symmetry::Shape>;

    Positions normalized = Symmetry::continuous::normalize(positions);
    const Symmetry::Top top = Symmetry::standardizeTop(normalized);
    if(top == Symmetry::Top::Asymmetric) {
      Symmetry::reorientAsymmetricTop(normalized);
    }

    return temple::accumulate(
      Symmetry::allShapes,
      CarryType {std::numeric_limits<double>::max(), Symmetry::Shape::Line},
      [&](const CarryType& bestPair, const Symmetry::Shape name) -> CarryType {
        if(Symmetry::size(name) != S) {
          return bestPair;
        }

        const double csm = Symmetry::continuous::pointGroup(
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

constexpr unsigned maxShapeSize = 6;

struct CShMPathDev final : public Recognizer {
  std::vector<Symmetry::Shape> validShapes;
  Eigen::MatrixXd minimumDistortionAngles;

  CShMPathDev() {
    for(const Symmetry::Shape shape : Symmetry::allShapes) {
      if(Symmetry::size(shape) <= maxShapeSize) {
        validShapes.push_back(shape);
      }
    }

    const unsigned N = validShapes.size();
    minimumDistortionAngles.resize(N, N);
    minimumDistortionAngles.setZero();
    for(unsigned i = 0; i < N; ++i) {
      Symmetry::Shape iShape = validShapes.at(i);

      for(unsigned j = i + 1; j < N; ++j) {
        Symmetry::Shape jShape = validShapes.at(j);

        if(Symmetry::size(iShape) == Symmetry::size(jShape)) {
          minimumDistortionAngles(i, j) = Symmetry::continuous::minimumDistortionAngle(
            iShape,
            jShape
          );
        }
      }
    }

    std::cout << "valid shapes:" << temple::stringify(
      temple::map(validShapes, [](auto x) { return Symmetry::name(x); })
    ) << "\n";
    std::cout << "minimum distortion angles:\n" << minimumDistortionAngles << "\n";
  }

  double minimumDistortionAngle(const Symmetry::Shape a, const Symmetry::Shape b) const {
    auto indexOfShape = [&](const Symmetry::Shape shape) -> unsigned {
      auto findIter = temple::find(validShapes, shape);
      if(findIter == std::end(validShapes)) {
        throw "Shape not found in valid shapes";
      }
      return findIter - std::begin(validShapes);
    };
    unsigned i, j;
    std::tie(i, j) = std::minmax(
      indexOfShape(a),
      indexOfShape(b)
    );

    return minimumDistortionAngles(i, j);
  }

  Symmetry::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    // Select shapes of matching size
    std::vector<Symmetry::Shape> matchingSizeShapes;
    for(const auto shape : Symmetry::allShapes) {
      if(Symmetry::size(shape) == S) {
        matchingSizeShapes.push_back(shape);
      }
    }

    // Calculate continuous shape measures for all selected shapes
    Positions normalized = Symmetry::continuous::normalize(positions);
    auto shapeMeasures = temple::map(
      matchingSizeShapes,
      [&](const Symmetry::Shape shape) -> double {
        return Symmetry::continuous::shape(normalized, shape);
      }
    );

    /* Calculate minimum distortion path deviations for all pairs, selecting
     * that pair for which the minimum distortion path deviation is minimal
     */
    using CarryType = std::tuple<unsigned, unsigned, double>;
    auto closestTuple = temple::accumulate(
      temple::adaptors::allPairs(temple::adaptors::range(matchingSizeShapes.size())),
      CarryType {0, 0, std::numeric_limits<double>::max()},
      [&](const CarryType& carry, const std::pair<unsigned, unsigned>& p) -> CarryType {
        const double cshm_a = shapeMeasures.at(p.first);
        const double cshm_b = shapeMeasures.at(p.second);

        const double deviation = (
          std::asin(std::sqrt(cshm_a) / 10) + std::asin(std::sqrt(cshm_b) / 10)
        ) / minimumDistortionAngle(
          matchingSizeShapes.at(p.first),
          matchingSizeShapes.at(p.second)
        ) - 1;

        if(deviation < std::get<2>(carry)) {
          return CarryType {p.first, p.second, deviation};
        }

        return carry;
      }
    );

    // Choose that shape from the best pair with minimal shape measure
    const unsigned a = std::get<0>(closestTuple);
    const unsigned b = std::get<1>(closestTuple);
    if(shapeMeasures.at(a) < shapeMeasures.at(b)) {
      return matchingSizeShapes.at(a);
    }

    return matchingSizeShapes.at(b);
  }

  std::string name() const final {
    return "CShM min dist path dev";
  }
};

struct Random final : public Recognizer {
  std::reference_wrapper<PRNG> prngRef;

  Random(PRNG& globalPRNG) : prngRef(globalPRNG) {}

  Symmetry::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    std::vector<Symmetry::Shape> viableSymmetries;
    for(const auto& name : Symmetry::allShapes) {
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

using RecognizersTuple = boost::mpl::list<PureAngularDeviationSquare, AngularDeviationGeometryIndexHybrid, CShMPathDev>;
constexpr std::size_t nRecognizers = boost::mpl::size<RecognizersTuple>::value;
using RecognizersList = std::vector<std::unique_ptr<Recognizer>>;

constexpr unsigned nDistortionValues = 11;
constexpr unsigned nRepeats = 100;
constexpr double maxDistortion = 1.0;

struct RScriptWriter {
  std::ofstream file;

  RScriptWriter() : file("recognition_data.R") {}

  void writeHeader(const RecognizersList& recognizers) {
    file << "nDistortionValues <- " << nDistortionValues << "\n";
    file << "maxDistortion <- " << maxDistortion << "\n";
    file << "nRepeats <- " << nRepeats << "\n";
    file << "shapeNames <- c(\"" << temple::condense(
      temple::map(Symmetry::allShapes, [](auto name) { return Symmetry::name(name); }),
      "\",\""
    ) << "\")\n";
    file << "symmetrySizes <- c(" << temple::condense(
      temple::map(Symmetry::allShapes, [](auto name) { return Symmetry::size(name); })
    ) << ")\n";
    file << "recognizers <- c(\"" << temple::condense(temple::map(recognizers, [](const auto& f) {return f->name();}), "\",\"") << "\")\n";
    file << "x.values <- c(0:" << (nDistortionValues - 1) << ") * " << maxDistortion << " / " << (nDistortionValues - 1) << "\n";
    file << "results <- array(numeric(), c(" << nRepeats << ", " << nDistortionValues << ", " << nRecognizers << ", " << Symmetry::allShapes.size() << "))\n";
  }

  void writeSeed(int seed) {
    file << "seed <- " << seed << "\n";
  }

  void addResults(const Symmetry::Shape name, const std::vector<std::vector<std::vector<Symmetry::Shape>>>& results) {
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
                  [](Symmetry::Shape n) -> unsigned {
                    return static_cast<std::underlying_type<Symmetry::Shape>::type>(n);
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

  std::vector<std::unique_ptr<Recognizer>> recognizerPtrs;
  boost::mpl::for_each<RecognizersTuple, boost::mpl::make_identity<boost::mpl::_1>>(
    [&recognizerPtrs](auto x) {
      using T = typename decltype(x)::type;
      recognizerPtrs.emplace_back(std::make_unique<T>());
    }
  );

  const unsigned R = recognizerPtrs.size();

  scriptFile.writeHeader(recognizerPtrs);

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

  for(const Symmetry::Shape name : Symmetry::allShapes) {
    const unsigned S = Symmetry::size(name);
    if(S > maxShapeSize) {
      // Skip anything bigger than maximum shape size for now
      break;
    }

    unsigned symmetriesOfSameSize = temple::accumulate(
      Symmetry::allShapes,
      0u,
      [&S](const unsigned carry, Symmetry::Shape n) -> unsigned {
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
        std::vector<Symmetry::Shape>
      >
    > recognizerCounts(R);

    for(unsigned i = 0; i < nDistortionValues; ++i) {
      for(unsigned r = 0; r < R; ++r) {
        recognizerCounts.at(r).push_back(std::vector<Symmetry::Shape> {});
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
