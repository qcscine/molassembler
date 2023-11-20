/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "boost/mpl/list.hpp"
#include "boost/mpl/for_each.hpp"
#include "boost/mpl/size.hpp"

#include "Molassembler/Detail/Cartesian.h"

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Shapes/ContinuousMeasures.h"
#include "Molassembler/Shapes/InertialMoments.h"
#include "Molassembler/Shapes/Properties.h"
#include "Molassembler/Shapes/TauCriteria.h"

#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Adaptors/Iota.h"
#include "Molassembler/Temple/Adaptors/Transform.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/Jsf.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include <iostream>
#include <fstream>
#include <random>

using namespace std::string_literals;
using namespace Scine;
using namespace Molassembler;

using Positions = Eigen::Matrix<double, 3, Eigen::Dynamic>;
using PRNG = Temple::JSF64;

struct Recognizer {
  virtual ~Recognizer() = default;

  /*! @brief Figure out which shape is present
   *
   * @param positions Positions of all particles of the shape. The first
   * column is the central particle of the shape.
   */
  virtual Shapes::Shape identify(const Positions& positions) const = 0;
  virtual std::string name() const = 0;
};

struct AngularDeviation {
  template<typename F>
  static Shapes::Shape identify(const Positions& positions, F&& f) {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Shapes::Shape>;

    return Temple::accumulate(
      Shapes::allShapes,
      CarryType {std::numeric_limits<double>::max(), Shapes::Shape::Line},
      [&](const CarryType& carry, const Shapes::Shape name) -> CarryType {
        if(Shapes::size(name) != S) {
          return carry;
        }
        const auto angleFunction = Shapes::angleFunction(name);

        /* Minimize angular deviations over all rotations of maximally
         * asymmetric shape case
         */
        const double penalty = Temple::accumulate(
          Shapes::Properties::generateAllRotations(name, Temple::iota<Shapes::Vertex>(S)),
          std::numeric_limits<double>::max(),
          [&](const double minAngularDeviation, const auto& rotation) -> double {
            const double angleDeviation = Temple::sum(
              Temple::Adaptors::transform(
                Temple::Adaptors::allPairs(
                  Temple::Adaptors::range(S)
                ),
                [&](const unsigned siteI, const unsigned siteJ) -> double {
                  return f(
                    Molassembler::Cartesian::angle(
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

  Shapes::Shape identify(const Positions& positions) const final {
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

  Shapes::Shape identify(const Positions& positions) const final {
    return AngularDeviation::identify(positions, SquareFunctor {});
  }

  std::string name() const final {
    return "Pure angular deviation square";
  }
};

struct AngularDeviationGeometryIndexHybrid final : public Recognizer {
  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Shapes::Shape>;

    /* Exclude symmetries using geometry indices if a relevant size */
    std::vector<Shapes::Shape> excludedSymmetries;
    if(S == 4 || S == 5) {
      const double tau = Shapes::tau(
        Temple::sorted(
          Temple::map(
            Temple::Adaptors::allPairs(Temple::Adaptors::range(S)),
            [&](const unsigned i, const unsigned j) -> double {
              return Molassembler::Cartesian::angle(
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
         * - τ₄' = 0 -> shape is square planar
         * - τ₄' = 0.24 -> shape is seesaw
         * - τ₄' = 1 -> shape is tetrahedral
         */
        if(tau < 0.12) {
          // shape is square planar
          excludedSymmetries.push_back(Shapes::Shape::Seesaw);
          excludedSymmetries.push_back(Shapes::Shape::Tetrahedron);
        } else if(0.12 <= tau && tau < 0.62) {
          excludedSymmetries.push_back(Shapes::Shape::Square);
          // shape is seesaw
          excludedSymmetries.push_back(Shapes::Shape::Tetrahedron);
        } else if(0.62 <= tau) {
          excludedSymmetries.push_back(Shapes::Shape::Square);
          excludedSymmetries.push_back(Shapes::Shape::Seesaw);
          // shape is tetrahedral
        }
      } else if(S == 5) {
        /* Thresholds:
         * - τ₅ = 0 -> shape is square pyramidal
         * - τ₅ = 1 -> shape is trigonal bipyramidal
         */

        if(tau < 0.5) {
          excludedSymmetries.push_back(Shapes::Shape::TrigonalBipyramid);
        } else if(tau > 0.5) {
          excludedSymmetries.push_back(Shapes::Shape::SquarePyramid);
        }
      }
    }

    return Temple::accumulate(
      Shapes::allShapes,
      CarryType {std::numeric_limits<double>::max(), Shapes::Shape::Line},
      [&](const CarryType& carry, const Shapes::Shape name) -> CarryType {
        if(Shapes::size(name) != S || Temple::makeContainsPredicate(excludedSymmetries)(name)) {
          return carry;
        }
        const auto angleFunction = Shapes::angleFunction(name);

        /* Minimize angular deviations over all rotations of maximally
         * asymmetric shape case
         */
        const double penalty = Temple::accumulate(
          Shapes::Properties::generateAllRotations(name, Temple::iota<Shapes::Vertex>(S)),
          std::numeric_limits<double>::max(),
          [&](const double minAngularDeviation, const auto& rotation) -> double {
            const double angleDeviation = Temple::sum(
              Temple::Adaptors::transform(
                Temple::Adaptors::allPairs(
                  Temple::Adaptors::range(S)
                ),
                [&](const unsigned siteI, const unsigned siteJ) -> double {
                  const double deviation = (
                    Molassembler::Cartesian::angle(
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

struct PureCSM final : public Recognizer {
  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    using CarryType = std::pair<double, Shapes::Shape>;

    Positions normalized = Shapes::Continuous::normalize(positions);
    const Shapes::Top top = Shapes::standardizeTop(normalized);
    if(top == Shapes::Top::Asymmetric) {
      Shapes::reorientAsymmetricTop(normalized);
    }

    return Temple::accumulate(
      Shapes::allShapes,
      CarryType {std::numeric_limits<double>::max(), Shapes::Shape::Line},
      [&](const CarryType& bestPair, const Shapes::Shape name) -> CarryType {
        if(Shapes::size(name) != S) {
          return bestPair;
        }

        const double csm = Shapes::Continuous::pointGroup(
          normalized,
          Shapes::pointGroup(name)
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

struct CShM final : public Recognizer {
  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;

    Positions normalized = Shapes::Continuous::normalize(positions);
    using CarryType = std::pair<double, Shapes::Shape>;
    return Temple::accumulate(
      Shapes::allShapes,
      CarryType {std::numeric_limits<double>::max(), Shapes::Shape::Line},
      [&](const CarryType& carry, const Shapes::Shape shape) {
        if(Shapes::size(shape) != S) {
          return carry;
        }

        double shapeMeasure = Shapes::Continuous::shape(normalized, shape).measure;

        if(shapeMeasure < carry.first) {
          return CarryType {shapeMeasure, shape};
        }

        return carry;
      }
    ).second;
  }

  std::string name() const final {
    return "CShM";
  }
};

struct BiasedCShM final : public Recognizer {
  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;

    Positions normalized = Shapes::Continuous::normalize(positions);
    using CarryType = std::pair<double, Shapes::Shape>;
    return Temple::accumulate(
      Shapes::allShapes,
      CarryType {std::numeric_limits<double>::max(), Shapes::Shape::Line},
      [&](const CarryType& carry, const Shapes::Shape shape) {
        if(Shapes::size(shape) != S) {
          return carry;
        }

        double shapeMeasure = Shapes::Continuous::shape(normalized, shape).measure;

        if(shape == Shapes::Shape::TrigonalPyramid) {
          shapeMeasure *= 4;
        } else if(shape == Shapes::Shape::Seesaw) {
          shapeMeasure *= 2;
        }

        if(shapeMeasure < carry.first) {
          return CarryType {shapeMeasure, shape};
        }

        return carry;
      }
    ).second;
  }

  std::string name() const final {
    return "Biased CShM";
  }
};

struct CShMPathDev final : public Recognizer {
  std::vector<Shapes::Shape> validShapes;
  Eigen::MatrixXd minimumDistortionAngles;

  CShMPathDev() {
    for(const Shapes::Shape shape : Shapes::allShapes) {
      if(Shapes::size(shape) <= maxShapeSize) {
        validShapes.push_back(shape);
      }
    }

    const unsigned N = validShapes.size();
    minimumDistortionAngles.resize(N, N);
    minimumDistortionAngles.setZero();
    for(unsigned i = 0; i < N; ++i) {
      Shapes::Shape iShape = validShapes.at(i);

      for(unsigned j = i + 1; j < N; ++j) {
        Shapes::Shape jShape = validShapes.at(j);

        if(Shapes::size(iShape) == Shapes::size(jShape)) {
          minimumDistortionAngles(i, j) = Shapes::Continuous::minimumDistortionAngle(
            iShape,
            jShape
          );
        }
      }
    }

    std::cout << "valid shapes:" << Temple::stringify(
      Temple::map(validShapes, [](auto x) { return Shapes::name(x); })
    ) << "\n";
    std::cout << "minimum distortion angles:\n" << minimumDistortionAngles << "\n";
  }

  double minimumDistortionAngle(const Shapes::Shape a, const Shapes::Shape b) const {
    auto indexOfShape = [&](const Shapes::Shape shape) -> unsigned {
      auto findIter = Temple::find(validShapes, shape);
      if(findIter == std::end(validShapes)) {
        throw "Shape not found in valid shapes";
      }
      return findIter - std::begin(validShapes);
    };
    unsigned i = indexOfShape(a);
    unsigned j = indexOfShape(b);

    if(i > j) {
      std::swap(i, j);
    }

    return minimumDistortionAngles(i, j);
  }

  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    // Select shapes of matching size
    std::vector<Shapes::Shape> matchingSizeShapes;
    for(const auto shape : Shapes::allShapes) {
      if(Shapes::size(shape) == S) {
        matchingSizeShapes.push_back(shape);
      }
    }

    // Calculate continuous shape measures for all selected shapes
    Positions normalized = Shapes::Continuous::normalize(positions);
    auto shapeMeasures = Temple::map(
      matchingSizeShapes,
      [&](const Shapes::Shape shape) -> double {
        return Shapes::Continuous::shape(normalized, shape).measure;
      }
    );

    /* Calculate minimum distortion path deviations for all pairs, selecting
     * that pair for which the minimum distortion path deviation is minimal
     */
    using CarryType = std::tuple<unsigned, unsigned, double>;
    auto closestTuple = Temple::accumulate(
      Temple::Adaptors::allPairs(Temple::Adaptors::range(matchingSizeShapes.size())),
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

struct CShMPathDevBiased final : public Recognizer {
  std::vector<Shapes::Shape> validShapes;
  Eigen::MatrixXd minimumDistortionAngles;

  CShMPathDevBiased() {
    for(const Shapes::Shape shape : Shapes::allShapes) {
      if(Shapes::size(shape) <= maxShapeSize) {
        validShapes.push_back(shape);
      }
    }

    const unsigned N = validShapes.size();
    minimumDistortionAngles.resize(N, N);
    minimumDistortionAngles.setZero();
    for(unsigned i = 0; i < N; ++i) {
      Shapes::Shape iShape = validShapes.at(i);

      for(unsigned j = i + 1; j < N; ++j) {
        Shapes::Shape jShape = validShapes.at(j);

        if(Shapes::size(iShape) == Shapes::size(jShape)) {
          minimumDistortionAngles(i, j) = Shapes::Continuous::minimumDistortionAngle(
            iShape,
            jShape
          );
        }
      }
    }

    std::cout << "valid shapes:" << Temple::stringify(
      Temple::map(validShapes, [](auto x) { return Shapes::name(x); })
    ) << "\n";
    std::cout << "minimum distortion angles:\n" << minimumDistortionAngles << "\n";
  }

  double minimumDistortionAngle(const Shapes::Shape a, const Shapes::Shape b) const {
    auto indexOfShape = [&](const Shapes::Shape shape) -> unsigned {
      auto findIter = Temple::find(validShapes, shape);
      if(findIter == std::end(validShapes)) {
        throw "Shape not found in valid shapes";
      }
      return findIter - std::begin(validShapes);
    };
    unsigned i;
    unsigned j;
    std::tie(i, j) = std::minmax(
      indexOfShape(a),
      indexOfShape(b)
    );

    return minimumDistortionAngles(i, j);
  }

  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    // Select shapes of matching size
    std::vector<Shapes::Shape> matchingSizeShapes;
    for(const auto shape : Shapes::allShapes) {
      if(Shapes::size(shape) == S) {
        matchingSizeShapes.push_back(shape);
      }
    }

    // Calculate continuous shape measures for all selected shapes
    Positions normalized = Shapes::Continuous::normalize(positions);
    auto shapeMeasures = Temple::map(
      matchingSizeShapes,
      [&](const Shapes::Shape shape) -> double {
        return Shapes::Continuous::shape(normalized, shape).measure;
      }
    );

    /* Calculate minimum distortion path deviations for all pairs, selecting
     * that pair for which the minimum distortion path deviation is minimal
     */
    using CarryType = std::tuple<unsigned, unsigned, double>;
    auto closestTuple = Temple::accumulate(
      Temple::Adaptors::allPairs(Temple::Adaptors::range(matchingSizeShapes.size())),
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

    const Shapes::Shape aShape = matchingSizeShapes.at(a);
    const Shapes::Shape bShape = matchingSizeShapes.at(b);

    /* Bias tetrahedral - trigonal pyramid towards tetrahedral 80:20 */
    if(std::minmax(aShape, bShape) == std::minmax(Shapes::Shape::Tetrahedron, Shapes::Shape::TrigonalPyramid)) {
      double tetrahedronShapeMeasure;
      double trigonalPyramidShapeMeasure;
      std::tie(tetrahedronShapeMeasure, trigonalPyramidShapeMeasure) = (aShape == Shapes::Shape::Tetrahedron)
        ? std::tie(shapeMeasures.at(a), shapeMeasures.at(b))
        : std::tie(shapeMeasures.at(b), shapeMeasures.at(a));

      const double angle = std::atan2(trigonalPyramidShapeMeasure, tetrahedronShapeMeasure);
      if(angle < 2.0 * M_PI / 5) {
        return Shapes::Shape::Tetrahedron;
      }

      return Shapes::Shape::TrigonalPyramid;
    }

    /* All other cases */
    if(shapeMeasures.at(a) < shapeMeasures.at(b)) {
      return matchingSizeShapes.at(a);
    }

    return matchingSizeShapes.at(b);
  }

  std::string name() const final {
    return "CShM min dist path dev biased";
  }
};

struct ShapeDistribution final : public Recognizer {
  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    const Positions normalized = Shapes::Continuous::normalize(
      positions.rowwise().reverse()
    );
    std::vector<Shapes::Shape> viableShapes;
    for(const Shapes::Shape shape : Shapes::allShapes) {
      if(Shapes::size(shape) == S) {
        viableShapes.push_back(shape);
      }
    }
    const unsigned shapesCount = viableShapes.size();
    std::vector<Shapes::Continuous::ShapeResult> shapeMeasureResults (shapesCount);
    std::vector<boost::optional<double>> randomCloudProbabilities (shapesCount);

    for(unsigned i = 0; i < shapesCount; ++i) {
      const Shapes::Shape candidateShape = viableShapes[i];
      shapeMeasureResults[i] = Shapes::Continuous::shapeCentroidLast(
        normalized,
        candidateShape
      );
      // Shape classification for size 2 is better based on continuous shape measures themselves
      if(Shapes::size(candidateShape) > 2) {
        randomCloudProbabilities[i] = Shapes::Continuous::probabilityRandomCloud(
          shapeMeasureResults[i].measure,
          candidateShape
        );
      }
    }

    // Ensure centroids are mapped against one another
    assert(
      Temple::all_of(
        shapeMeasureResults,
        [](const auto& shapeResult) -> bool {
          assert(!shapeResult.mapping.empty());
          return shapeResult.mapping.back() == shapeResult.mapping.size() - 1;
        }
      )
    );

    // Prefer probabilities for comparison
    if(Temple::all_of(randomCloudProbabilities)) {
      const auto minElementIter = std::min_element(
        std::begin(randomCloudProbabilities),
        std::end(randomCloudProbabilities)
      );
      const unsigned minimalShapeIndex = minElementIter - std::begin(randomCloudProbabilities);
      const Shapes::Shape minimalShape = viableShapes.at(minimalShapeIndex);
      return minimalShape;
    }

    // Fall back to minimal shape measure
    const auto minElementIter = std::min_element(
      std::begin(shapeMeasureResults),
      std::end(shapeMeasureResults),
      [](const auto& a, const auto& b) -> bool {
        return a.measure < b.measure;
      }
    );
    const unsigned minimalShapeIndex = minElementIter - std::begin(shapeMeasureResults);
    const Shapes::Shape minimalShape = viableShapes.at(minimalShapeIndex);
    return minimalShape;
  }

  std::string name() const final {
    return "CShM distr";
  }
};

struct Random final : public Recognizer {
  std::reference_wrapper<PRNG> prngRef;

  Random(PRNG& globalPRNG) : prngRef(globalPRNG) {}

  Shapes::Shape identify(const Positions& positions) const final {
    const unsigned S = positions.cols() - 1;
    std::vector<Shapes::Shape> viableSymmetries;
    for(const auto& name : Shapes::allShapes) {
      if(Shapes::size(name) == S) {
        viableSymmetries.push_back(name);
      }
    }
    return Temple::Random::pick(viableSymmetries, prngRef.get());
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

using RecognizersTuple = boost::mpl::list<PureAngularDeviationSquare, AngularDeviationGeometryIndexHybrid, PureCSM, CShM, BiasedCShM, ShapeDistribution, CShMPathDev>;
constexpr std::size_t nRecognizers = boost::mpl::size<RecognizersTuple>::value;
using RecognizersList = std::vector<std::unique_ptr<Recognizer>>;

constexpr unsigned nDistortionValues = 11;
constexpr unsigned nRepeats = 100;
constexpr double maxDistortion = 1.0;

struct PythonScriptWriter {
  std::ofstream file;

  PythonScriptWriter() : file("recognition_data.py") {}

  void writeHeader(const RecognizersList& recognizers) {
    file << "n_distortion_values = " << nDistortionValues << "\n";
    file << "max_distortion = " << maxDistortion << "\n";
    file << "repeats = " << nRepeats << "\n";
    file << "shapes = [\"" << Temple::condense(
      Temple::map(Shapes::allShapes, [](auto name) { return Shapes::name(name); }),
      "\",\""
    ) << "\"]\n";
    file << "symmetrySizes = [" << Temple::condense(
      Temple::map(Shapes::allShapes, [](auto name) { return Shapes::size(name); })
    ) << "]\n";
    file << "recognizers = [\"" << Temple::condense(Temple::map(recognizers, [](const auto& f) {return f->name();}), "\",\"") << "\"]\n";
    file << "results = {}\n";
  }

  void writeSeed(int seed) {
    file << "seed = " << seed << "\n";
  }

  void addResults(const Shapes::Shape name, const std::vector<std::vector<std::vector<Shapes::Shape>>>& results) {
    const unsigned shapeIndex = Shapes::nameIndex(name);
    for(unsigned recognizerIndex = 0; recognizerIndex < nRecognizers; ++recognizerIndex) {
      const std::string array = "[" +(
        Temple::condense(
          Temple::map(
            results.at(recognizerIndex),
            [](const auto& distortionResult) {
              return Temple::condense(
                Temple::map(
                  distortionResult,
                  [](Shapes::Shape n) -> unsigned {
                    return static_cast<std::underlying_type<Shapes::Shape>::type>(n);
                  }
                )
              );
            }
          )
        )
      ) + "]";

      file << "results[(" << recognizerIndex << ", " << shapeIndex << ")] = " << array << "\n";
    }
  }
};

struct RScriptWriter {
  std::ofstream file;

  RScriptWriter() : file("recognition_data.R") {}

  void writeHeader(const RecognizersList& recognizers) {
    file << "nDistortionValues <- " << nDistortionValues << "\n";
    file << "maxDistortion <- " << maxDistortion << "\n";
    file << "nRepeats <- " << nRepeats << "\n";
    file << "shapeNames <- c(\"" << Temple::condense(
      Temple::map(Shapes::allShapes, [](auto name) { return Shapes::name(name); }),
      "\",\""
    ) << "\")\n";
    file << "symmetrySizes <- c(" << Temple::condense(
      Temple::map(Shapes::allShapes, [](auto name) { return Shapes::size(name); })
    ) << ")\n";
    file << "recognizers <- c(\"" << Temple::condense(Temple::map(recognizers, [](const auto& f) {return f->name();}), "\",\"") << "\")\n";
    file << "x.values <- c(0:" << (nDistortionValues - 1) << ") * " << maxDistortion << " / " << (nDistortionValues - 1) << "\n";
    file << "results <- array(numeric(), c(" << nRepeats << ", " << nDistortionValues << ", " << nRecognizers << ", " << Shapes::allShapes.size() << "))\n";
  }

  void writeSeed(int seed) {
    file << "seed <- " << seed << "\n";
  }

  void addResults(const Shapes::Shape name, const std::vector<std::vector<std::vector<Shapes::Shape>>>& results) {
    const unsigned shapeIndex = Shapes::nameIndex(name);
    for(unsigned recognizerIndex = 0; recognizerIndex < nRecognizers; ++recognizerIndex) {
      const std::string array = (
        "array(c("
        + Temple::condense(
          Temple::map(
            results.at(recognizerIndex),
            [](const auto& distortionResult) {
              return Temple::condense(
                Temple::map(
                  distortionResult,
                  [](Shapes::Shape n) -> unsigned {
                    return static_cast<std::underlying_type<Shapes::Shape>::type>(n);
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

      file << "results[,," << (recognizerIndex + 1) << ", " << (shapeIndex + 1) << "] <- " << array << "\n";
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
  PythonScriptWriter scriptFile {};

  std::vector<std::unique_ptr<Recognizer>> recognizerPtrs;
  boost::mpl::for_each<RecognizersTuple, boost::mpl::make_identity<boost::mpl::_1>>(
    [&recognizerPtrs](auto x) {
      using T = typename decltype(x)::type;
      recognizerPtrs.emplace_back(std::make_unique<T>());
    }
  );

  const unsigned R = recognizerPtrs.size();

  scriptFile.writeHeader(recognizerPtrs);

  Temple::JSF64 prng;
  if(options_variables_map.count("seed") > 0) {
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

  for(const Shapes::Shape name : Shapes::allShapes) {
    const unsigned S = Shapes::size(name);
    if(S > maxShapeSize) {
      // Skip anything bigger than maximum shape size for now
      break;
    }

    const unsigned shapesOfSameSize = Temple::accumulate(
      Shapes::allShapes,
      0U,
      [&S](const unsigned carry, Shapes::Shape n) -> unsigned {
        if(Shapes::size(n) == S) {
          return carry + 1;
        }

        return carry;
      }
    );

    if(shapesOfSameSize == 1) {
      continue;
    }

    Positions basePositions (3, S + 1);
    basePositions.col(0) = Eigen::Vector3d::Zero();
    basePositions.rightCols(S) = coordinates(name);

    std::vector< // each recognizer has its own vector
      std::vector< // each distortion value
        std::vector<Shapes::Shape>
      >
    > recognizerCounts(R);

    for(unsigned i = 0; i < nDistortionValues; ++i) {
      for(unsigned r = 0; r < R; ++r) {
        recognizerCounts.at(r).push_back(std::vector<Shapes::Shape> {});
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

    std::cout << "Shape: " << Shapes::name(name) << "\n";
  }
}
