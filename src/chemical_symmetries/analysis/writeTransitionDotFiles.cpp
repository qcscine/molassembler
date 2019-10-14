/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "chemical_symmetries/Symmetries.h"
#include "chemical_symmetries/DynamicProperties.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>

#include "temple/GroupBy.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Stringify.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"

#include <Eigen/Geometry>

using namespace Scine;

class RGBGradient {
public:
  using RGBType = std::array<double, 3>;

private:
  const RGBType _from, _to;
  const double _min, _max;

public:
  RGBGradient(
    const RGBType& from,
    const RGBType& to,
    double min,
    double max
  ) : _from(from),
      _to(to),
      _min(min),
      _max(max) {};

  RGBType operator () (const double value) const {
    if(!(
      _min <= value
      && value <= _max
    )) {
      throw std::logic_error("Value given to gradient is not in min-max interval");
    }

    double advancement = 0;
    if(std::fabs(_max - _min) >= 1e-10) {
      advancement = (value - _min) / (_max - _min);
    }

    return RGBType {{
      _from[0] + advancement * (static_cast<double>(_to[0]) - _from[0]),
      _from[1] + advancement * (static_cast<double>(_to[1]) - _from[1]),
      _from[2] + advancement * (static_cast<double>(_to[2]) - _from[2])
    }};
  }

  std::string getHexString (const double value) const {
    auto rgb = this -> operator() (value);

    std::stringstream stream;
    stream << "#";
    for(const auto colorValue : rgb) {
      stream << std::setfill('0') << std::setw(2) << std::hex
        << static_cast<int>(colorValue);
    }

    return stream.str();
  }
};

std::string getGraphvizNodeName(const Symmetry::Shape shape) {
  auto stringName = Symmetry::name(shape);

  stringName.erase(
    std::remove_if(
      stringName.begin(),
      stringName.end(),
      [](const char& singleChar) -> bool {
        return (
          singleChar == ' '
          || singleChar == '-'
        );
      }
    ),
    stringName.end()
  );
  return stringName;
}

constexpr char br = '\n';

template<class SelectionPredicateFunction>
void writeSymmetryTransitionDotFile(
  const std::string& filename,
  SelectionPredicateFunction&& predicate,
  const bool showEdgesWithHighMultiplicity = true,
  const bool explainTransitions = false
) {
  using namespace Symmetry::properties;

  std::ofstream dotFile(filename.c_str());

  // Write the entire graphviz file

  dotFile << "digraph g {\n";

  // Global stuff
  dotFile << R"(  graph [fontname = "Arial", nodesep="1.5", ranksep="1.2"];)" << br
    << R"(  node [fontname = "Arial", style = "filled", fillcolor="white"];)" << br
    << R"(  edge [fontname = "Arial", penwidth=2, labelfontsize="10"];)" << br;

  std::set<Symmetry::Shape> redNodes {
    Symmetry::Shape::Line,
    Symmetry::Shape::Bent,
    Symmetry::Shape::EquilateralTriangle
  };

  // Write clusters for every symmetry name node
  unsigned currentSize = 1;
  for(const auto& shape : Symmetry::allShapes) {
    if(Symmetry::size(shape) > currentSize) {
      if(currentSize > 1) {
        dotFile << "  }" << br;
      }
      currentSize = Symmetry::size(shape);
      dotFile << "  subgraph cluster_size" << currentSize << " {" << br;
      dotFile << R"(    color="white";)" << br;
    }

    dotFile << "    " << getGraphvizNodeName(shape) << R"( [label=")"
      << Symmetry::name(shape) << R"(")";


    if(redNodes.count(shape) == 1) {
      dotFile << R"(, fillcolor="tomato", fontcolor="white")";
    }

    dotFile << "];" << br;
  }

  // Final cluster closing bracket
  dotFile << "  }" << br << br;

  for(const auto& sourceSymmetry : Symmetry::allShapes) {
    std::map<
      Symmetry::Shape,
      Symmetry::properties::SymmetryTransitionGroup
    > distortionsMap;

    for(const auto& targetSymmetry : Symmetry::allShapes) {
      if(predicate(sourceSymmetry, targetSymmetry)) { // Ligand gain
        distortionsMap[targetSymmetry] = selectBestTransitionMappings(
          symmetryTransitionMappings(
            sourceSymmetry,
            targetSymmetry
          )
        );
      }
    }

    if(!distortionsMap.empty()) {
      double maxDistortion = temple::max(
        temple::adaptors::transform(
          distortionsMap,
          [](const auto& mapIteratorPair) -> double {
            const auto& ligandGainReturnStruct = mapIteratorPair.second;

            return (
              ligandGainReturnStruct.angularDistortion
              + ligandGainReturnStruct.chiralDistortion
            );
          }
        )
      );

      RGBGradient gradient {
        {{0u, 100u, 0u}}, // HTML dark green
        {{255u, 99u, 71u}}, // HTML color tomato
        0.0, // smallest ist green
        maxDistortion // largest is red
      };

      for(const auto& symmetryDistortionsPair : distortionsMap) {
        const auto& targetSymmetry = symmetryDistortionsPair.first;
        const auto& mappingData = symmetryDistortionsPair.second;

        const unsigned multiplicity = mappingData.indexMappings.size();

        // In case you want transitions explained
        if(explainTransitions) {
          std::cout << "Transitions of distortion "
            << (
              mappingData.angularDistortion
              + mappingData.chiralDistortion
            ) << " from "
            << Symmetry::name(sourceSymmetry)
            << " to " << Symmetry::name(targetSymmetry) << ":\n";

          for(const auto& mapping : mappingData.indexMappings) {
            std::cout << "mapping {"
              << temple::condense(mapping)
              << "}" << std::endl;
          }
        }

        // Write edge
        if(
          (
            !showEdgesWithHighMultiplicity
            && multiplicity <= 3
          ) || showEdgesWithHighMultiplicity
        ) {
          dotFile << "  " << getGraphvizNodeName(sourceSymmetry)
            << " -> " << getGraphvizNodeName(targetSymmetry)
            << " [";

          // Begin edge modifiers
          if(multiplicity <= 3) {
            std::vector<std::string> repeatColor (
              multiplicity,
              gradient.getHexString(mappingData.angularDistortion + mappingData.chiralDistortion)
            );

            dotFile << "color=\""
              << temple::condense(repeatColor, ":invis:")
              << "\"";
          } else {
            dotFile << "color=\""
              << gradient.getHexString(mappingData.angularDistortion + mappingData.chiralDistortion)
              << R"(", style="dashed")";
          }

          dotFile << ", label=\""
            << temple::Math::round(
              mappingData.angularDistortion + mappingData.chiralDistortion,
              2
            );


          if(multiplicity > 3) {
            dotFile << " (" << multiplicity << ")";
          }

          // close label
          dotFile << "\"";

          // End edge modifiers
          dotFile << "];\n";
        }
      }

    }

  }

  // Final graph closing bracket
  dotFile << "}" << br;
}

void writeLigandLossDotFile(
  const std::string& filename,
  const bool showEdgesWithHighMultiplicity = true,
  const bool explainTransitions = false
) {
  using namespace Symmetry::properties;

  std::ofstream dotFile(filename.c_str());

  // Write the entire graphviz file

  dotFile << "digraph g {\n";

  // Global stuff
  dotFile << R"(  graph [fontname = "Arial", nodesep="1.5", ranksep="1.2"];)" << br
    << R"(  node [fontname = "Arial", style = "filled", fillcolor="white"];)" << br
    << R"(  edge [fontname = "Arial", penwidth=2, labelfontsize="10"];)" << br
    << R"(  rankdir="LR";)" << br;

  std::set<Symmetry::Shape> redNodes {
    Symmetry::Shape::Line,
    Symmetry::Shape::Bent,
    Symmetry::Shape::EquilateralTriangle
  };

  // Write clusters for every symmetry name node
  unsigned currentSize = 1;
  for(const Symmetry::Shape shape : Symmetry::allShapes) {
    if(Symmetry::size(shape) > currentSize) {
      if(currentSize > 1) {
        dotFile << "  }" << br;
      }
      currentSize = Symmetry::size(shape);
      dotFile << "  subgraph cluster_size" << currentSize << " {" << br;
      dotFile << R"(    color="white";)" << br;
    }

    dotFile << "    " << getGraphvizNodeName(shape) << R"( [label=")"
      << Symmetry::name(shape) << R"(")";


    if(redNodes.count(shape) == 1) {
      dotFile << R"(, fillcolor="tomato", fontcolor="white")";
    }

    dotFile << "];" << br;
  }

  // Final cluster closing bracket
  dotFile << "  }" << br << br;

  for(const auto& sourceSymmetry : Symmetry::allShapes) {
    for(const auto& targetSymmetry : Symmetry::allShapes) {

      if(Symmetry::size(sourceSymmetry) == Symmetry::size(targetSymmetry) + 1) {

        std::vector<
          std::pair<
            unsigned,
            Symmetry::properties::SymmetryTransitionGroup
          >
        > allMappings;

        for(unsigned i = 0; i < Symmetry::size(sourceSymmetry); ++i) {
          allMappings.emplace_back(
            i,
            selectBestTransitionMappings(
              ligandLossTransitionMappings(
                sourceSymmetry,
                targetSymmetry,
                i
              )
            )
          );
        }

        // Analyze all mappings - which indices have "identical" target mappings?
        auto groups = temple::groupByEquality(
          allMappings,
          [&](const auto& firstMappingPair, const auto& secondMappingPair) -> bool {
            return (
              temple::floating::isCloseRelative(
                firstMappingPair.second.angularDistortion,
                secondMappingPair.second.angularDistortion,
                Symmetry::properties::floatingPointEqualityThreshold
              ) && temple::floating::isCloseRelative(
                firstMappingPair.second.chiralDistortion,
                secondMappingPair.second.chiralDistortion,
                Symmetry::properties::floatingPointEqualityThreshold
              )
            );
          }
        );

        for(const auto& mappingGroup : groups) {
          const auto equivalentPositions = temple::map(
            mappingGroup,
            [](const auto& mappingPair) -> unsigned {
              return mappingPair.first;
            }
          );

          const auto& mappingData = mappingGroup.front().second;

          const unsigned multiplicity = mappingData.indexMappings.size();

          // In case you want transitions explained
          if(explainTransitions) {
            std::cout << "Transitions of distortion "
              << (
                mappingData.angularDistortion
                + mappingData.chiralDistortion
              ) << " from "
              << Symmetry::name(sourceSymmetry)
              << " to " << Symmetry::name(targetSymmetry) << ":\n";

            for(const auto& mapping : mappingData.indexMappings) {
              std::cout << "mapping {"
                << temple::condense(mapping)
                << "}" << std::endl;
            }
          }

          // Write edge
          if(
            (
              !showEdgesWithHighMultiplicity
              && multiplicity <= 3
            ) || showEdgesWithHighMultiplicity
          ) {
            dotFile << "  " << getGraphvizNodeName(sourceSymmetry)
              << " -> " << getGraphvizNodeName(targetSymmetry)
              << " [";

            // Begin edge modifiers
            if(multiplicity <= 3) {
              std::vector<std::string> repeatColor (
                multiplicity,
                "black"
              );

              dotFile << "color=\""
                << temple::condense(repeatColor, ":invis:")
                << "\"";
            } else {
              dotFile << R"(color="black", style="dashed")";
            }

            dotFile << ", label=\""
              << temple::Math::round(
                mappingData.angularDistortion + mappingData.chiralDistortion,
                2
              );


            if(multiplicity > 3) {
              dotFile << " (" << multiplicity << ")";
            }

            // Add equivalent positions to label
            dotFile << " {"
              << temple::condense(equivalentPositions)
              << "}";

            // close label
            dotFile << "\"";

            // End edge modifiers
            dotFile << "];\n";
          }
        }
      }
    }
  }

  // Final graph closing bracket
  dotFile << "}" << br;
}

int main() {
  writeSymmetryTransitionDotFile(
    "ligand_gain_pathways.dot",
    [](const auto& sourceSymmetry, const auto& targetSymmetry) -> bool {
      return (
        Symmetry::size(sourceSymmetry) + 1 == Symmetry::size(targetSymmetry)
      );
    }
  );
  writeSymmetryTransitionDotFile(
    "ligand_rearrangement_pathways.dot",
    [](const auto& sourceSymmetry, const auto& targetSymmetry) -> bool {
      return (
        Symmetry::size(sourceSymmetry) == Symmetry::size(targetSymmetry)
        && sourceSymmetry != targetSymmetry
      );
    }
  );
  writeLigandLossDotFile("ligand_loss_pathways.dot");

  return 0;
}
