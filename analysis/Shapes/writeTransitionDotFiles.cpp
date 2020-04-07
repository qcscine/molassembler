/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Shapes/Data.h"
#include "molassembler/Shapes/PropertyCaching.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>

#include "molassembler/Temple/GroupBy.h"
#include "molassembler/Temple/Adaptors/Transform.h"
#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Stringify.h"
#include "molassembler/Temple/constexpr/FloatingPointComparison.h"
#include "molassembler/Temple/constexpr/Numeric.h"

#include <Eigen/Geometry>

using namespace Scine;

class RgbGradient {
public:
  using RgbType = std::array<double, 3>;

private:
  const RgbType from_, to_;
  const double min_, max_;

public:
  RgbGradient(
    const RgbType& from,
    const RgbType& to,
    double min,
    double max
  ) : from_(from),
      to_(to),
      min_(min),
      max_(max) {};

  RgbType operator () (const double value) const {
    if(!(
      min_ <= value
      && value <= max_
    )) {
      throw std::logic_error("Value given to gradient is not in min-max interval");
    }

    double advancement = 0;
    if(std::fabs(max_ - min_) >= 1e-10) {
      advancement = (value - min_) / (max_ - min_);
    }

    return RgbType {{
      from_[0] + advancement * (static_cast<double>(to_[0]) - from_[0]),
      from_[1] + advancement * (static_cast<double>(to_[1]) - from_[1]),
      from_[2] + advancement * (static_cast<double>(to_[2]) - from_[2])
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

std::string getGraphvizNodeName(const shapes::Shape shape) {
  auto stringName = shapes::name(shape);

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
  using namespace shapes::properties;

  std::ofstream dotFile(filename.c_str());

  // Write the entire graphviz file

  dotFile << "digraph g {\n";

  // Global stuff
  dotFile << R"(  graph [fontname = "Arial", nodesep="1.5", ranksep="1.2"];)" << br
    << R"(  node [fontname = "Arial", style = "filled", fillcolor="white"];)" << br
    << R"(  edge [fontname = "Arial", penwidth=2, labelfontsize="10"];)" << br;

  std::set<shapes::Shape> redNodes {
    shapes::Shape::Line,
    shapes::Shape::Bent,
    shapes::Shape::EquilateralTriangle
  };

  // Write clusters for every symmetry name node
  unsigned currentSize = 1;
  for(const auto& shape : shapes::allShapes) {
    if(shapes::size(shape) > currentSize) {
      if(currentSize > 1) {
        dotFile << "  }" << br;
      }
      currentSize = shapes::size(shape);
      dotFile << "  subgraph cluster_size" << currentSize << " {" << br;
      dotFile << R"(    color="white";)" << br;
    }

    dotFile << "    " << getGraphvizNodeName(shape) << R"( [label=")"
      << shapes::name(shape) << R"(")";


    if(redNodes.count(shape) == 1) {
      dotFile << R"(, fillcolor="tomato", fontcolor="white")";
    }

    dotFile << "];" << br;
  }

  // Final cluster closing bracket
  dotFile << "  }" << br << br;

  for(const auto& sourceSymmetry : shapes::allShapes) {
    std::map<
      shapes::Shape,
      shapes::properties::ShapeTransitionGroup
    > distortionsMap;

    for(const auto& targetSymmetry : shapes::allShapes) {
      if(predicate(sourceSymmetry, targetSymmetry)) { // Ligand gain
        distortionsMap[targetSymmetry] = selectBestTransitionMappings(
          shapeTransitionMappings(
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

      RgbGradient gradient {
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
            << shapes::name(sourceSymmetry)
            << " to " << shapes::name(targetSymmetry) << ":\n";

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
            << std::round(mappingData.angularDistortion + mappingData.chiralDistortion);


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
  using namespace shapes::properties;

  std::ofstream dotFile(filename.c_str());

  // Write the entire graphviz file

  dotFile << "digraph g {\n";

  // Global stuff
  dotFile << R"(  graph [fontname = "Arial", nodesep="1.5", ranksep="1.2"];)" << br
    << R"(  node [fontname = "Arial", style = "filled", fillcolor="white"];)" << br
    << R"(  edge [fontname = "Arial", penwidth=2, labelfontsize="10"];)" << br
    << R"(  rankdir="LR";)" << br;

  std::set<shapes::Shape> redNodes {
    shapes::Shape::Line,
    shapes::Shape::Bent,
    shapes::Shape::EquilateralTriangle
  };

  // Write clusters for every symmetry name node
  unsigned currentSize = 1;
  for(const shapes::Shape shape : shapes::allShapes) {
    if(shapes::size(shape) > currentSize) {
      if(currentSize > 1) {
        dotFile << "  }" << br;
      }
      currentSize = shapes::size(shape);
      dotFile << "  subgraph cluster_size" << currentSize << " {" << br;
      dotFile << R"(    color="white";)" << br;
    }

    dotFile << "    " << getGraphvizNodeName(shape) << R"( [label=")"
      << shapes::name(shape) << R"(")";


    if(redNodes.count(shape) == 1) {
      dotFile << R"(, fillcolor="tomato", fontcolor="white")";
    }

    dotFile << "];" << br;
  }

  // Final cluster closing bracket
  dotFile << "  }" << br << br;

  for(const auto& sourceSymmetry : shapes::allShapes) {
    for(const auto& targetSymmetry : shapes::allShapes) {

      if(shapes::size(sourceSymmetry) == shapes::size(targetSymmetry) + 1) {

        std::vector<
          std::pair<
            unsigned,
            shapes::properties::ShapeTransitionGroup
          >
        > allMappings;

        for(shapes::Vertex i {0}; i < shapes::size(sourceSymmetry); ++i) {
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
                shapes::properties::floatingPointEqualityThreshold
              ) && temple::floating::isCloseRelative(
                firstMappingPair.second.chiralDistortion,
                secondMappingPair.second.chiralDistortion,
                shapes::properties::floatingPointEqualityThreshold
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
              << shapes::name(sourceSymmetry)
              << " to " << shapes::name(targetSymmetry) << ":\n";

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
              << std::round(mappingData.angularDistortion + mappingData.chiralDistortion);


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
        shapes::size(sourceSymmetry) + 1 == shapes::size(targetSymmetry)
      );
    }
  );
  writeSymmetryTransitionDotFile(
    "ligand_rearrangement_pathways.dot",
    [](const auto& sourceSymmetry, const auto& targetSymmetry) -> bool {
      return (
        shapes::size(sourceSymmetry) == shapes::size(targetSymmetry)
        && sourceSymmetry != targetSymmetry
      );
    }
  );
  writeLigandLossDotFile("ligand_loss_pathways.dot");

  return 0;
}
