#include "Symmetries.h"
#include "Properties.h"
#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>

#include "constexpr_magic/Math.h"
#include "template_magic/Containers.h"
#include "template_magic/Numeric.h"

#include <Eigen/Geometry>

/* TODO
 * - Something still isn't right: Seesaw to TrigonalBipyramidal should be a 0!
 *   The way the source symmetry is indexed affects the total angle distortion.
 *   It's because I only fit onto the first N indices of the next symmetry,
 *   assuming the last index is the new one (this is true for sq py to oct, but
 *   not for seesaw to trig bipy).
 */

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

  RGBType operator () (const double& value) const {
    if(!(
      _min <= value
      && value <= _max
    )) {
      throw std::logic_error("Value given to gradient is not in min-max interval");
    }

    /* say I have 4-12 as min-max. value is 9
     * fractional advancement in index is (9 - 4) / (12 - 4)
     * so (value - min) / (max - min)
     *
     * red range is 0 to 255
     *
     * target is (255 - 0) * advancement + 0
     * so (to - from) * advancement + from
     *
     */


    double advancement = 0;
    if(std::fabs(_max - _min) >= 1e-10) {
      advancement = (value - _min) / (_max - _min);
    }

    /*std::cout << "advancement " << advancement << " along R [" 
      << _from.at(0) << ", " << _to.at(0) << "] is " 
      << (_from[0] + advancement * (static_cast<double>(_to[0]) - _from[0])) << std::endl;*/

    return RGBType {
      _from[0] + advancement * (static_cast<double>(_to[0]) - _from[0]),
      _from[1] + advancement * (static_cast<double>(_to[1]) - _from[1]),
      _from[2] + advancement * (static_cast<double>(_to[2]) - _from[2])
    };
  }

  std::string getHexString (const double& value) const {
    auto rgb = this -> operator() (value);

    std::stringstream stream;
    stream << "#";
    for(const auto& value : rgb) {
      stream << std::setfill('0') << std::setw(2) << std::hex << static_cast<int>(value);
    }

    return stream.str();
  }
};

std::string getGraphvizNodeName(const Symmetry::Name& symmetryName) {
  auto stringName = Symmetry::name(symmetryName);

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

void writeLigandGainPathwaysDotfile(
  const std::string& filename,
  const bool& showEdgesWithHighMultiplicity = true,
  const bool& explainTransitions = false
) {
  using namespace Symmetry::properties;

  std::ofstream dotFile(filename.c_str());

  // Write the entire graphviz file

  dotFile << "digraph g {\n";

  // Global stuff
  dotFile << R"(  graph [fontname = "Arial", nodesep="1.5", ranksep="1.2"];)" << br
    << R"(  node [fontname = "Arial", style = "filled", fillcolor="white"];)" << br 
    << R"(  edge [fontname = "Arial", penwidth=2, labelfontsize="10"];)" << br;

  std::set<Symmetry::Name> redNodes {
    Symmetry::Name::Linear,
    Symmetry::Name::Bent,
    Symmetry::Name::TrigonalPlanar
  };

  // Write clusters for every symmetry name node
  unsigned currentSize = 1;
  for(const auto& symmetryName : Symmetry::allNames) {
    if(Symmetry::size(symmetryName) > currentSize) {
      if(currentSize > 1) {
        dotFile << "  }" << br;
      }
      currentSize = Symmetry::size(symmetryName);
      dotFile << "  subgraph cluster_size" << currentSize << " {" << br;
      dotFile << R"(    color="white";)" << br;
    }

    dotFile << "    " << getGraphvizNodeName(symmetryName) << R"( [label=")"
      << Symmetry::name(symmetryName) << R"(")";
    

    if(redNodes.count(symmetryName) == 1) {
      dotFile << R"(, fillcolor="tomato", fontcolor="white")";
    } 

    dotFile << "];" << br;
  }

  // Final cluster closing bracket
  dotFile << "  }" << br << br;

  for(const auto& symmetryName : Symmetry::allNames) {
    auto size = Symmetry::size(symmetryName);

    std::map<
      Symmetry::Name,
      std::vector<DistortionInfo>
    > distortionsMap;

    for(const auto& oneLargerSymmetry : Symmetry::allNames) {
      if(Symmetry::size(oneLargerSymmetry) == size + 1) { // Ligand gain
        distortionsMap[oneLargerSymmetry] = ligandGainDistortions(
          symmetryName,
          oneLargerSymmetry
        );
      }
    }

    if(distortionsMap.size() > 0) {

      double maxDistortion = TemplateMagic::max(
        TemplateMagic::mapValues(
          distortionsMap,
          [](const auto& distortionsList) -> double {
            return TemplateMagic::min(
              TemplateMagic::map( // get distortion member of each DistortionInfo
                distortionsList,
                [](const auto& distortion) -> double {
                  return distortion.totalDistortion;
                }
              )
            );
          }
        )
      );

      RGBGradient gradient {
        {0u, 100u, 0u}, // HTML dark green
        {255u, 99u, 71u}, // HTML color tomato
        0.0, // smallest ist green
        maxDistortion // largest is red
      };

      for(const auto& symmetryDistortionsPair : distortionsMap) {
        const auto& targetSymmetry = symmetryDistortionsPair.first;
        const auto& distortionsList = symmetryDistortionsPair.second;

        // Find distortions of lowest total distortion
        double lowestDistortion = 100;

        for(const auto& distortion : distortionsList) {
          if(distortion.totalDistortion < lowestDistortion) {
            lowestDistortion = distortion.totalDistortion;
          } 
        }

        std::vector<DistortionInfo> lowestDistortions = TemplateMagic::copyIf(
          distortionsList,
          [&](const auto& distortionObject) -> bool {
            return distortionObject.totalDistortion == lowestDistortion;
          }
        );

        // Throw out all but those with the lowest chiral distortion
        double lowestChiralDistortion = 100;
        for(const auto& distortion : lowestDistortions) {
          if(distortion.chiralDistortion < lowestChiralDistortion) {
            lowestChiralDistortion = distortion.chiralDistortion;
          }
        }

        assert(!lowestDistortions.empty());

        assert( // Make sure we're not deleting everything
          !TemplateMagic::all_of(
            TemplateMagic::map(
              lowestDistortions,
              [&](const auto& distortion) -> bool {
                return(
                  distortion.chiralDistortion > lowestChiralDistortion
                );
              }
            )
          )
        );

        TemplateMagic::inplaceRemoveIf(
          lowestDistortions,
          [&](const auto& distortion) -> bool {
            return distortion.chiralDistortion > lowestChiralDistortion;
          }
        );
        
        assert(!lowestDistortions.empty());

        const unsigned multiplicity = lowestDistortions.size();

        // In case you want transitions explained
        if(explainTransitions) {
          std::cout << "Transitions of distortion " << lowestDistortion << " from "
            << Symmetry::name(symmetryName)
            << " to " << Symmetry::name(targetSymmetry) << ":\n";

          for(const auto& distortion : lowestDistortions) {
            std::cout << "mapping vector {" << TemplateMagic::condenseIterable(distortion.indexMapping) 
              << "}, chiral : " << distortion.chiralDistortion << std::endl;
          }
        }

        // Write edge
        if(
          (
            !showEdgesWithHighMultiplicity
            && multiplicity <= 3
          ) || showEdgesWithHighMultiplicity
        ) {
          dotFile << "  " << getGraphvizNodeName(symmetryName)
            << " -> " << getGraphvizNodeName(targetSymmetry)
            << " [";

          // Begin edge modifiers
          if(multiplicity <= 3) {
            std::vector<std::string> repeatColor (
              multiplicity,
              gradient.getHexString(lowestDistortion)
            );

            dotFile << "color=\"" << TemplateMagic::condenseIterable(repeatColor, ":invis:") << "\"";
          } else {
            dotFile << "color=\"" << gradient.getHexString(lowestDistortion) << "\"";
            dotFile << ", style=\"dashed\"";
          }

          dotFile << ", label=\"" << ConstexprMagic::Math::round(lowestDistortion, 2);


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

int main() {
  writeLigandGainPathwaysDotfile("ligand_gain_pathways.dot");

  return 0;
}
