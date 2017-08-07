#include "Symmetries.h"
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

std::vector<unsigned> iota(unsigned nValues) {
  std::vector<unsigned> values (nValues);

  std::iota(
    values.begin(),
    values.end(),
    0
  );

  return values;
}

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

std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const Symmetry::Name& symmetryName,
  unsigned rotationFunctionIndex
) {
  std::vector<unsigned> retv;

  for(
    const auto& index : 
    Symmetry::rotations(symmetryName).at(rotationFunctionIndex)
  ) {
    retv.push_back(
      indices.at(index)
    );
  }
  
  return retv;
}

Eigen::Vector3d getCoordinates(
  const Symmetry::Name& symmetryName,
  const boost::optional<unsigned>& indexInSymmetryOption
) {
  assert(
    (
      indexInSymmetryOption 
      && indexInSymmetryOption.value() < Symmetry::size(symmetryName)
    ) || !indexInSymmetryOption
  );

  if(indexInSymmetryOption) {
    return Symmetry::symmetryData.at(symmetryName).coordinates.at(
      indexInSymmetryOption.value()
    );
  }

  return {0, 0, 0};
}

double getTetrahedronVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
) {
  return (
    i - l
  ).dot(
    (j - l).cross(k - l)
  );
}

double calculateAngleDistortion(
  const Symmetry::Name& from,
  const Symmetry::Name& to,
  const std::vector<unsigned> indexMapping
) {
  assert(indexMapping.size() == Symmetry::size(to));

  double angleDistortion = 0;

  for(unsigned i = 0; i < Symmetry::size(from); ++i) {
    for(unsigned j = i + 1; j < Symmetry::size(from); ++j) {
      angleDistortion += std::fabs(
        Symmetry::angleFunction(from)(i, j)
        - Symmetry::angleFunction(to)(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angleDistortion;
}

boost::optional<unsigned> propagateIndexOptionalThroughMapping(
  const boost::optional<unsigned>& indexOptional,
  const std::vector<unsigned>& indexMapping
) {
  if(indexOptional) {
    return indexMapping.at(indexOptional.value());
  }

  return boost::none;
}


double calculateChiralDistortion(
  const Symmetry::Name& from,
  const Symmetry::Name& to,
  const std::vector<unsigned> indexMapping
) {
  assert(indexMapping.size() == Symmetry::size(to));

  double chiralDistortion = 0;

  for(const auto& tetrahedron : Symmetry::tetrahedra(from)) {
    chiralDistortion += std::fabs(
      getTetrahedronVolume(
        getCoordinates(from, tetrahedron.at(0)),
        getCoordinates(from, tetrahedron.at(1)),
        getCoordinates(from, tetrahedron.at(2)),
        getCoordinates(from, tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(0), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(1), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(2), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(3), indexMapping)
        )
      )
    );
  }

  return chiralDistortion;
}


std::set<
  std::vector<unsigned>
> generateAllRotations(
  const Symmetry::Name& symmetryName,
  const std::vector<unsigned>& indices
) {
  using IndicesList = std::vector<unsigned>;

  std::set<IndicesList> allRotations = {indices};

  unsigned linkLimit = Symmetry::rotations(symmetryName).size();

  std::vector<unsigned> chain = {0};
  std::vector<IndicesList> chainStructures = {indices};
  unsigned depth = 0;

  // begin loop
  while(chain.front() < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    auto generated = applyRotation(
      chainStructures.back(),
      symmetryName,
      chain.back()
    );

    // is it something new?
    if(allRotations.count(generated) == 0) {
      // add it to the set
      allRotations.insert(generated);

      // add it to the chain
      chainStructures.push_back(generated);
      chain.emplace_back(0);

      // increase depth, add a link
      depth++;
    } else {
      // if we are not at the maximum instruction
      if(chain.at(depth) < linkLimit - 1) {
        chain.at(depth)++;
      } else {
        // collapse the chain until we are at an incrementable position
        while(
          depth > 0
          && chain.at(depth) == linkLimit - 1
        ) {
          chain.pop_back();
          chainStructures.pop_back();
          depth--;
        }

        chain.at(depth)++;
      }
    }
  }

  return allRotations;
}

std::vector<unsigned> symPosMapping(
  const std::vector<unsigned>& mapping
) {
  std::vector<unsigned> symmetryPositions (mapping.size());

  for(unsigned i = 0; i < mapping.size(); ++i) {
    symmetryPositions.at(
      mapping.at(i)
    ) = i;
  }
  
  return symmetryPositions;
}

struct DistortionInfo {
  std::vector<unsigned> indexMapping;
  double totalDistortion;
  double chiralDistortion;

  DistortionInfo(
    const std::vector<unsigned>& passIndexMapping,
    const double& passTotalDistortion,
    const double& passChiralDistortion
  ) : indexMapping(passIndexMapping),
      totalDistortion(passTotalDistortion),
      chiralDistortion(passChiralDistortion)
  {}
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
        std::vector<DistortionInfo> distortions;

        // Create a vector of indices for the new symmetry to use as mapping
        std::vector<unsigned> indexMapping (Symmetry::size(oneLargerSymmetry));
        std::iota(
          indexMapping.begin(),
          indexMapping.end(),
          0
        );

        std::set<
          std::vector<unsigned>
        > encounteredSymmetryMappings;

        do {
          if(encounteredSymmetryMappings.count(symPosMapping(indexMapping)) == 0) {
            distortions.emplace_back(
              indexMapping,
              calculateAngleDistortion(symmetryName, oneLargerSymmetry, indexMapping),
              calculateChiralDistortion(symmetryName, oneLargerSymmetry, indexMapping)
            );

            auto allRotations = generateAllRotations(
              oneLargerSymmetry,
              symPosMapping(indexMapping)
            );

            encounteredSymmetryMappings.insert(
              allRotations.begin(),
              allRotations.end()
            );
          }
        } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

        distortionsMap[oneLargerSymmetry] = distortions;
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
