/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Molecule/MolGraphWriter.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"

#include "Utils/Geometry/ElementInfo.h"

#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Shapes/Data.h"

using namespace std::string_literals;

namespace Scine {
namespace Molassembler {

/* Constructor */
MolGraphWriter::MolGraphWriter(
  const PrivateGraph* passGraphPtr,
  const StereopermutatorList* passPermutatorListPtr
) : graphPtr(passGraphPtr), stereopermutatorListPtr(passPermutatorListPtr) {}

/* Information */
std::vector<std::string> MolGraphWriter::edgeTooltips(const AtomIndex /* source */, const AtomIndex /* target */) const {
  return {};
}

std::vector<std::string> MolGraphWriter::atomStereopermutatorTooltips(const AtomStereopermutator& permutator) const {
  return {
    Shapes::name(permutator.getShape()),
    permutator.info()
  };
}

std::vector<std::string> MolGraphWriter::bondStereopermutatorTooltips(const BondStereopermutator& /* permutator */) const {
  return {};
}

// Global options
void MolGraphWriter::operator() (std::ostream& os) const {
  os << "graph [fontname = \"Arial\", layout=\"neato\"];\n"
    << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
    << "edge [fontname = \"Arial\"];\n";
}

std::map<std::string, std::string> MolGraphWriter::vertexAttributes(const PrivateGraph::Vertex v) const {
  std::map<std::string, std::string> attributes;
  attributes.emplace("label", vertexLabel(v));
  auto colors = fillFontColors(v);
  attributes.emplace("fillcolor", colors.first);
  attributes.emplace("fontcolor", colors.second);

  if(stereopermutatorListPtr != nullptr) {
    if(auto stereopermutatorOption = stereopermutatorListPtr->option(v)) {
      auto tooltips = atomStereopermutatorTooltips(*stereopermutatorOption);
      if(!tooltips.empty()) {
        attributes.emplace("tooltip", Temple::condense(tooltips, "&#10;"s));
      }
    }
  }

  const Utils::ElementType e = graphPtr->elementType(v);
  if(e == Utils::ElementType::H) {
    attributes.emplace("fontsize", "10");
    attributes.emplace("width", ".3");
    attributes.emplace("fixedsize", "true");
  }
  return attributes;
}

std::pair<std::string, std::string>
MolGraphWriter::fillFontColors(const PrivateGraph::Vertex v) const {
  const Utils::ElementType e = graphPtr->elementType(v);
  const std::string symbolString = Utils::ElementInfo::symbol(e);
  return {
    elementBGColorMap().at(symbolString),
    elementTextColorMap().at(symbolString)
  };
}

std::string MolGraphWriter::vertexLabel(const PrivateGraph::Vertex v) const {
  const Utils::ElementType e = graphPtr->elementType(v);

  // Do not mark element type for hydrogen or carbon, those are commonplace
  if(e == Utils::ElementType::H || e == Utils::ElementType::C) {
    return std::to_string(v);
  }

  const std::string symbolString = Utils::ElementInfo::symbol(e);
  return symbolString + std::to_string(v);
}

// Vertex options
void MolGraphWriter::operator() (
  std::ostream& os,
  const AtomIndex v
) const {
  os << "[" << Temple::condense(
    Temple::map(
      vertexAttributes(v),
      [](const auto& strPair) {
        return strPair.first + "=\"" + strPair.second +"\"";
      }
    )
  ) << "]";
}

std::map<std::string, std::string> MolGraphWriter::edgeAttributes(const PrivateGraph::Edge& e) const {
  std::map<std::string, std::string> attributes;
  attributes.emplace("penwidth", "3");

  const PrivateGraph::Vertex i = graphPtr->source(e);
  const PrivateGraph::Vertex j = graphPtr->target(e);

  auto tooltips = edgeTooltips(i, j);

  if(stereopermutatorListPtr != nullptr) {
    if(auto permutator = stereopermutatorListPtr->option(BondIndex {i, j})) {
      tooltips.push_back(permutator->info());
    }
  }

  if(!tooltips.empty()) {
    attributes.emplace("edgetooltip", Temple::condense(tooltips, "&#10;"s));
  }

  const std::string color = edgeColor(e);
  const BondType type = graphPtr->bondType(e);
  switch(type) {
    case(BondType::Single): {
      attributes.emplace("color", color);
      break;
    }
    case(BondType::Double): {
      attributes.emplace("color", color + ":invis:" + color);
      break;
    };
    case(BondType::Triple): {
      attributes.emplace("color", color + ":invis:" + color + ":invis:" + color);
      break;
    };
    case(BondType::Quadruple): {
      attributes.emplace("color", color);
      attributes.emplace("label", "4");
      break;
    };
    case(BondType::Quintuple): {
      attributes.emplace("color", color);
      attributes.emplace("label", "5");
      break;
    };
    case(BondType::Sextuple): {
      attributes.emplace("color", color);
      attributes.emplace("label", "6");
      break;
    };
    case(BondType::Eta): {
      attributes.emplace("color", color);
      attributes.emplace("style", "dotted");
      break;
    };
  }

  return attributes;
}

std::string MolGraphWriter::edgeColor(const PrivateGraph::Edge& e) const {
  if(stereopermutatorListPtr != nullptr) {
    const BondIndex b {graphPtr->source(e), graphPtr->target(e)};
    if(auto permutatorOption = stereopermutatorListPtr->option(b)) {
      if(permutatorOption->numAssignments() > 1) {
        return "tomato";
      }

      return "steelblue";
    }
  }

  return "black";
}

// Edge options
void MolGraphWriter::operator() (
  std::ostream& os,
  const PrivateGraph::Edge& e
) const {
  os << "[" << Temple::condense(
    Temple::map(
      edgeAttributes(e),
      [](const auto& strPair) {
        return strPair.first + "=\"" + strPair.second +"\"";
      }
    )
  ) << "]";
}


// Color maps
const std::map<std::string, std::string>& MolGraphWriter::elementBGColorMap() {
  // This is similar to PyMOL element-wise coloring
  static const std::map<std::string, std::string> map = {
    {"Ac", "#6faafa"},
    {"Al", "#bfa5a5"},
    {"Am", "#545cf2"},
    {"Sb", "#9d62b5"},
    {"Ar", "#7fd0e2"},
    {"As", "#bd7fe2"},
    {"At", "#744f44"},
    {"Ba", "#00c800"},
    {"Bk", "#8a4fe2"},
    {"Be", "#c2ff00"},
    {"Bi", "#9d4fb5"},
    {"Bh", "#e00037"},
    {"B", "#ffb5b5"},
    {"Br", "#a52929"},
    {"Cd", "#ffd88f"},
    {"Ca", "#3cff00"},
    {"Cf", "#a036d3"},
    {"C", "gray"},
    {"Ce", "#ffffc7"},
    {"Cs", "#57168f"},
    {"Cl", "#1ef01e"},
    {"Cr", "#8a99c7"},
    {"Co", "#f08f9f"},
    {"Cu", "#c77f33"},
    {"Cm", "#775ce2"},
    {"D", "#e5e5e5"},
    {"Db", "#d0004f"},
    {"Dy", "#1effc7"},
    {"Es", "#b21ed3"},
    {"Er", "#00e574"},
    {"Eu", "#61ffc7"},
    {"Fm", "#b21eba"},
    {"F", "#b2ffff"},
    {"Fr", "#410066"},
    {"Gd", "#44ffc7"},
    {"Ga", "#c28f8f"},
    {"Ge", "#668f8f"},
    {"Au", "#ffd023"},
    {"Hf", "#4cc2ff"},
    {"Hs", "#e5002e"},
    {"He", "#d8ffff"},
    {"Ho", "#00ff9c"},
    {"H", "#e5e5e5"},
    {"In", "#a57472"},
    {"I", "#940094"},
    {"Ir", "#165487"},
    {"Fe", "#e06633"},
    {"Kr", "#5cb7d0"},
    {"La", "#6fd3ff"},
    {"Lr", "#c70066"},
    {"Pb", "#575961"},
    {"Li", "#cc7fff"},
    {"Lu", "#00aa24"},
    {"Mg", "#8aff00"},
    {"Mn", "#9c7ac7"},
    {"Mt", "#ea0026"},
    {"Md", "#b20ca5"},
    {"Hg", "#b7b7d0"},
    {"Mo", "#54b5b5"},
    {"Nd", "#c7ffc7"},
    {"Ne", "#b2e2f5"},
    {"Np", "#007fff"},
    {"Ni", "#4fd04f"},
    {"Nb", "#72c2c8"},
    {"N", "#3333ff"},
    {"No", "#bd0c87"},
    {"Os", "#266695"},
    {"O", "#ff4c4c"},
    {"Pd", "#006984"},
    {"P", "#ff7f00"},
    {"Pt", "#d0d0e0"},
    {"Pu", "#006aff"},
    {"Po", "#aa5c00"},
    {"K", "#8f3fd3"},
    {"Pr", "#d8ffc7"},
    {"Pm", "#a2ffc7"},
    {"Ra", "#007c00"},
    {"Rn", "#418295"},
    {"Re", "#267caa"},
    {"Rh", "#097c8c"},
    {"Rb", "#6f2eaf"},
    {"Ru", "#248f8f"},
    {"Rf", "#cc0059"},
    {"Sm", "#8fffc7"},
    {"Sc", "#e5e5e5"},
    {"Sg", "#d80044"},
    {"Se", "#ffa000"},
    {"Si", "#f0c79f"},
    {"Ag", "#bfbfbf"},
    {"Na", "#aa5cf2"},
    {"Sr", "#00ff00"},
    {"S", "#e5c53f"},
    {"Ta", "#4ca5ff"},
    {"Tc", "#3a9d9d"},
    {"Te", "#d37a00"},
    {"Tb", "#2fffc7"},
    {"Tl", "#a5544c"},
    {"Th", "#00baff"},
    {"Tm", "#00d351"},
    {"Sn", "#667f7f"},
    {"Ti", "#bfc2c7"},
    {"W", "#2194d5"},
    {"U", "#008fff"},
    {"V", "#a5a5aa"},
    {"Xe", "#419daf"},
    {"Yb", "#00bf37"},
    {"Y", "#94ffff"},
    {"Zn", "#7c7faf"},
    {"Zr", "#94e0e0"},
  };
  return map;
}

const std::map<std::string, std::string>& MolGraphWriter::elementTextColorMap() {
  static const std::map<std::string, std::string> map {
    {"Ac", "white"},
    {"Al", "white"},
    {"Am", "white"},
    {"Sb", "white"},
    {"Ar", "white"},
    {"As", "white"},
    {"At", "white"},
    {"Ba", "white"},
    {"Bk", "white"},
    {"Be", "black"},
    {"Bi", "white"},
    {"Bh", "white"},
    {"B", "black"},
    {"Br", "white"},
    {"Cd", "black"},
    {"Ca", "white"},
    {"Cf", "white"},
    {"C", "white"},
    {"Ce", "black"},
    {"Cs", "white"},
    {"Cl", "white"},
    {"Cr", "white"},
    {"Co", "white"},
    {"Cu", "white"},
    {"Cm", "white"},
    {"D", "black"},
    {"Db", "white"},
    {"Dy", "white"},
    {"Es", "white"},
    {"Er", "white"},
    {"Eu", "black"},
    {"Fm", "white"},
    {"F", "black"},
    {"Fr", "white"},
    {"Gd", "black"},
    {"Ga", "white"},
    {"Ge", "white"},
    {"Au", "black"},
    {"Hf", "white"},
    {"Hs", "white"},
    {"He", "black"},
    {"Ho", "white"},
    {"H", "black"},
    {"In", "white"},
    {"I", "white"},
    {"Ir", "white"},
    {"Fe", "white"},
    {"Kr", "white"},
    {"La", "black"},
    {"Lr", "white"},
    {"Pb", "white"},
    {"Li", "white"},
    {"Lu", "white"},
    {"Mg", "black"},
    {"Mn", "white"},
    {"Mt", "white"},
    {"Md", "white"},
    {"Hg", "white"},
    {"Mo", "white"},
    {"Nd", "black"},
    {"Ne", "black"},
    {"Np", "white"},
    {"Ni", "white"},
    {"Nb", "white"},
    {"N", "white"},
    {"No", "white"},
    {"Os", "white"},
    {"O", "white"},
    {"Pd", "white"},
    {"P", "white"},
    {"Pt", "black"},
    {"Pu", "white"},
    {"Po", "white"},
    {"K", "white"},
    {"Pr", "black"},
    {"Pm", "black"},
    {"Ra", "white"},
    {"Rn", "white"},
    {"Re", "white"},
    {"Rh", "white"},
    {"Rb", "white"},
    {"Ru", "white"},
    {"Rf", "white"},
    {"Sm", "black"},
    {"Sc", "black"},
    {"Sg", "white"},
    {"Se", "white"},
    {"Si", "black"},
    {"Ag", "black"},
    {"Na", "white"},
    {"Sr", "white"},
    {"S", "black"},
    {"Ta", "white"},
    {"Tc", "white"},
    {"Te", "white"},
    {"Tb", "black"},
    {"Tl", "white"},
    {"Th", "white"},
    {"Tm", "white"},
    {"Sn", "white"},
    {"Ti", "black"},
    {"W", "white"},
    {"U", "white"},
    {"V", "white"},
    {"Xe", "white"},
    {"Yb", "white"},
    {"Y", "black"},
    {"Zn", "white"},
    {"Zr", "black"},
  };
  return map;
}

} // namespace Molassembler
} // namespace Scine
