/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Molecule/MolGraphWriter.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"

#include "Utils/Geometry/ElementInfo.h"

#include "Molassembler/Temple/Stringify.h"
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
Utils::ElementType MolGraphWriter::getElementType(
  const AtomIndex vertexIndex
) const {
  return graphPtr->elementType(vertexIndex);
}

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

// Vertex options
void MolGraphWriter::operator() (
  std::ostream& os,
  const AtomIndex vertexIndex
) const {
  const std::string symbolString = Utils::ElementInfo::symbol(
    graphPtr->elementType(vertexIndex)
  );

  os << "[";

  // Add element name and index label
  os << R"(label = ")" << symbolString << vertexIndex << R"(")";

  // Coloring
  // C++17 if-init
  auto bgColorFindIter = MolGraphWriter::elementBGColorMap().find(symbolString);
  if(bgColorFindIter != MolGraphWriter::elementBGColorMap().end()) {
    os << R"(, fillcolor=")" << bgColorFindIter->second << R"(")";
  } else { // default
    os << R"(, fillcolor="white")";
  }

  auto textColorFindIter = MolGraphWriter::elementTextColorMap().find(symbolString);
  if(textColorFindIter != MolGraphWriter::elementTextColorMap().end()) {
    os << R"(, fontcolor=")" << textColorFindIter->second << R"(")";
  } else { // default
    os << R"(, fontcolor="orange")";
  }

  // Font sizing
  if(symbolString == "H") {
    os << ", fontsize=10, width=.3, fixedsize=true";
  }

  // Any angles this atom is the central atom in
  std::vector<std::string> tooltipStrings {
  };

  if(stereopermutatorListPtr != nullptr) {
    if(auto stereopermutatorOption = stereopermutatorListPtr->option(vertexIndex)) {
      auto additionalTooltips = atomStereopermutatorTooltips(*stereopermutatorOption);
      std::move(
        std::begin(additionalTooltips),
        std::end(additionalTooltips),
        std::back_inserter(tooltipStrings)
      );
    }
  }

  if(!tooltipStrings.empty()) {
    os << R"(, tooltip=")"
      << Temple::condense(tooltipStrings, "&#10;"s)
      << R"(")";
  }

  os << "]";
}

// Edge options
void MolGraphWriter::operator() (
  std::ostream& os,
  const PrivateGraph::Edge& edgeIndex
) const {
  const PrivateGraph& inner = *graphPtr;
  os << "[";

  // Bond Type display options
  BondType bondType = inner.bondType(edgeIndex);
  if(bondTypeDisplayString().count(bondType) != 0) {
    os << bondTypeDisplayString().at(bondType);
  }

  os << ", penwidth=3";

  std::vector<std::string> tooltips = edgeTooltips(
    inner.source(edgeIndex),
    inner.target(edgeIndex)
  );

  const BondIndex b {inner.source(edgeIndex), inner.target(edgeIndex)};
  if(stereopermutatorListPtr != nullptr) {
    if(auto permutatorOption = stereopermutatorListPtr->option(b)) {
      if(permutatorOption->numAssignments() > 1) {
        os << R"(, color="tomato")";
      } else {
        os << R"(, color="steelblue")";
      }

      tooltips.push_back(permutatorOption->info());
    }
  }

  if(!tooltips.empty()) {
    os << R"(, edgetooltip=")"
      << Temple::condense(tooltips, "&#10;"s)
      << R"(")";
  }

  os << "]";
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
    {"C", "#33ff33"},
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

const std::map<BondType, std::string>& MolGraphWriter::bondTypeDisplayString() {
  static const std::map<BondType, std::string> map {
    {BondType::Single, R"(color = "black")"},
    {BondType::Double, R"(color = "black:invis:black")"},
    {BondType::Triple, R"(color = "black:invis:black:invis:black")"},
    {BondType::Quadruple, R"(label = "4")"},
    {BondType::Quintuple, R"(label = "5")"},
    {BondType::Sextuple, R"(label = "6")"},
    {BondType::Eta, R"(style = "dotted")"}
  };
  return map;
}

} // namespace Molassembler
} // namespace Scine
