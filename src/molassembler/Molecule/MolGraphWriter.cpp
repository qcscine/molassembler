/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Molecule/MolGraphWriter.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/StereopermutatorList.h"

#include "Utils/Geometry/ElementInfo.h"

#include "molassembler/Temple/Stringify.h"
#include "molassembler/Shapes/Data.h"

using namespace std::string_literals;

namespace Scine {
namespace molassembler {

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
    shapes::name(permutator.getShape()),
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
  auto bgColorFindIter = MolGraphWriter::elementBGColorMap.find(symbolString);
  if(bgColorFindIter != MolGraphWriter::elementBGColorMap.end()) {
    os << R"(, fillcolor=")" << bgColorFindIter->second << R"(")";
  } else { // default
    os << R"(, fillcolor="white")";
  }

  auto textColorFindIter = MolGraphWriter::elementTextColorMap.find(symbolString);
  if(textColorFindIter != MolGraphWriter::elementTextColorMap.end()) {
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
      << temple::condense(tooltipStrings, "&#10;"s)
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
  if(bondTypeDisplayString.count(bondType) != 0u) {
    os << bondTypeDisplayString.at(bondType);
  }

  os << ", penwidth=3";

  std::vector<std::string> tooltips = edgeTooltips(
    inner.source(edgeIndex),
    inner.target(edgeIndex)
  );

  const BondIndex b {inner.source(edgeIndex), inner.target(edgeIndex)};
  if(auto permutatorOption = stereopermutatorListPtr->option(b)) {
    if(permutatorOption->numAssignments() > 1) {
      os << R"(, color="tomato")";
    } else {
      os << R"(, color="steelblue")";
    }

    tooltips.push_back(permutatorOption->info());
  }

  if(!tooltips.empty()) {
    os << R"(, edgetooltip=")"
      << temple::condense(tooltips, "&#10;"s)
      << R"(")";
  }

  os << "]";
}


// Color maps
const std::map<std::string, std::string> MolGraphWriter::elementBGColorMap {
  {"H", "white"},
  {"C", "gray"},
  {"N", "blue"},
  {"O", "red"}
};

const std::map<std::string, std::string> MolGraphWriter::elementTextColorMap {
  {"H", "black"},
  {"C", "white"},
  {"N", "white"},
  {"O", "white"}
};

const std::map<BondType, std::string> MolGraphWriter::bondTypeDisplayString {
  {BondType::Single, R"(color = "black")"},
  {BondType::Double, R"(color = "black:invis:black")"},
  {BondType::Triple, R"(color = "black:invis:black:invis:black")"},
  {BondType::Quadruple, R"(label = "4")"},
  {BondType::Quintuple, R"(label = "5")"},
  {BondType::Sextuple, R"(label = "6")"},
  {BondType::Eta, R"(style = "dotted")"}
};

} // namespace molassembler
} // namespace Scine
