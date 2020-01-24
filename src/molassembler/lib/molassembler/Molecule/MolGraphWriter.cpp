/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Molecule/MolGraphWriter.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/StereopermutatorList.h"

#include "Utils/Geometry/ElementInfo.h"

#include "temple/Stringify.h"
#include "shapes/Data.h"

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

void MolGraphWriter::writeBondStereopermutatorNodes(std::ostream& os) const {
  for(
    const auto& bondStereopermutator :
    stereopermutatorListPtr->bondStereopermutators()
  ) {
    std::string state;
    if(bondStereopermutator.assigned()) {
      state = std::to_string(bondStereopermutator.assigned().value());
    } else {
      state = "u";
    }

    state += "/"s + std::to_string(bondStereopermutator.numAssignments());

    std::string graphNodeName = "BS"
      + std::to_string(bondStereopermutator.edge().first)
      + std::to_string(bondStereopermutator.edge().second);

    std::vector<std::string> tooltipStrings {bondStereopermutator.info()};

    auto additionalTooltips = this->bondStereopermutatorTooltips(bondStereopermutator);
    std::move(
      std::begin(additionalTooltips),
      std::end(additionalTooltips),
      std::back_inserter(tooltipStrings)
    );

    os << "  " << graphNodeName << R"( [label=")" << state
      << R"(", fillcolor="steelblue", shape="box", fontcolor="white", )"
      << R"(tooltip=")"
      << temple::condense(tooltipStrings, "&#10;"s)
      << R"("];)" << "\n";

    // Add connections to the vertices (although those don't exist yet)
    os << "  " << graphNodeName << " -- " << bondStereopermutator.edge().first
      << R"( [color="gray", dir="forward", len="2"];)"
      << "\n";
    os << "  " << graphNodeName << " -- " << bondStereopermutator.edge().second
      << R"( [color="gray", dir="forward", len="2"];)"
      << "\n";
  }
}

// Global options
void MolGraphWriter::operator() (std::ostream& os) const {
  os << "graph [fontname = \"Arial\", layout=\"neato\"];\n"
    << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
    << "edge [fontname = \"Arial\"];\n";

  writeBondStereopermutatorNodes(os);
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

  if(auto stereopermutatorOption = stereopermutatorListPtr->option(vertexIndex)) {
    auto additionalTooltips = atomStereopermutatorTooltips(*stereopermutatorOption);
    std::move(
      std::begin(additionalTooltips),
      std::end(additionalTooltips),
      std::back_inserter(tooltipStrings)
    );
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

  auto tooltips = edgeTooltips(inner.source(edgeIndex), inner.target(edgeIndex));
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
