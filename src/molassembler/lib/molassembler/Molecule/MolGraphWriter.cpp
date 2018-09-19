// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Molecule/MolGraphWriter.h"

#include "Delib/ElementInfo.h"

namespace molassembler {

/* Constructor */
MolGraphWriter::MolGraphWriter(const InnerGraph* passGraphPtr) : graphPtr(passGraphPtr) {}

/* Information */
Delib::ElementType MolGraphWriter::getElementType(
  const AtomIndex vertexIndex
) const {
  return graphPtr->elementType(vertexIndex);
}

// Global options
void MolGraphWriter::operator() (std::ostream& os) const {
  os << "graph [fontname = \"Arial\"];\n"
    << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
    << "edge [fontname = \"Arial\"];\n";
}

// Vertex options
void MolGraphWriter::operator() (
  std::ostream& os,
  const AtomIndex vertexIndex
) const {
  const std::string symbolString = Delib::ElementInfo::symbol(
    getElementType(vertexIndex)
  );

  os << "[";

  // Add element name and index label
  os << R"(label = ")" << symbolString << vertexIndex << R"(")";

  // Coloring
  if(elementBGColorMap.count(symbolString) != 0u) {
    os << R"(, fillcolor=")" << elementBGColorMap.at(symbolString) << R"(")";
  } else { // default
    os << R"(, fillcolor="white")";
  }
  if(elementTextColorMap.count(symbolString) != 0u) {
    os << R"(, fontcolor=")" << elementTextColorMap.at(symbolString) << R"(")";
  } else { // default
    os << R"(, fontcolor="orange")";
  }

  // Font sizing
  if(symbolString == "H") {
    os << ", fontsize=10, width=.3, fixedsize=true";
  }

  os << "]";
}

// Edge options
void MolGraphWriter::operator() (
  std::ostream& os,
  const InnerGraph::Edge& edgeIndex
) const {
  const InnerGraph& inner = *graphPtr;
  os << "[";

  // Bond Type display options
  BondType bondType = inner.bondType(edgeIndex);
  if(bondTypeDisplayString.count(bondType) != 0u) {
    os << bondTypeDisplayString.at(bondType);
  }

  // If one of the bonded atoms is a hydrogen, shorten the bond
  if(
    inner.elementType(inner.target(edgeIndex)) == Delib::ElementType::H
    || inner.elementType(inner.source(edgeIndex)) == Delib::ElementType::H
  ) {
    os << ", len=0.5";
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
