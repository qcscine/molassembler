#ifndef INCLUDE_MOLASSEMBLER_MOL_GRAPH_WRITER_H
#define INCLUDE_MOLASSEMBLER_MOL_GRAPH_WRITER_H

/*! @file
 *
 * @brief Write graphviz representations of Molecules for state visualization
 *
 * Implements a graphviz writing helper class for the visualization of a
 * molecular graph.
 */

#include "molassembler/Graph/InnerGraph.h"

namespace molassembler {

// Helper class to write the Graph as Graphviz output
struct MolGraphWriter {
  static const std::map<std::string, std::string> elementBGColorMap;
  static const std::map<std::string, std::string> elementTextColorMap;
  static const std::map<BondType, std::string> bondTypeDisplayString;

  /* State */
  // We promise to be good and not change anything
  const InnerGraph* const graphPtr;

  explicit MolGraphWriter(const InnerGraph* passGraphPtr);

  /* Information */
  Delib::ElementType getElementType(InnerGraph::Vertex vertexIndex) const;

  // Global options
  void operator() (std::ostream& os) const;

  // Vertex options
  void operator() (std::ostream& os, InnerGraph::Vertex vertexIndex) const;

  // Edge options
  void operator() (std::ostream& os, const InnerGraph::Edge& edgeIndex) const;
};

} // namespace molassembler

#endif
