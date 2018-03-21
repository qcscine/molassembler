#ifndef INCLUDE_MOLECULE_MANIP_MOL_GRAPH_WRITER_H
#define INCLUDE_MOLECULE_MANIP_MOL_GRAPH_WRITER_H

/*! @file
 *
 * Implements a graphviz writing helper class for the visualization of a
 * molecular graph.
 */

#include "common_typedefs.h"

namespace molassembler {

// Helper class to write the Graph as Graphviz output
struct MolGraphWriter {
  static const std::map<std::string, std::string> elementBGColorMap;
  static const std::map<std::string, std::string> elementTextColorMap;
  static const std::map<BondType, std::string> bondTypeDisplayString;

  /* State */
  // We promise to be good and not change anything
  const GraphType* const graphPtr;

  explicit MolGraphWriter(const GraphType* passGraphPtr);

  /* Information */
  Delib::ElementType getElementType(const AtomIndexType& vertexIndex) const;

  // Global options
  void operator() (std::ostream& os) const;

  // Vertex options
  void operator() (std::ostream& os, const AtomIndexType& vertexIndex) const;

  // Edge options
  void operator() (std::ostream& os, const EdgeIndexType& edgeIndex) const;
};

} // namespace molassembler

#endif
