/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Write graphviz representations of Molecules for state visualization
 *
 * Implements a graphviz writing helper class for the visualization of a
 * molecular graph.
 */

#ifndef INCLUDE_MOLASSEMBLER_MOL_GRAPH_WRITER_H
#define INCLUDE_MOLASSEMBLER_MOL_GRAPH_WRITER_H

#include "molassembler/Graph/InnerGraph.h"

namespace Scine {

namespace molassembler {

class StereopermutatorList;
class AtomStereopermutator;
class BondStereopermutator;

// Helper class to write the Graph as Graphviz output
struct MolGraphWriter {
//!@name Static data
//!@{
  static const std::map<std::string, std::string> elementBGColorMap;
  static const std::map<std::string, std::string> elementTextColorMap;
  static const std::map<BondType, std::string> bondTypeDisplayString;
//!@}

//!@name Closures
//!@{
  const InnerGraph* const graphPtr;
  const StereopermutatorList* const stereopermutatorListPtr;
//!@}

  MolGraphWriter(
    const InnerGraph* passGraphPtr,
    const StereopermutatorList* passPermutatorListPtr
  );

  virtual ~MolGraphWriter() = default;

  /* Information */
  Delib::ElementType getElementType(InnerGraph::Vertex vertexIndex) const;

  void writeBondStereopermutatorNodes(std::ostream& os) const;

  // Global options
  void operator() (std::ostream& os) const;

  // Vertex options
  void operator() (std::ostream& os, InnerGraph::Vertex vertexIndex) const;

  // Edge options
  void operator() (std::ostream& os, const InnerGraph::Edge& edgeIndex) const;

  virtual std::vector<std::string> edgeTooltips(AtomIndex source, AtomIndex target) const;

  virtual std::vector<std::string> atomStereopermutatorTooltips(
    const AtomStereopermutator& permutator
  ) const;

  virtual std::vector<std::string> bondStereopermutatorTooltips(
    const BondStereopermutator& permutator
  ) const;
};

} // namespace molassembler

} // namespace Scine

#endif
