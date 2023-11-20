/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Write graphviz representations of Molecules for state visualization
 *
 * Implements a graphviz writing helper class for the visualization of a
 * molecular graph.
 */

#ifndef INCLUDE_MOLASSEMBLER_MOL_GRAPH_WRITER_H
#define INCLUDE_MOLASSEMBLER_MOL_GRAPH_WRITER_H

#include "Molassembler/Graph/PrivateGraph.h"

namespace Scine {
namespace Molassembler {

class StereopermutatorList;
class AtomStereopermutator;
class BondStereopermutator;

//! Helper class to write the Graph as Graphviz output
struct MolGraphWriter {
//!@name Static data
//!@{
  static const std::map<std::string, std::string>& elementBGColorMap();
  static const std::map<std::string, std::string>& elementTextColorMap();
//!@}

/*!@name Members
 *
 * Note: The members are pointers due to the manner in which this object is
 * passed to boost.
 *
 *!@{
 */
  //! Non-null pointer to private graph
  const PrivateGraph* const graphPtr;
  //! Maybe pointer to stereopermutator list, maybe nullptr
  const StereopermutatorList* const stereopermutatorListPtr;
//!@}

  MolGraphWriter(
    const PrivateGraph* passGraphPtr,
    const StereopermutatorList* passPermutatorListPtr
  );

  virtual ~MolGraphWriter() = default;

  /* Information */
  void operator() (std::ostream& os) const;
  void operator() (std::ostream& os, PrivateGraph::Vertex v) const;
  void operator() (std::ostream& os, const PrivateGraph::Edge& e) const;

  /* Display customization */
  //! All attributes to display for a vertex
  virtual std::map<std::string, std::string> vertexAttributes(PrivateGraph::Vertex v) const;

  //! Label string for a vertex
  virtual std::string vertexLabel(PrivateGraph::Vertex v) const;

  //! Fill and font color for a vertex
  virtual std::pair<std::string, std::string> fillFontColors(PrivateGraph::Vertex v) const;

  //! All attributes to display for an edge
  virtual std::map<std::string, std::string> edgeAttributes(const PrivateGraph::Edge& e) const;

  //! Display attribute pair for bond type (arbitrary)
  virtual std::string edgeColor(const PrivateGraph::Edge& e) const;

  //! Tooltips for an edge
  virtual std::vector<std::string> edgeTooltips(AtomIndex source, AtomIndex target) const;

  //! Tooltips for an atom stereopermutator
  virtual std::vector<std::string> atomStereopermutatorTooltips(
    const AtomStereopermutator& permutator
  ) const;

  //! Tooltips for a bond stereopermutator
  virtual std::vector<std::string> bondStereopermutatorTooltips(
    const BondStereopermutator& permutator
  ) const;
};

} // namespace Molassembler
} // namespace Scine

#endif
