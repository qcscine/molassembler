/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
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
  static const std::map<std::string, std::string> elementBGColorMap;
  static const std::map<std::string, std::string> elementTextColorMap;
  static const std::map<BondType, std::string> bondTypeDisplayString;
//!@}

//!@name Members
//!@{
  /*! Note: The members are pointers due to the manner in which this object is
   * passed to boost.
   */
  const PrivateGraph* const graphPtr;
  const StereopermutatorList* const stereopermutatorListPtr;
//!@}

  MolGraphWriter(
    const PrivateGraph* passGraphPtr,
    const StereopermutatorList* passPermutatorListPtr
  );

  virtual ~MolGraphWriter() = default;

  /* Information */
  Utils::ElementType getElementType(PrivateGraph::Vertex vertexIndex) const;

  // Global options
  void operator() (std::ostream& os) const;

  // Vertex options
  void operator() (std::ostream& os, PrivateGraph::Vertex vertexIndex) const;

  // Edge options
  void operator() (std::ostream& os, const PrivateGraph::Edge& edgeIndex) const;

  virtual std::vector<std::string> edgeTooltips(AtomIndex source, AtomIndex target) const;

  virtual std::vector<std::string> atomStereopermutatorTooltips(
    const AtomStereopermutator& permutator
  ) const;

  virtual std::vector<std::string> bondStereopermutatorTooltips(
    const BondStereopermutator& permutator
  ) const;
};

} // namespace Molassembler

} // namespace Scine

#endif
