#ifndef INCLUDE_MOLASSEMBLER_SUBGRAPHS_H
#define INCLUDE_MOLASSEMBLER_SUBGRAPHS_H

#include "boost/graph/mcgregor_common_subgraphs.hpp"

#include "Molecule.h"

namespace molassembler {

namespace subgraphs {

//! Maybe consider using boost::bimap for this instead?
using IndexMap = std::unordered_map<AtomIndexType, AtomIndexType>;

/*! Mappings for the maximum common subgraph between two molecules
 *
 * Finds an index mapping from a to b representing the maximum common subgraph
 * (if present).
 *
 * \warning For subgraph comparison, only element and bond types are considered
 * Stereocenters and Stereopermutations are not graph-local properties suitable
 * to a common substructure matching.
 */
std::vector<IndexMap> maximum(
  const GraphType& a,
  const GraphType& b
);

//! Unique mappings for the maximum common subgraph between two molecules
std::vector<IndexMap> maximumUnique(
  const GraphType& a,
  const GraphType& b
);

} // namespace subgraphs

} // namespace molassembler

#endif
