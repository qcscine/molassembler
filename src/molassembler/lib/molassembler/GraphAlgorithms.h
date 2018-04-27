#ifndef INCLUDE_GRAPH_ALGORITHMS_H
#define INCLUDE_GRAPH_ALGORITHMS_H

#include "common_typedefs.h"
#include <boost/graph/breadth_first_search.hpp>

#include "Delib/ElementInfo.h"
#include "temple/Containers.h"
#include "temple/UnorderedSets.h"
#include "temple/TinySet.h"

/*! @file
 *
 * Contains a number of graph-level algorithms where connectivity alone is
 * relevant.
 */

namespace molassembler {

// Forward-declare CycleData
class CycleData;

//! Core graph-level algorithms (not requiring stereocenter information)
namespace GraphAlgorithms {

struct LinkInformation {
  //! An (asc) ordered pair of the substituent indices that are linked
  std::pair<AtomIndexType, AtomIndexType> indexPair;

  /*! The in-order atom sequence of the cycle atom indices.
   *
   * NOTE: The front and back indices are repeated.
   */
  std::vector<AtomIndexType> cycleSequence;

  /* TODO this operator isn't quite correct, in principle, the cycle sequence
   * is not the completely reduced form / has degrees of freedom, and a
   * lexicographical comparison isn't correct for it
   */
  bool operator == (const LinkInformation& other) const {
    return (
      indexPair == other.indexPair
      && cycleSequence == other.cycleSequence
    );
  }
};

std::vector<LinkInformation> substituentLinks(
  const GraphType& graph,
  const CycleData& cycleData,
  const AtomIndexType& source,
  const std::vector<AtomIndexType>& activeAdjacents
);

GraphType findAndSetEtaBonds(GraphType&& graph);

/*!
 * Returns the number of connected components of the graph. This is a central
 * property as the library enforces this number to be always one for any given
 * Molecule. The data representation of a Molecule should not contain two
 * disconnected graphs!
 */
unsigned numConnectedComponents(const GraphType& graph);

//! Data class to return removal safety information on the graph
struct RemovalSafetyData {
  std::set<EdgeIndexType> bridges;
  std::set<AtomIndexType> articulationVertices;
};

/*!
 * Determines which vertices of the graph are articulation points and which
 * edges are bridges. Removing an articulation vertex would incur the removal of
 * a bridge, whose removal always disconnects a graph. Therefore, neither can
 * be removed without prior graph changes.
 */
RemovalSafetyData getRemovalSafetyData(const GraphType& graph);

} // namespace GraphAlgorithms

} // namespace molassembler

#endif
