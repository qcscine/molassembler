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

// Forward-declare Cycles from CycleData
class Cycles;

//! Core graph-level algorithms (not requiring stereocenter information)
namespace GraphAlgorithms {

struct LinkInformation {
  //! An (asc) ordered pair of the ligand site indices that are linked
  std::pair<unsigned, unsigned> indexPair;

  /*! The in-order atom sequence of the cycle atom indices.
   *
   * NOTE: The front and back indices are repeated.
   */
  std::vector<AtomIndexType> cycleSequence;

  /*! Performs a lexicographical comparison on cycleSequence
   *
   * TODO this operator isn't quite correct, in principle, the cycle sequence
   * is not the completely reduced form / has degrees of freedom, and a
   * lexicographical comparison isn't correct for it
   */
  bool operator == (const LinkInformation& other) const;
};

std::vector<LinkInformation> substituentLinks(
  const GraphType& graph,
  const Cycles& cycleData,
  const AtomIndexType source,
  const std::vector<
    std::vector<AtomIndexType>
  >& ligands,
  const std::vector<AtomIndexType>& activeAdjacents
);

namespace detail {

//! Predicate to determin if a ligand is haptic
bool isHapticLigand(
  const std::vector<AtomIndexType>& ligand,
  const GraphType& graph
);

//! Implementation of ligand site grouping algorithm
void findLigands(
  const GraphType& graph,
  const AtomIndexType centralIndex,
  std::function<void(const std::vector<AtomIndexType>&)> callback
);

} // namespace detail

/*! Differentiate adjacent vertices of a central index into ligand site groups
 *
 * A ligand site group is made up of all immediately group-internally-adjacent
 * substituents of a central index. The reverse subdivision starting from a
 * ligand may be more intuitive: A ligand may be multidentate and have varying
 * hapticity at any denticity point. It consists of bonding atoms (those
 * connecting to the central metal) and non-bonding atoms (which may make up a
 * linker or other extraneous groups). Bonding atoms can be subdivided into
 * connected components that are separated by non-bonding atoms, each of which
 * make up a possibly haptic group. These are called ligand site groups because
 * they each take up a site of the central index's coordination geometry.
 */
std::vector<
  std::vector<AtomIndexType>
> ligandSiteGroups(const GraphType& graph, AtomIndexType centralIndex);

//! For each atom, determines ligands and sets eta bond type for haptic ligands
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
