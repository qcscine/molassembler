/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Purely graph-based algorithms
 *
 * Contains a number of graph-level algorithms where connectivity alone is
 * relevant.
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_GRAPH_ALGORITHMS_H

#include "molassembler/Graph/PrivateGraph.h"

#include <tuple>
#include <vector>

namespace Scine {
namespace molassembler {

// Forward-declarations
class Cycles;
struct LinkInformation;
class AtomStereopermutator;

//! Core graph-level algorithms (not requiring stereopermutator information)
namespace graph_algorithms {

/*! @brief Find links between two adjacent stereopermutators, returns unordered links
 *
 * @complexity{@math{\Theta(R)} where @math{R} is the number of relevant cycles
 * containing both stereopermutators' central indices}
 */
std::vector<LinkInformation> siteLinks(
  const PrivateGraph& graph,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
);

/*! @brief Find links between substituent atoms, returns ordered links
 *
 * @complexity{@math{O(S^2)} where @math{S} is the number of substituents of the
 * central vertex}
 */
std::vector<LinkInformation> siteLinks(
  const PrivateGraph& graph,
  AtomIndex source,
  const std::vector<
    std::vector<AtomIndex>
  >& sites,
  const std::vector<AtomIndex>& excludeAdjacents
);

namespace detail {

/*! @brief Predicate to determine if a site is haptic
 *
 * @complexity{@math{\Theta(1)}}
 * @note This is not as simple as siteAtoms.size() > 1.
 */
bool isHapticSite(
  const std::vector<AtomIndex>& siteAtoms,
  const PrivateGraph& graph
);

/*! @brief Implementation of ligand site grouping algorithm
 *
 * @complexity{@math{\Theta(S)} where @math{S} is the number of substituents of
 * the central vertex}
 */
void findSites(
  const PrivateGraph& graph,
  AtomIndex centralIndex,
  const std::function<void(const std::vector<AtomIndex>&)>& callback
);

} // namespace detail

/*!
 * @brief Differentiate adjacent vertices of a central index into sites
 *
 * A site is made up of all immediately group-internally-adjacent
 * substituents of a central index. The reverse subdivision starting from a
 * ligand may be more intuitive: A ligand may be multidentate and have varying
 * hapticity at any shape position. It consists of bonding atoms (those
 * connecting to the central metal) and non-bonding atoms (which may make up a
 * linker or other extraneous groups). Bonding atoms can be subdivided into
 * connected components that are separated by non-bonding atoms, each of which
 * make up a possibly haptic group. These are called sites because
 * they each take up a site of the central index's coordination shape.
 *
 * @complexity{@math{\Theta(S)} where @math{S} is the number of substituents of
 * the central vertex}
 *
 * @post Each site's list of constituting atom indices is sorted.
 */
std::vector<
  std::vector<AtomIndex>
> ligandSiteGroups(
  const PrivateGraph& graph,
  AtomIndex centralIndex,
  const std::vector<AtomIndex>& excludeAdjacents = {}
);

/*!
 * @brief For each atom, determine sites and set eta bond for haptic sites
 *
 * Cycle through all atoms, determine the sites and set the eta bond type for
 * atoms constituting haptic ligand sites.
 *
 * @complexity{@math{\Theta(N)}}
 */
void updateEtaBonds(PrivateGraph& graph);

/**
 * @brief Calculates the graph distance from an atom to all other atoms of the
 *   molecule
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param a The atom from which to calculate distances
 * @param graph The graph
 *
 * @return A list of graph distances for each atom of the graph
 */
std::vector<unsigned> distance(AtomIndex a, const PrivateGraph& graph);

} // namespace graph_algorithms
} // namespace molassembler
} // namespace Scine

#endif
