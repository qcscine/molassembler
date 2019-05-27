/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Purely graph-based algorithms
 *
 * Contains a number of graph-level algorithms where connectivity alone is
 * relevant.
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_GRAPH_ALGORITHMS_H

#include "molassembler/Graph/InnerGraph.h"

#include <tuple>
#include <vector>

namespace Scine {

namespace molassembler {

// Forward-declarations
class Cycles;
struct LinkInformation;
class AtomStereopermutator;

//! Core graph-level algorithms (not requiring stereopermutator information)
namespace GraphAlgorithms {

//! Find links between two adjacent stereopermutators, returns unordered
std::vector<LinkInformation> siteLinks(
  const InnerGraph& graph,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
);

//! Find links between substituent atoms, returns ordered links
std::vector<LinkInformation> siteLinks(
  const InnerGraph& graph,
  AtomIndex source,
  const std::vector<
    std::vector<AtomIndex>
  >& sites,
  const std::vector<AtomIndex>& excludeAdjacents
);

namespace detail {

/*!
 * @brief Predicate to determine if a site is haptic
 * @note This is not as simple as siteAtoms.size() > 1.
 */
bool isHapticSite(
  const std::vector<AtomIndex>& siteAtoms,
  const InnerGraph& graph
);

//! Implementation of ligand site grouping algorithm
void findSites(
  const InnerGraph& graph,
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
 * hapticity at any symmetry position. It consists of bonding atoms (those
 * connecting to the central metal) and non-bonding atoms (which may make up a
 * linker or other extraneous groups). Bonding atoms can be subdivided into
 * connected components that are separated by non-bonding atoms, each of which
 * make up a possibly haptic group. These are called sites because
 * they each take up a site of the central index's coordination symmetry.
 *
 * @post Each site's list of constituting atom indices is sorted.
 */
std::vector<
  std::vector<AtomIndex>
> ligandSiteGroups(
  const InnerGraph& graph,
  AtomIndex centralIndex,
  const std::vector<AtomIndex>& excludeAdjacents = {}
);

/*!
 * @brief For each atom, determine sites and set eta bond for haptic sites
 *
 * Cycle through all atoms, determine the sites and set the eta bond type for
 * atoms constituting haptic ligand sites.
 */
void updateEtaBonds(InnerGraph& graph);

/**
 * @brief Calculates the graph distance from an atom to all other atoms of the
 *   molecule
 *
 * @param a The atom from which to calculate distances
 * @param graph The graph
 *
 * @return A list of graph distances for each atom of the graph
 */
std::vector<unsigned> distance(AtomIndex a, const InnerGraph& graph);

} // namespace GraphAlgorithms

} // namespace molassembler

} // namespace Scine

#endif
