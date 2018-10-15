// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_GRAPH_ALGORITHMS_H

#include "molassembler/Graph/InnerGraph.h"

#include <tuple>
#include <vector>

/*! @file
 *
 * @brief Purely graph-based algorithms
 *
 * Contains a number of graph-level algorithms where connectivity alone is
 * relevant.
 */

namespace molassembler {

// Forward-declarations
class Cycles;
struct LinkInformation;

//! Core graph-level algorithms (not requiring stereopermutator information)
namespace GraphAlgorithms {

std::vector<LinkInformation> substituentLinks(
  const InnerGraph& graph,
  const Cycles& cycleData,
  AtomIndex source,
  const std::vector<
    std::vector<AtomIndex>
  >& ligands,
  const std::vector<AtomIndex>& excludeAdjacents
);

namespace detail {

//! Predicate to determin if a ligand is haptic
bool isHapticLigand(
  const std::vector<AtomIndex>& ligand,
  const InnerGraph& graph
);

//! Implementation of ligand site grouping algorithm
void findLigands(
  const InnerGraph& graph,
  AtomIndex centralIndex,
  const std::function<void(const std::vector<AtomIndex>&)>& callback
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
  std::vector<AtomIndex>
> ligandSiteGroups(
  const InnerGraph& graph,
  AtomIndex centralIndex,
  const std::vector<AtomIndex>& excludeAdjacents = {}
);

//! For each atom, determines ligands and sets eta bond type for haptic ligands
void findAndSetEtaBonds(InnerGraph& graph);

/*!
 * Returns the number of connected components of the graph. This is a central
 * property as the library enforces this number to be always one for any given
 * Molecule. The data representation of a Molecule should not contain two
 * disconnected graphs!
 */
[[deprecated]]
unsigned numConnectedComponents(const InnerGraph& graph);

std::vector<unsigned> distance(AtomIndex a, const InnerGraph& graph);

} // namespace GraphAlgorithms

} // namespace molassembler

#endif
