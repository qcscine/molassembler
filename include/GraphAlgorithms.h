#ifndef INCLUDE_GRAPH_ALGORITHMS_H
#define INCLUDE_GRAPH_ALGORITHMS_H

#include "AdjacencyList.h"
#include <iostream>
#include <deque> 

/* TODO
 * - Test -> We have a infinite loop of some sort, maybe an out of bounds error
 * - Cycle detection -> see Tucker, Alan. "Applied Combinatorics p.49
 * 
 * NOTES
 * - Yes, it has to be a deque. A queue, although seemingly the minimal 
 *   functionality needed here for the FIFO structure, does not offer traversal,
 *   which we need here.
 */

namespace MoleculeManip {

namespace GraphAlgorithms {

/*!
 * Connected Components algorithm. Returns a vector of unsigned numbers 
 * that maps AtomIndexType -> Connected component group ID.
 * \param adjacencies The AdjacencyList instance to process
 * \returns A vector of unsigned numbers that maps AtomIndexType -> 
 * ConnectedComponentType.
 */
std::vector<unsigned> _connected_components(const AdjacencyList& adjacencies) {
  if(adjacencies.size() == 0) return {};

  std::vector<unsigned> visited (
    adjacencies.size(),
    0
  );
  std::vector<unsigned>::iterator visIter = visited.begin();
  std::deque<AtomIndexType> toVisit = {0};
  unsigned counter = 1;
  AtomIndexType current;

  while(toVisit.size() > 0) {

    // take the first element
    current = toVisit[0];
    toVisit.pop_front();

    // mark it
    visited[current] = counter;

    /* add it's connections to toVisit if they are unvisited and not in the 
     *  deque already
     */
    for(const auto& connectedIndex : adjacencies[current]) {
      if(
        visited[connectedIndex] == 0 
        && std::find(
          toVisit.begin(),
          toVisit.end(),
          connectedIndex
        ) == toVisit.end()
      ) {
        toVisit.push_back(connectedIndex);
      }
    }

    // if toVisit is empty, increment the counter and add the next 
    // non-visited element to toVisit
    if(
      toVisit.size() == 0 &&
      visIter != visited.end()
    ) {
      // move visIter to the next unvisited element
      while(
        visIter != visited.end() 
        && *visIter != 0
      ) {
        visIter++;
      }
      // if we're not at the end, then
      if(visIter != visited.end()) {
        // increment the counter
        counter += 1;
        // and add it's index to toVisit
        toVisit.push_back(
          visIter - visited.begin()
        );
      }
    }
  }

  // all atoms must be visited
  assert(std::accumulate(
    visited.begin(),
    visited.end(),
    true,
    [](const bool& carry, const unsigned& element) {
      return (
        carry 
        && element != 0
      );
    }
  ));

  return visited;
}

/*!
 * Returns the number of connected components in an AdjacencyList.
 * \param adjacencies The AdjacencyList instance to process
 * \returns The number of connected components in the AdjacencyList.
 */
unsigned num_connected_components(const AdjacencyList& adjacencies) {
  auto visited = _connected_components(adjacencies);
  if(visited.size() == 0) return 0;
  else return *max_element(visited.begin(), visited.end());
}

/*!
 * Constructs a list of AtomIndexType groups that are connected in an 
 * AdjacencyList.
 * \param adjacencies The AdjacencyList instance to process
 * \returns An unordered vector of vectors containing atom indices that are
 * connected.
 */
std::vector<
  std::vector<AtomIndexType>
> connected_component_groups(const AdjacencyList& adjacencies) {
  auto visited = _connected_components(adjacencies);

  // guard against 0-length
  if(visited.size() == 0) return std::vector<
    std::vector<AtomIndexType>
  >();

  unsigned num_groups = *max_element(visited.begin(), visited.end());
  std::vector<
    std::vector<AtomIndexType>
  > groups (num_groups);

  // go through visited
  // it could be e.g. [1, 2, 1, 3, 4, 3]
  for(auto it = visited.begin(); it != visited.end(); it++) {
    // make space if the group number is bigger
    if(groups.size() < *it) {
      groups.resize(*it);
    }

    groups[(*it) - 1].push_back( 
      it - visited.begin() 
    );
  }

  return groups;
} 

} // eo namespace GraphAlgorithms

}

#endif
