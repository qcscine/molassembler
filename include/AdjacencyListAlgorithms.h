#ifndef INCLUDE_ADJACENCYLIST_ALGORITHMS_H
#define INCLUDE_ADJACENCYLIST_ALGORITHMS_H

#include "TreeAlgorithms.h"
#include <iostream>
#include <deque> 
#include <set>

/* TODO
 * - Test -> We have a infinite loop of some sort, maybe an out of bounds error
 * - Cycle detection -> see Tucker, Alan. "Applied Combinatorics" p.49
 * 
 * NOTES
 * - Yes, it has to be a deque. A queue, although seemingly the minimal 
 *   functionality needed here for the FIFO structure, does not offer traversal
 *   without removal, which we need here.
 */

namespace MoleculeManip {

namespace AdjacencyListAlgorithms {

/*!
 * Connected Components algorithm. Returns a vector of unsigned numbers 
 * that maps AtomIndexType -> Connected component group ID.
 * \param adjacencies The AdjacencyList instance to process
 * \returns A vector of unsigned numbers that maps AtomIndexType -> 
 * ConnectedComponentType.
 */
std::vector<unsigned> _connectedComponents(const AdjacencyList& adjacencies) {
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
unsigned numConnectedComponents(const AdjacencyList& adjacencies) {
  auto visited = _connectedComponents(adjacencies);
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
> connectedComponentGroups(const AdjacencyList& adjacencies) {
  auto visited = _connectedComponents(adjacencies);

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

/* Prerequisite: AdjacencyList has a single connected component. */
std::vector<
  std::vector<
    AtomIndexType
  >
> detectCycles(const AdjacencyList& adjacencies) {
  /* Rough outline of algorithm:
   *
   * Create a top-down tree from the AdjacencyList using a DFS-like algorithm.
   * Detect ring closure atoms (by detecting vertices that have been visited
   * before). The stored key in a node is the atom's index. Store node pointers
   * to repeated nodes.
   *
   * Then traverse the tree upwards from repeated nodes simultaneously until we hit 
   * a common ancestor. The common path is then the cycle.
   *
   * e.g. 
   *  
   *    0
   *   / \
   *  1 – 2
   *
   * is expanded into a tree, starting from 0, as:
   *
   * root/top ->   0 – 1 – 2 – 0   <- leaf/bottom
   *               ^           ^
   * repeated:     a           b
   *
   * then, as long as one of both node pointers has a parentOption and the sets 
   * of node keys of both pointers' traversal up the tree have an empty 
   * intersection, we move the pointers up the tree. If along the way we 
   * re-encounter the duplicate key, that pointer's traversal path is the cycle
   * set.
   *
   * step: 0 – 1 – 2 – 0    aSet = {}, bSet = {2}
   *       ^       ^
   *       a       b
   *
   * step: 0 – 1 – 2 – 0    aSet = {}, bSet = {2, 1}
   *       ^   ^
   *       a   b
   *
   * step: 0 – 1 – 2 – 0    aSet = {}, bSet = {2, 1, 0}
   *       ^^                                        ^
   *       ab                                        REPEAT FOUND
   *
   * -> Due to found repeating key, {2, 1, 0} is the cycle.
   */
  auto treeStruct = BasicTree::makeTree(adjacencies);
  std::cout << treeStruct.rootPtr << std::endl;

  std::vector<
    std::vector<
      AtomIndexType
    >
  > cycles;

  for(auto& candidateNodePtr : treeStruct.duplicateNodes) {
    AtomIndexType currentIndex = candidateNodePtr -> key;
    // get matching ptr from visited
    assert(treeStruct.nodes[candidateNodePtr -> key]);
    auto matchingPtr = treeStruct.nodes[candidateNodePtr -> key].value();

    // traverse the tree up to common ancestor

    std::pair<
      std::vector<AtomIndexType>,
      std::vector<AtomIndexType>
    > backtrackingPaths = { {candidateNodePtr -> key}, {matchingPtr -> key} };

    std::pair<
      std::set<AtomIndexType>,
      std::set<AtomIndexType>
    > backtrackingPathSets = {};

    std::vector<AtomIndexType> intersection;
    std::set_intersection(
      backtrackingPathSets.first.begin(),
      backtrackingPathSets.first.end(),
      backtrackingPathSets.second.begin(),
      backtrackingPathSets.second.end(),
      std::back_inserter(intersection)
    );

    std::experimental::optional<
      std::vector<AtomIndexType>
    > optionFoundWhileBacktracking;

    //std::cout << "Constructing path for candidate: " << candidateNodePtr -> key << std::endl;

    while(
      (
        candidateNodePtr -> parentOption
        || matchingPtr -> parentOption
      ) && intersection.size() == 0
    ) {
      /*std::cout << "Step" << std::endl
        << candidateNodePtr << std::endl
        << matchingPtr << std::endl;*/
      /* Special case:
       * Entire ring is in single branch. In this case, if one of the
       * backtracking pointers encounters the current index as a key, then that
       * backtracking chain is a ring!
       */
      // move up the tree one step for both pointers if you can
      if(candidateNodePtr -> parentOption) {
        candidateNodePtr = candidateNodePtr -> parentOption.value();
        if(candidateNodePtr -> key == currentIndex) {
          optionFoundWhileBacktracking = backtrackingPaths.first;
          break;
        }
        backtrackingPaths.first.push_back(candidateNodePtr -> key);
        backtrackingPathSets.first.insert(candidateNodePtr -> key);
      }
      if(matchingPtr -> parentOption) {
        matchingPtr = matchingPtr -> parentOption.value();
        if(matchingPtr -> key == currentIndex) {
          optionFoundWhileBacktracking = backtrackingPaths.second;
          break;
        }
        backtrackingPaths.second.push_back(matchingPtr -> key);
        backtrackingPathSets.second.insert(matchingPtr -> key);
      }

      // recompute the intersection
      intersection.clear();
      std::set_intersection(
        backtrackingPathSets.first.begin(),
        backtrackingPathSets.first.end(),
        backtrackingPathSets.second.begin(),
        backtrackingPathSets.second.end(),
        std::back_inserter(intersection)
      );
    }
    
    if(optionFoundWhileBacktracking) {
      cycles.push_back(optionFoundWhileBacktracking.value());
    } else {
      // the intersection is the "pivot" between both backtracking paths
      /*for(const auto& index: intersection) {
        std::cout << index << ", ";
      }
      std::cout << std::endl;*/
      assert(intersection.size() == 1);

      /* ring chain is then:
       * backtrackingPaths.first (forwards) up to but not including intersection[0]
       * + backtrackingPaths.second (backwards) from intersection[0] up to but not
       *    including backtrackingPaths.second[0] (this would be duplicate with
       *    backtrackingPaths.first[0], the common ring index)
       */

      std::vector<AtomIndexType> cycleChain;

      std::copy(
        backtrackingPaths.first.begin(),
        std::find( // position of the pivot
          backtrackingPaths.first.begin(),
          backtrackingPaths.first.end(),
          intersection[0]
        ),
        std::back_inserter(cycleChain)
      );

      std::copy(
        std::find( // position of pivot in reverse iteration
          backtrackingPaths.second.rbegin(),
          backtrackingPaths.second.rend(),
          intersection[0]
        ),
        backtrackingPaths.second.rend() - 1, // beginning + 1 (skip duplicate)
        std::back_inserter(cycleChain)
      );

      // add to cycles
      cycles.push_back(cycleChain);
    }
  }

  return cycles;
}

} // eo namespace AdjacencyListAlgorithms

}

#endif
