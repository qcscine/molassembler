#ifndef INCLUDE_ADJACENCYLIST_ALGORITHMS_H
#define INCLUDE_ADJACENCYLIST_ALGORITHMS_H

#include "AdjacencyList.h"
#include "Tree.h"
#include "template_magic/templateMagic.h"

#include <deque> 

/* TODO
 * 
 * NOTES
 * - Yes, it has to be a deque. A queue, although seemingly the minimal 
 *   functionality needed here for the FIFO structure, does not offer traversal
 *   without removal, which we need here.
 *
 */

namespace MoleculeManip {

namespace AdjacencyListAlgorithms {

// WARNING: Assumes atom indices are monotonous starting from 0!
template<typename Function>
void BFSVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  Function&& function
) {
  std::vector<bool> visited (
    adjacencyList.size(), 
    false
  );
  std::deque<AtomIndexType> toVisit {initial};

  while(!TemplateMagic::all_of(visited) && toVisit.size() != 0) {
    auto current = toVisit.front();
    toVisit.pop_front();

    visited[current] = true;

    std::copy_if(
      adjacencyList[current].begin(),
      adjacencyList[current].end(),
      std::back_inserter(toVisit),
      [&visited, &toVisit](const AtomIndexType& idx) {
        return (
          !visited[idx] 
          && !TemplateMagic::makeContainsPredicate(toVisit)(idx)
        );
      }
    );

    // allow bool false return values to break
    if(!function(current)) break;
  }
}

// WARNING: Assumes atom indices are monotonous starting from 0!
template<typename Function>
void DFSVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  Function&& function
) {
  std::vector<bool> visited (
    adjacencyList.size(), 
    false
  );
  std::deque<AtomIndexType> toVisit {initial};

  while(!TemplateMagic::all_of(visited) && toVisit.size() != 0) {
    auto current = toVisit.front();
    toVisit.pop_front();

    visited[current] = true;

    std::copy_if(
      adjacencyList[current].begin(),
      adjacencyList[current].end(),
      std::front_inserter(toVisit),
      [&visited, &toVisit](const AtomIndexType& idx) {
        return (
          !visited[idx]
          && !TemplateMagic::makeContainsPredicate(toVisit)(idx)
        );
      }
    );

    // allow bool false return values to break
    if(!function(current)) break;
  }
}

/* Tree-related algorithms */
using NodeType = BasicTree::Node<AtomIndexType>;
std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies,
  const AtomIndexType& startingFrom
);

std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies
);

// Independent algorithms

/*!
 * Connected Components algorithm. Returns a vector of unsigned numbers 
 * that maps AtomIndexType -> Connected component group ID.
 * \param adjacencies The AdjacencyList instance to process
 * \returns A vector of unsigned numbers that maps AtomIndexType -> 
 * ConnectedComponentType.
 */
std::vector<unsigned> _connectedComponents(const AdjacencyList& adjacencies);

/*!
 * Returns the number of connected components in an AdjacencyList.
 * \param adjacencies The AdjacencyList instance to process
 * \returns The number of connected components in the AdjacencyList.
 */
unsigned numConnectedComponents(const AdjacencyList& adjacencies);

/*!
 * Constructs a list of AtomIndexType groups that are connected in an 
 * AdjacencyList.
 * \param adjacencies The AdjacencyList instance to process
 * \returns An unordered vector of vectors containing atom indices that are
 * connected.
 */
std::vector<
  std::vector<AtomIndexType>
> connectedComponentGroups(const AdjacencyList& adjacencies);

} // eo namespace AdjacencyListAlgorithms

}

#endif
