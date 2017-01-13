#ifndef INCLUDE_ADJACENCYLIST_ALGORITHMS_H
#define INCLUDE_ADJACENCYLIST_ALGORITHMS_H

#include "AdjacencyList.h"
#include "Tree.h"
#include "Traits.h"
#include "template_magic/templateMagic.h"

#include <deque> 
#include <type_traits>

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

namespace detail {

// Call visitor if Function is Binary
template<typename Function>
std::enable_if_t<
  Traits::is_callable<Function(AtomIndexType, unsigned)>::value,
  bool
> callVisitor(
  Function&& function,
  const AtomIndexType& current,
  const unsigned& depth
) {
  return function(current, depth);
}

// Call visitor if Function is Unary
template<typename Function>
std::enable_if_t<
  Traits::is_callable<Function(AtomIndexType)>::value,
  bool
> callVisitor(
  Function&& function,
  const AtomIndexType& current,
  const unsigned& depth __attribute__ ((unused))
) {
  return function(current);
}

} // eo namespace detail

// WARNING: Assumes atom indices are monotonous starting from 0!
template<
  template<class = std::deque<AtomIndexType>
  > class Inserter,
  class UnaryFunction
>
void DequeVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  UnaryFunction&& function
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
      Inserter<>(toVisit),
      [&](const AtomIndexType& idx) {
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

template<typename Function>
void DFSVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  Function&& function
) {
  DequeVisit<std::front_insert_iterator, Function>(
    adjacencyList,
    initial,
    std::forward<Function>(function)
  );
}

template<typename Function>
void BFSVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  Function&& function
) {
  DequeVisit<std::back_insert_iterator, Function>(
    adjacencyList,
    initial,
    std::forward<Function>(function)
  );
}


// WARNING: Assumes atom indices are monotonous starting from 0!
template<
  template<class = std::deque<AtomIndexType> 
  > class Inserter,
  class UnaryOrBinaryFunction
>
void DequeVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  UnaryOrBinaryFunction&& function,
  const unsigned& maxDepth
) {
  std::vector<bool> visited (
    adjacencyList.size(), 
    false
  );
  std::deque<AtomIndexType> toVisit {initial};

  /* keys in an adjacencyList are unique, so we can directly map them to their
   * depth from the starting point
   */
  std::map<AtomIndexType, unsigned> depthMap {
    {initial, 0}
  };

  while(!TemplateMagic::all_of(visited) && toVisit.size() != 0) {
    auto current = toVisit.front();
    toVisit.pop_front();

    visited[current] = true;

    std::copy_if(
      adjacencyList[current].begin(),
      adjacencyList[current].end(),
      Inserter<>(toVisit),
      [&](const AtomIndexType& idx) {
        bool toCopy = (
          !visited[idx] 
          && !TemplateMagic::makeContainsPredicate(toVisit)(idx)
          && depthMap[current] + 1 != maxDepth
        );

        if(toCopy) {
          depthMap[idx] = depthMap[current] + 1;
        }

        return toCopy;
      }
    );

    // allow bool false return values to break
    // if(!function(current, depthMap[current])) break;
    /* - allow bool false return values to break
     * - use callVisitor to distinguish between unary and binary function
     *   instantiation
     */
    if(!(
      detail::callVisitor(
        std::forward<UnaryOrBinaryFunction>(function),
        current,
        depthMap[current]
      )
    )) break;
  }
}


template<typename BinaryFunction>
void BFSVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  BinaryFunction&& function,
  const unsigned& maxDepth
) {
  DequeVisit<std::back_insert_iterator, BinaryFunction>(
    adjacencyList,
    initial,
    std::forward<BinaryFunction>(function),
    maxDepth
  );
}

template<typename BinaryFunction>
void DFSVisit(
  const AdjacencyList& adjacencyList,
  const AtomIndexType& initial,
  BinaryFunction&& function,
  const unsigned& maxDepth
) {
  DequeVisit<std::front_insert_iterator, BinaryFunction>(
    adjacencyList,
    initial,
    std::forward<BinaryFunction>(function),
    maxDepth
  );
}

/* Tree-related algorithms */
using NodeType = Tree::Node<AtomIndexType>;
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
