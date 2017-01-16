#ifndef INCLUDE_ADJACENCYLIST_ALGORITHMS_H
#define INCLUDE_ADJACENCYLIST_ALGORITHMS_H

#include "AdjacencyList.h"
#include "Tree.h"
#include "Traits.h"
#include "template_magic/templateMagic.h"

#include <deque> 
#include <type_traits>

/* TODO
 * - Add variant of makeTree that creates a tree of limited depth from a 
 *   starting atom
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

/*!
 * Traverses an AdjacencyList in a BFS or DFS like fashion, depending on 
 * whether Inserter is std::back_insert_iterator or std::front_insert_iterator,
 * respectively. It then calls the UnaryFunction (not template constrained) 
 * with the current atom's index.
 */
// WARNING: Assumes atom indices are monotonous starting from 0!
template<
  template<class = std::deque<AtomIndexType>
  > class Inserter,
  class UnaryFunction
>
void TraverseAdjacencyList(
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
  TraverseAdjacencyList<std::front_insert_iterator, Function>(
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
  TraverseAdjacencyList<std::back_insert_iterator, Function>(
    adjacencyList,
    initial,
    std::forward<Function>(function)
  );
}


// WARNING: Assumes atom indices are monotonous starting from 0!
/*!
 * Traverses an AdjacencyList in a BFS or DFS like fashion, depending on 
 * whether Inserter is std::back_insert_iterator or std::front_insert_iterator,
 * respectively. It then calls the UnaryOrBinaryFunction (which can be any 
 * callable with the argument types AtomIndexType, unsigned (optionally, for 
 * a binary callable) with the current node's index and the current depth from
 * the starting position specified by the parameter initial. Specifying 
 * maxDepth = 0 gives no limitation on depth (behavior is thus strictly equal
 * to traversal without depth indicators and if these are unneeded, the other
 * traversal function above will perform better; if depth indicators are 
 * desired however without depth limitation, this gives you that option).
 */
template<
  template<class = std::deque<AtomIndexType> 
  > class Inserter,
  class UnaryOrBinaryFunction
>
void TraverseAdjacencyList(
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
  TraverseAdjacencyList<std::back_insert_iterator, BinaryFunction>(
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
  TraverseAdjacencyList<std::front_insert_iterator, BinaryFunction>(
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
  const AtomIndexType& startingFrom,
  const unsigned& maxDepth
);

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
