#ifndef INCLUDE_TREE_ALGORITHMS_H
#define INCLUDE_TREE_ALGORITHMS_H

#include "Tree.h"

#include "template_magic/Traits.h"

/*! @file
 *
 * Defines a number of algorithms for tree traversal and information gathering.
 */

/* TODO
 * - Add variants that offer traversal up to depth limit
 */

namespace TreeAlgorithms {

template<
  typename T,
  template<
    class = std::deque<
      std::weak_ptr<
        Tree::Node<T>
      >
    >
  > class Inserter,
  typename UnaryFunction
>
void DequeVisit(
  std::shared_ptr<Tree::Node<T>>& rootPtr,
  UnaryFunction&& visitor
) {
  std::deque<
    std::weak_ptr<Tree::Node<T>>
  > toVisit {rootPtr};

  while(!toVisit.empty()) {
    if(std::shared_ptr<Tree::Node<T>> current = toVisit.front().lock()) {
      toVisit.pop_front(); // remove before copy-in!

      /* Two Inserter instantiation cases:
       * std::back_inserter -> children of current node are visited after other
       *  nodes of equal depth
       * std::front_inserter -> children of current node are visited before 
       *  other nodes of equal depth
       */

      std::copy(
        current -> children.begin(),
        current -> children.end(),
        Inserter<>(toVisit)
      );

      if(!visitor(current)) { // if the visitor returns false, stop 
        break;
      }
    } else {
      toVisit.pop_front();
    }
  }
}

template<typename T, typename UnaryFunction>
void BFSVisit(
  std::shared_ptr<Tree::Node<T>>& rootPtr,
  UnaryFunction&& visitor
) {
  DequeVisit<T, std::back_insert_iterator, UnaryFunction>(
    rootPtr,
    std::forward<UnaryFunction>(visitor)
  );
}

template<typename T, typename UnaryFunction>
void DFSVisit(
  std::shared_ptr<Tree::Node<T>>& rootPtr,
  UnaryFunction&& visitor
) {
  DequeVisit<T, std::front_insert_iterator, UnaryFunction>(
    rootPtr,
    std::forward<UnaryFunction>(visitor)
  );
}

namespace detail {

// Call visitor if Function is Binary
template<typename T, typename Function>
std::enable_if_t<
  // Traits::is_callable<Function(std::shared_ptr<Tree::Node<T>>, unsigned)>::value,
  TemplateMagic::traits::isCallableValue<
    Function(
      std::shared_ptr<
        Tree::Node<T>
      >,
      unsigned
    )
  >,
  bool
> callVisitor(
  Function&& function,
  std::shared_ptr<Tree::Node<T>>& nodePtr,
  const unsigned& depth
) {
  return function(nodePtr, depth);
}

// Call visitor if Function is Unary
template<typename T, typename Function>
std::enable_if_t<
  //Traits::is_callable<Function(std::shared_ptr<Tree::Node<T>>)>::value,
  TemplateMagic::traits::isCallableValue<
    Function(
      std::shared_ptr<
        Tree::Node<T>
      >
    )
  >,
  bool
> callVisitor(
  Function&& function,
  std::shared_ptr<Tree::Node<T>>& nodePtr,
  const unsigned& depth __attribute__ ((unused))
) {
  return function(nodePtr);
}

} // namespace detail

template<
  typename T,
  bool insertFront,
  typename UnaryOrBinaryFunction
>
void DequeVisit(
  std::shared_ptr<Tree::Node<T>>& rootPtr,
  UnaryOrBinaryFunction&& visitor,
  const unsigned& maxDepth
) {
  using PairType = std::pair<
    std::weak_ptr<Tree::Node<T>>,
    unsigned
  >;

  std::deque<PairType> toVisit {
    {rootPtr, 0}
  };

  while(!toVisit.empty()) {
    auto currentPair = toVisit.front();
    toVisit.pop_front();
    if(std::shared_ptr<Tree::Node<T>> current = currentPair.first.lock()) {

      /* Two cases:
       * insertFront -> DFS like
       * !insertFront -> BFS like
       */

      if(currentPair.second != maxDepth) {
        for(const auto& childPtr : current -> children) {
          if(insertFront) {
            toVisit.push_front(PairType({childPtr, currentPair.second + 1}));
          } else {
            toVisit.push_back(PairType({childPtr, currentPair.second + 1}));
          }
        }
      }

      if(!(
        detail::callVisitor(
          std::forward<UnaryOrBinaryFunction>(visitor),
          current,
          currentPair.second
        )
      )) { // if the visitor returns false, stop 
        break;
      }
    } else {
      toVisit.pop_front();
    }
  }
}

template<typename T, typename UnaryOrBinaryFunction>
void BFSVisit(
  std::shared_ptr<Tree::Node<T>>& rootPtr,
  UnaryOrBinaryFunction&& visitor,
  const unsigned& maxDepth
) {
  DequeVisit<T, false, UnaryOrBinaryFunction>(
    rootPtr,
    std::forward<UnaryOrBinaryFunction>(visitor),
    maxDepth
  );
}

template<typename T, typename UnaryOrBinaryFunction>
void DFSVisit(
  std::shared_ptr<Tree::Node<T>>& rootPtr,
  UnaryOrBinaryFunction&& visitor,
  const unsigned& maxDepth
) {
  DequeVisit<T, true, UnaryOrBinaryFunction>(
    rootPtr,
    std::forward<UnaryOrBinaryFunction>(visitor),
    maxDepth
  );
}

} // namespace TreeAlgorithms

#endif
