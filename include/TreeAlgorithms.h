#ifndef INCLUDE_TREE_ALGORITHMS_H
#define INCLUDE_TREE_ALGORITHMS_H

#include "Tree.h"

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

  while(toVisit.size() > 0) {
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


}

#endif
