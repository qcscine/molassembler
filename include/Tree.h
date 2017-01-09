#ifndef INCLUDE_BASIC_TREE_H
#define INCLUDE_BASIC_TREE_H

#include <boost/optional.hpp>
#include <memory>
#include <vector>
#include <iostream>
#include <deque>
#include <sstream>
#include <map>


namespace BasicTree {

template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::map<T1, T2>& map) {
  for(const auto& mapping : map) {
    os << mapping.first << " => " << mapping.second << std::endl;
  }

  return os;
}

/* TODO
 * - boost::optional<
 *     std::weak_ptr<
 *       Node
 *     >
 *   > parentOption is overkill, a simple weak_ptr will do just fine, but will
 *   alter all dependent code.
 */
template<typename T>
struct Node {
/* Public members */
  std::weak_ptr<
    Node
  > parentWeakPtr;

  std::vector<
    std::shared_ptr<
      Node
    >
  > children;

  T key;

/* Public member functions */
  /* Constructors */
  Node(const T& passKey) : key(passKey) {}
  Node(
    std::shared_ptr<
      Node<T>
    >& parentPtr,
    const T& passKey
  ) : parentWeakPtr(parentPtr), key(passKey) {}

  //! Add a child using an existing node
  void addChild(const std::shared_ptr<Node>& nodePtr) {
    children.push_back(nodePtr);
    // cannot notify child that we are parent since this object cannot
    // reference a shared_ptr to itself
  }

  //! Create a child with a new key
  std::shared_ptr<Node>& addChild(const T& key) {
    children.push_back(
      std::make_shared<Node>(
        key
      )
    );
    return children.back();
  }


  /* Information */
  //! Execute a lambda for every Node in the tree
  void forEach(
    std::function<
      void(const Node<T>&)
    > function
  ) const {
    function(*this);

    std::deque<
      std::shared_ptr<
        Node<T>
      >
    > nodesToVisit;

    std::copy(
      children.begin(),
      children.end(),
      std::back_inserter(nodesToVisit)
    );

    while(nodesToVisit.size() != 0) {
      auto current = nodesToVisit.front();
      nodesToVisit.pop_front();

      std::copy(
        current -> children.begin(),
        current -> children.end(),
        std::back_inserter(nodesToVisit)
      );

      function(*current);
    }
  }
  bool isRoot() const {
    /* if the weak pointer is "expired", i.e. use_count is zero, then this 
     * Node has no parent. This should, given algorithm correctness, only be 
     * the case where this is the root pointer.
     *
     * During destruction, child nodes are destroyed first, since they do not 
     * own their parents through shared_ptr. Nodes with children are destroyed
     * as soon as their children have been destroyed.
     */
    return parentWeakPtr.expired(); 
  }

  bool isLeaf() const {
    return children.size() == 0;
  }

  //! Converts this tree into graphviz format
  std::string toString() const {
    std::stringstream graphViz;
    graphViz << "digraph tree {\n";

    forEach(
      [&graphViz](const Node<T>& node) -> void {
        for(const auto& childPtr : node.children) {
          graphViz << "  \"" << node.key << "\" -> \"" << childPtr -> key 
            << "\";\n";
        }
      }
    );

    graphViz << "}\n";

    return graphViz.str();
  }
};

template<typename T>
std::ostream& operator << (
    std::ostream& os,
    const std::shared_ptr<
      Node<T> 
    >& rootPtr
) {
  os << rootPtr -> toString();
  return os;
}

} // eo namespace

#endif
