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

template<typename T>
struct Node {
/* Public members */
  boost::optional<
    std::shared_ptr<
      Node
    >
  > parentOption;
  std::vector<
    std::shared_ptr<
      Node
    >
  > children;

  T key;

/* Public member functions */
  /* Constructors */
  Node(const T& passKey) : key(passKey) {};

  //! Add a child using an existing node
  void addChild(const std::shared_ptr<Node>& nodePtr) {
    children.push_back(nodePtr);
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
    return !parentOption;
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
