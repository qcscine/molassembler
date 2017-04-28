#ifndef INCLUDE_BASIC_TREE_H
#define INCLUDE_BASIC_TREE_H

#include <boost/optional.hpp>
#include <memory>
#include <vector>
#include <deque>
#include <sstream>
#include <map>

namespace Tree {

template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::map<T1, T2>& map) {
  for(const auto& mapping : map) {
    os << mapping.first << " => " << mapping.second << std::endl;
  }

  return os;
}

/* TODO
 * - nice initializer? No clue what form or how to do it
 */
template<typename T>
struct Node : std::enable_shared_from_this<Node<T>> {
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
    // C++17 change to weak_from_this
    nodePtr -> parentWeakPtr = this -> shared_from_this();
    children.push_back(nodePtr);
  }

  //! Create a child with a new key
  std::shared_ptr<Node>& addChild(const T& key) {
    children.push_back(
      std::make_shared<Node>(
        key
      )
    );
    children.back() -> parentWeakPtr = this -> shared_from_this();
    return children.back();
  }

  unsigned depth() { // get depth
    unsigned depth = 0;

    std::shared_ptr<Node> iterPtr = this -> shared_from_this();

    while(
      !iterPtr->isRoot() 
    ) {
      iterPtr = (iterPtr->parentWeakPtr).lock();
      depth += 1;
    }

    return depth;
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

    while(!nodesToVisit.empty()) {
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
    return children.empty();
  }

  /*! Converts this tree into graphviz format
   * TODO assumes that an ordered pair of a node's key and it's parent's key
   * are unique in a tree, which is true for trees made from an adjacencyList,
   * but not nearly all trees.
   */
  std::string toString() const {
    std::stringstream graphViz, connections;
    std::map<
      std::pair<T, boost::optional<T>>, // own key, parent key
      std::string // ID
    > nodeIDMap;
    std::array<unsigned, 2> charIndices {97, 97}; // 'a', 'a'

    auto incrementCharIndices = [&charIndices]() {
      if(charIndices[1] == 122) {
        charIndices[0]++;
        charIndices[1] = 97;
      } else {
        charIndices[1]++;
      }
    };

    auto getNewID = [&charIndices, &incrementCharIndices]() {
      std::string ID;
      ID += char(charIndices[0]);
      ID += char(charIndices[1]);
      incrementCharIndices();
      return ID;
    };

    graphViz << "digraph tree {\n";

    forEach(
      [
        &graphViz,
        &connections,
        &nodeIDMap,
        &getNewID
      ](const Node<T>& node) -> void {
        if(auto parentPtr = node.parentWeakPtr.lock()) {
          auto newID = getNewID();
          // add pair to nodeIDMap
          nodeIDMap[
            std::make_pair(node.key, parentPtr -> key)
          ] = newID;

          // get parent ID, for that we need its parent's key too
          boost::optional<T> parentParentKey = boost::none;
          if(auto parentParentPtr = parentPtr -> parentWeakPtr.lock()) {
            parentParentKey = parentParentPtr -> key;
          }

          graphViz << "  " << newID << "[label=\"" << node.key << "\"];\n";

          connections << "  "
            << nodeIDMap.at(std::make_pair(parentPtr -> key, parentParentKey))
            << " -> " << newID <<  ";\n";
        } else { // for root pointer
          auto newID = getNewID();
          nodeIDMap[
            std::make_pair(node.key, boost::none)
          ] = newID;
          
          graphViz << "  " << newID << "[label=\"" << node.key << "\"];\n";
          // do not add a connection
        }
      }
    );

    graphViz << connections.str() << "}\n";

    return graphViz.str();
  }

  bool operator == (const Node& other) {
    if(isLeaf() && other.isLeaf()) {
      return key == other.key;
    } else if(!isLeaf() && !other.isLeaf()) {
      return (
        key == other.key
        && children.size() == other.children.size()
        && std::is_permutation(
          children.begin(),
          children.end(),
          other.children.begin(),
          [](
            const std::shared_ptr<Node<T>>& a,
            const std::shared_ptr<Node<T>>& b
          ) {
            return *a == *b;
          }
        )
      );
    } else {
      return false;
    }
  }

  bool operator != (const Node& other) {
    return !(
      *this == other
    );
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

// Tree construction helper functions
template<typename T>
std::shared_ptr<Node<T>> nodePtr(
  const T& key
) {
  return std::make_shared<Node<T>>(key);
}

template<typename T>
std::shared_ptr<Node<T>> nodePtr(
  const T& key,
  const std::initializer_list<
    std::shared_ptr<Node<T>>
  >& children
) {
  auto returnNode = nodePtr(key);

  for(const auto& child : children) {
    returnNode -> addChild(child);
  }

  return returnNode;
}

template<typename T>
std::shared_ptr<Node<T>> nodePtr(
  const T& key,
  const std::initializer_list<T>& children
) {
  auto returnNode = nodePtr(key);

  for(const auto& child : children) {
    returnNode -> addChild(child);
  }

  return returnNode;
}

} // namespace Tree

#endif
