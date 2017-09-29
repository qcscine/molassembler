#ifndef INCLUDE_CONSTEXPR_MAGIC_BTREE_H
#define INCLUDE_CONSTEXPR_MAGIC_BTREE_H

#include "DynamicArray.h"
#include "Math.h"
#include "Optional.h"

#include <sstream>
#include <map>
#include <set>

/*! @file
 *
 * Implements a BTree which doesn't store key-value pairs, only keys.
 */

namespace ConstexprMagic {

namespace BTreeProperties {

/*!
 * Calculates the maximum number of nodes in a B-tree of a specific minDegree and
 * height.
 */
constexpr size_t maxNodesInTree(const size_t& height, const size_t& minDegree) {
  return static_cast<double>(
    Math::pow(2 * minDegree, static_cast<unsigned>(height + 1)) - 1
  ) / (
    2 * minDegree - 1
  );
}

/*!
 * Calculates the minimal height needed for a B-Tree of a specific minDegree to
 * be able to hold a certain number of keys
 */
constexpr size_t minHeight(const size_t& numKeys, const size_t& minDegree) {
  return Math::ceil(
    Math::log(
      static_cast<double>(numKeys + 1),
      static_cast<double>(2 * minDegree)
    ) - 1
  );
}

} // namespace BTreeImpl

/*! A B-Tree that merely stores keys, not key-value pairs.
 *
 * This class is a B-Tree that stores keys in its nodes, not key-value pairs. 
 * (This can be turned into an associative container with a simple modification,
 * though: Set KeyType to std::pair<KeyType, ValueType> and supply custom
 * LessThanComparator and EqualityComparators that merely compare using the
 * KeyType. Then return the ValueType from a lookup which has a default
 * ValueType.)
 *
 * @tparam KeyType - The desired stored type
 * @tparam minDegree - The minimum degree t of the B-Tree. Must be >= 2, since
 *   Nodes store a minimum of t-1 keys. A Node of degree t can store up to 2t-1 
 *   keys has up to 2t children.
 * @tparam numElements - The intended maximimum amount of stored elements. The 
 *   class must allocate all nodes that may possibly be stored at any time, so
 *   a tree height that allows the storage of at least numElements is chosen and
 *   then space is allocated for the case that this height is filled with nodes.
 * @tparam LessThanComparator - A binary functor that takes two keys and returns
 *   whether a is smaller than b. Must implement strict weak ordering. Defaults
 *   to std::less<KeyType>.
 * @tparam EqualityComparator - A binary functor that takes two keys and returns
 *   whether they are equal. Defaults to std::equal_to<KeyType>.
 */
template<
  typename KeyType,
  size_t minDegree,
  size_t numElements,
  class LessThanComparator = std::less<KeyType>,
  class EqualityComparator = std::equal_to<KeyType>
> class BTree {
private:
  //! Type for Nodes
  struct Node {
    static constexpr unsigned minKeys = minDegree - 1;
    static constexpr unsigned maxKeys = 2 * minDegree - 1;

    DynamicArray<KeyType, maxKeys> keys;
    DynamicArray<Node*, maxKeys + 1> children;

    constexpr Node() {}

    constexpr bool isLeaf() const {
      return children.size() == 0;
    }

    constexpr bool isFull() const {
      return keys.size() == maxKeys;
    }
  };

public:
  //! The height needed to be able to hold at least numElements keys
  static constexpr size_t maxHeight = BTreeProperties::minHeight(
    numElements,
    minDegree
  );

  //! The maximum number of nodes that the tree can hold
  static constexpr size_t maxNodes = BTreeProperties::maxNodesInTree(
    maxHeight,
    minDegree
  );

  //! The maximum number of keys that the tree can hold
  static constexpr size_t maxKeys = maxNodes * (2 * minDegree - 1);

private:
  //! Pointer to root node
  Node* _rootPtr; 

  //! Array holding all tree nodes
  DynamicArray<Node, maxNodes> _nodes;

  //! Array holding pointers to any 'deleted' tree nodes
  DynamicArray<Node*, maxNodes> _garbage;

  //! Less-than comparator instance
  LessThanComparator _lt;

  //! Equality comparator instance
  EqualityComparator _eq;

  //! 'Allocates' a new node and returns a pointer to it
  constexpr Node* _newNode() {
    if(_garbage.size() > 0) {
      Node* nodePtr = _garbage.back();
      _garbage.pop_back();

      // Refresh the node
      *nodePtr = Node {};

      return nodePtr;
    }

    if(_nodes.size() == maxNodes) {
      // We cannot get any new nodes! Nothing on the garbage
      throw "The maximum number of nodes has been reached for a BTree!";
    }

    _nodes.push_back(Node {});
    return &_nodes.back();
  }

  //! Marks a node as 'deleted' for recycling in _newNode
  constexpr void _markNodeDeleted(Node* nodePtr) {
    _garbage.push_back(nodePtr);
  }

  //! Recursive search for an element in a subtree rooted at node
  constexpr Node* _search(Node* node, const KeyType& key) const {
    auto keyLB = lowerBound<KeyType, LessThanComparator>(
      node->keys.begin(),
      node->keys.end(),
      key,
      _lt
    );

    // In case the lower bound is actually our sought key, return this node
    if(keyLB != node->keys.end() && _eq(*keyLB, key)) {
      return node;
    }

    // If we haven't found the key and this node is a leaf, search fails
    if(node->isLeaf()) {
      return nullptr;
    }

    /* Otherwise descend to the child at the same index as the lower bound in
     * keys (this also works in case LB is at keys.end() since children.size()
     * == keys.size() + 1) and look there
     */
    return _search(
      node->children.at(
        keyLB - node->keys.begin()
      ),
      key
    );
  }

  constexpr void _splitChild(Node* node, const unsigned& i) {
    // i is the child index in node's keys being split since that node is full

    // The node being split is afterwards considered the "left" node
    Node* left = node->children.at(i);

    // Allocate a new "right" node
    if(_nodes.size() == maxNodes) {
      throw "Inserting into full BTree";
    }

    auto right = _newNode();

    // Move keys
    right->keys = left->keys.splice(minDegree);

    // In case left is not a leaf, move the children too
    if(!left->isLeaf()) {
      right->children = left->children.splice(minDegree);
    }

    // Insert the original median key into the non-full node
    node->keys.insertAt(
      node->keys.begin() + i,
      left->keys.back()
    );

    // Have to remove it from left, too
    left->keys.pop_back();

    // And assign right as the child to the right of the inserted median key
    node->children.insertAt(
      node->children.begin() + i + 1,
      right
    );
  }

  constexpr void _insertNonFull(Node* node, const KeyType& key) {
    if(node->isLeaf()) {
      auto keyLB = lowerBound<KeyType, LessThanComparator>(
        node->keys.begin(),
        node->keys.end(),
        key,
        _lt
      );

      if(keyLB != node->keys.end() && _eq(*keyLB, key)) {
        throw "Inserting an already-existent key!";
      }

      node->keys.insertAt(keyLB, key);
    } else {
      // Where to go?
      auto keyLB = lowerBound<KeyType, LessThanComparator>(
        node->keys.begin(),
        node->keys.end(),
        key,
        _lt
      );

      auto childPos = keyLB - node->keys.begin();

      // In case the purported child is full, split it!
      if(node->children.at(childPos)->isFull()) {
        _splitChild(node, childPos);

        /* Keys has an additional key from the split, check if index has to be
         * incremented
         */
        if(_lt(node->keys.at(childPos), key)) {
          ++childPos;
        }
      }

      // The target child cannot be full anymore, so we can call
      _insertNonFull(node->children.at(childPos), key);
    }
  }

  constexpr bool _isRootNode(const Node* const node) const {
    return node == _rootPtr;
  }

  //! Returns the smallest leaf node in the sub-tree rooted at nodePtr
  constexpr Node* _smallestLeafNode(Node* nodePtr) const {
    while(!nodePtr->isLeaf()) {
      nodePtr = nodePtr->children.front();
    }

    return nodePtr;
  }

  //! Returns the largest leaf node in the sub-tree rooted at nodePtr
  constexpr Node* _largestLeafNode(Node* nodePtr) const {
    while(!nodePtr->isLeaf()) {
      nodePtr = nodePtr->children.back();
    }

    return nodePtr;
  }

  //! Recursively deletes a key from a sub-tree rooted at node
  constexpr void _delete(Node* node, const KeyType& key) {
    auto keyLB = lowerBound<KeyType, LessThanComparator>(
      node->keys.begin(),
      node->keys.end(),
      key,
      _lt
    );
    
    unsigned indexOfLB = keyLB - node->keys.begin();

    if(keyLB != node->keys.end() && _eq(*keyLB, key)) { 
      // Key to remove is in this node's keys
      if(node->isLeaf()) {
        // Case 1
        node->keys.removeAt(keyLB);
      } else {
        // Case 2
        if(node->children.at(indexOfLB)->keys.size() >= minDegree) {
          // Case 2a: Predecessor of key is maximum in subtree to the left
          
          // Predecessor key is largest key of largest leaf node in left subtree
          KeyType predecessor = _largestLeafNode(
            node->children.at(indexOfLB)
          )->keys.back();

          // Replace the key to be deleted by its predecessor
          *keyLB = predecessor;

          // Recursively delete the predecessor
          _delete(node->children.at(indexOfLB), predecessor);
        } else if(node->children.at(indexOfLB + 1)->keys.size() >= minDegree) {
          // Case 2b: Successor of key is minimum in subtree to the right

          // The successor key is the leftmost / smallest one
          KeyType successor = _smallestLeafNode(
            node->children.at(indexOfLB + 1)
          )->keys.front();

          // Replace the key to be deleted by its successor
          *keyLB = successor;

          // Recursively delete the successor
          _delete(node->children.at(indexOfLB + 1), successor);
        } else {
          /* Case 2c: Merge the key to delete, all of the right child into the
           * left child. The current node loses both k and the pointer to the
           * right child
           */

          Node* leftChild = node->children.at(indexOfLB);
          Node* rightChild = node->children.at(indexOfLB + 1);

          // Add the key to the left child
          leftChild->keys.push_back(key);

          // Merge the right child into the left child
          leftChild->keys.copyIn(rightChild->keys);
          if(!leftChild->isLeaf()) {
            leftChild->children.copyIn(rightChild->children);
          }

          // Remove the key and child pointer to rightChild from left
          node->keys.removeAt(keyLB);
          node->children.removeAt(
            node->children.begin() + indexOfLB + 1
          );

          // Delete the right child
          _markNodeDeleted(rightChild);

          // Delete the key recursively from the left child
          _delete(leftChild, key);
        }
      }
    } else { 
      /* Case 3: The key to delete is not in this node's keys and we have to
       * descend in the tree. Need to ensure that any node we descend to has
       * at least minDegree keys!
       */
      Node* targetChild = node->children.at(indexOfLB);

      if(targetChild->keys.size() == minDegree - 1) {
        // Case 3a Move some keys around from left or right siblings 

        if(
          indexOfLB != 0 
          && node->children.at(indexOfLB - 1)->keys.size() >= minDegree
        ) {
          Node* leftSibling = node->children.at(indexOfLB - 1);

          // Move key at LB into targetChild
          targetChild->keys.insertAt(
            targetChild->keys.begin(),
            node->keys.at(indexOfLB - 1)
          );

          // Last key of left sibling replaces key at LB
          node->keys.at(indexOfLB - 1) = leftSibling->keys.back();
          leftSibling->keys.pop_back();

          // In case it is not a leaf, we move the child pointer too
          if(!targetChild->isLeaf()) {
            targetChild->children.insertAt(
              targetChild->children.begin(),
              leftSibling->children.back()
            );

            leftSibling->children.pop_back();
          }
        } else if(
          indexOfLB < node->keys.size()
          && node->children.at(indexOfLB + 1)->keys.size() >= minDegree
        ) {
          Node* rightSibling = node->children.at(indexOfLB + 1);

          // Move key at LB into targetChild
          targetChild->keys.push_back(
            node->keys.at(indexOfLB)
          );

          // First key of right sibling replaces key at LB
          *keyLB = rightSibling->keys.front();
          rightSibling->keys.removeAt(
            rightSibling->keys.begin()
          );

          // In case the target is not a leaf, we move the child pointer too
          if(!targetChild->isLeaf()) {
            targetChild->children.push_back(
              rightSibling->children.front()
            );

            rightSibling->children.removeAt(
              rightSibling->children.begin()
            );
          }
        } else {
          // Case 3b

          if(indexOfLB != 0) { // Merge with left sibling
            Node* leftSibling = node->children.at(indexOfLB - 1);

            // Move key down to left sibling
            leftSibling->keys.push_back(
              node->keys.at(indexOfLB - 1)
            );
            node->keys.removeAt(
              node->keys.begin() + indexOfLB - 1
            );

            // Merge keys and children of targetChild into leftSibling
            leftSibling->keys.copyIn(targetChild->keys);
            if(!targetChild->isLeaf()) {
              leftSibling->children.copyIn(targetChild->children);
            }

            node->children.removeAt(
              node->children.begin() + indexOfLB
            );

            _markNodeDeleted(targetChild);

            --indexOfLB;
          } else { // Merge with right sibling
            Node* rightSibling = node->children.at(indexOfLB + 1);

            targetChild->keys.push_back(
              node->keys.at(indexOfLB)
            );
            node->keys.removeAt(
              node->keys.begin() + indexOfLB
            );

            targetChild->keys.copyIn(rightSibling->keys);
            if(!targetChild->isLeaf()) {
              targetChild->children.copyIn(rightSibling->children);
            }

            node->children.removeAt(
              node->children.begin() + indexOfLB + 1
            );

            _markNodeDeleted(rightSibling);
          }
        }
      }

      _delete(node->children.at(indexOfLB), key);
    }
  }

  //! Checks whether the node is a valid B-Tree node, and throws if anything is off
  constexpr void _validate(const Node& node) const {
    // The node should not be in the garbage
    auto foundIter = _garbage.begin();

    while(foundIter != _garbage.end()) {
      if(*foundIter == &node) {
        break;
      }

      ++foundIter;
    }

    // Skip if the current node is in the garbage
    if(foundIter != _garbage.end()) {
      throw "An active node is marked as garbage!";
    }

    // A non-root node has min. t-1 keys
    if(
      !_isRootNode(&node) 
      && node.keys.size() < minDegree - 1
    ) {
      throw "Not every internal node has min. t-1 keys!";
    }

    // Every internal node with n keys has n+1 children
    if(
      !_isRootNode(&node) 
      && !node.isLeaf() 
      && node.keys.size() != node.children.size() - 1
    ) {
      throw "Not every internal node with n keys has n+1 children!";
    }

    // Every key list is ordered
    if(!isTotallyOrdered(node.keys)) {
      throw "Not all key lists are totally ordered!";
    }

    /* Children 'between' keys have keys that are within the interval set by
     * the parents
     */
    if(!node.isLeaf()) {
      for(unsigned i = 1; i < node.children.size(); ++i) {
        if(
          !_lt(
            node.children.at(i - 1)->keys.back(),
            node.keys.at(i - 1)
          ) || !_lt(
            node.keys.at(i - 1),
            node.children.at(i)->keys.front()
          )
        ) {
          throw "Not all children's keys are bounded by the parent!";
        }
      }
    }
  }

public:
  constexpr BTree() {
    _rootPtr = _newNode();
  }

  /*!
   * Inserts a new key into the tree. This key may *not* already be in the tree.
   * Complexity is O(t log_t N), where t is the minimum degree of the tree and
   * N the number of contained elements.
   */
  constexpr void insert(const KeyType& key) {
    Node* r = _rootPtr;

    if(r->isFull()) { // Root is full, must be split
      if(_nodes.size() == maxNodes) {
        throw "Inserting into full BTree!";
      }

      auto s = _newNode();

      _rootPtr = s;

      s->children.push_back(r);
      _splitChild(s, 0);
      _insertNonFull(s, key);
    } else {
      _insertNonFull(r, key);
    }
  }

  /*!
   * Check whether a key is stored in the tree. The complexity of this operation
   * is O(t log_t N), where t is the minimum degree of the tree and N the
   * number of contained elements.
   */
  constexpr bool contains(const KeyType& key) const {
    Node* foundPtr = _search(_rootPtr, key);
    return foundPtr != nullptr;
  }

  /*!
   * Deletes a key from the tree. This key must exist in the tree. The
   * complexity of this operation is O(t log_t N), where t is the minimum
   * degree of the tree and N the number of contained elements.
   */
  constexpr void remove(const KeyType& key) {
    _delete(_rootPtr, key);

    // In case the root node is keyless but has a child, shrink the tree
    if(_rootPtr->keys.size() == 0 && !_rootPtr->isLeaf()) {
      Node* emptyRoot = _rootPtr;

      _rootPtr = _rootPtr->children.front();

      _markNodeDeleted(emptyRoot);
    }
  }

  //! Dumps a graphViz representation of the B-Tree.
  std::string dumpGraphviz() const {
    using namespace std::string_literals;

    std::stringstream graph;
    graph << "digraph g {\n"
      << "  node [shape=record, height=.1]\n\n";

    std::set<const Node*> garbagePtrs {_garbage.begin(), _garbage.end()};
    std::map<const Node*, std::string> nodeNames;
    unsigned nodeCounter = 0;

    auto getName = [&](const Node* const nodePtr) -> std::string {
      if(nodeNames.count(nodePtr) == 0) {
        nodeNames[nodePtr] = "node"s + std::to_string(nodeCounter);
        ++nodeCounter;
      }
      
      return nodeNames.at(nodePtr);
    };

    for(auto it = _nodes.begin(); it != _nodes.end(); ++it) {
      const auto& node = *it;

      // Ensure the current node is not in the garbage
      if(garbagePtrs.count(&node) == 0) {
        graph << "  " << getName(&node) << "[label=\"";

        if(node.isLeaf()) {
          for(unsigned i = 0; i < node.keys.size(); ++i) {
            graph << node.keys.at(i);
            if(i != node.keys.size() - 1) {
              graph << "|";
            }
          }
        } else {
          graph << "<f0>|";
          for(unsigned i = 1; i < node.children.size(); ++i) {
            graph << node.keys.at(i-1) << "|<f" << i << ">";
            if(i != node.children.size() - 1) {
              graph << "|";
            }
          }
        }

        graph << "\"];\n";

        // Write all connections
        if(!node.isLeaf()) {
          for(unsigned i = 0; i < node.children.size(); ++i) {
            graph << "  \"" << getName(&node) << "\":f" << i 
              << " -> \"" << getName(node.children.at(i)) << "\";\n";
          }
        }
      }
    }

    graph << "}";

    return graph.str();
  }

  //! Validates the state of the tree by DFS traversal. Throws if anything is off.
  constexpr void validate() const {
    DynamicArray<const Node*, maxNodes> stack {_rootPtr};

    while(stack.size() > 0) {
      const Node* nodePtr = stack.back();
      stack.pop_back();

      for(auto& childPtr : nodePtr->children) {
        stack.push_back(childPtr);
      }

      _validate(*nodePtr);
    }
  }

  /* TODO
   * - key iterators (in order traversal)
   * - node iterators?
   */
};

} // namespace ConstexprMagic

#endif
