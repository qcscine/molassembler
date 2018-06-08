#ifndef INCLUDE_CONSTEXPR_MAGIC_BTREE_H
#define INCLUDE_CONSTEXPR_MAGIC_BTREE_H

#include "DynamicArray.h"
#include "Math.h"
#include "Optional.h"

#include <sstream>

/*! @file
 *
 * Implements a BTree which doesn't store key-value pairs, only keys.
 */


/* TODO
 * - node constIterators?
 */

namespace temple {

namespace BTreeProperties {

/*!
 * Calculates the maximum number of nodes in a B-tree of a specific minDegree and
 * height.
 */
constexpr size_t maxNodesInTree(const size_t height, const size_t minDegree) PURITY_STRONG;
constexpr size_t maxNodesInTree(const size_t height, const size_t minDegree) {
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
constexpr size_t minHeight(const size_t numKeys, const size_t minDegree) PURITY_STRONG;
constexpr size_t minHeight(const size_t numKeys, const size_t minDegree) {
  return Math::ceil(
    Math::log(
      static_cast<double>(numKeys + 1),
      static_cast<double>(2 * minDegree)
    ) - 1
  );
}

/*!
 * Calculates the maximal height bound needed for a B-Tree of a specific
 * minDegree to be able to hold a certain number of keys. The actual maximal
 * height may be lower.
 */
constexpr size_t maxHeightBound(const size_t numKeys, const size_t minDegree) PURITY_STRONG;
constexpr size_t maxHeightBound(const size_t numKeys, const size_t minDegree) {
  return Math::floor(
    Math::log(
      static_cast<double>(numKeys + 1) / 2,
      static_cast<double>(minDegree)
    )
  );
}

} // namespace BTreeImpl

/*! A B-Tree that merely stores keys, not key-value pairs.
 *
 * This class is a B-Tree that stores keys in its nodes, not key-value pairs.
 * (This can be turned into an associative container with a simple modification,
 * though: Set KeyType to std::pair<KeyType, ValueType> and supply custom
 * LessThanComparator and EqualityComparators that merely compare using the
 * KeyType. Then return the ValueType from a lookup with a default-constructed
 * ValueType.)
 *
 * @tparam KeyType The desired stored type
 * @tparam minDegree The minimum degree t of the B-Tree. Must be >= 2, since
 *   Nodes store a minimum of t-1 keys. A Node of degree t can store up to 2t-1
 *   keys and has up to 2t children.
 * @tparam numElements The intended maximimum amount of stored elements. The
 *   class must allocate all nodes that may possibly be stored at any time, so
 *   a tree height that allows the storage of at least numElements is chosen and
 *   then space is allocated for the case that this height is filled with nodes.
 *   This has the consequence that the instantiated tree can typically store
 *   quite a bit more elements than originally intended (see the static member
 *   maxKeys). Play with the minimum order a bit if you need space-efficiency.
 * @tparam LessThanComparator A pure binary functor that takes two keys and
 *   returns whether a is smaller than b. Must implement strict weak ordering.
 *   Defaults to std::less<KeyType>.
 * @tparam EqualityComparator A pure binary functor that takes two keys and
 *   returns whether they are equal. Defaults to std::equal_to<KeyType>.
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
    DynamicArray<unsigned, maxKeys + 1> children;

    constexpr Node() {}

    constexpr bool isLeaf() const PURITY_WEAK {
      return children.size() == 0;
    }

    constexpr bool isFull() const PURITY_WEAK {
      return keys.size() == maxKeys;
    }
  };

public:
  //! The height needed to be able to hold at least numElements keys
  static constexpr size_t maxHeight = BTreeProperties::maxHeightBound(
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
  unsigned _rootPtr;

  //! Array holding all tree nodes
  DynamicArray<Node, maxNodes> _nodes;

  //! Array holding pointers to any 'deleted' tree nodes
  DynamicArray<unsigned, maxNodes> _garbage;

  //! Less-than comparator instance
  LessThanComparator _lt;

  //! Equality comparator instance
  EqualityComparator _eq;

  //! Number of contained elements
  unsigned _count = 0;

  //! 'Allocates' a new node and returns a pointer to it
  constexpr unsigned _newNode() {
    if(_garbage.size() > 0) {
      unsigned newNodeIndex = _garbage.back();
      _garbage.pop_back();

      // Refresh the node
      _nodes.at(newNodeIndex) = Node {};

      return newNodeIndex;
    }

    if(_nodes.size() == maxNodes) {
      // We cannot get any new nodes! Nothing on the garbage
      throw "The maximum number of nodes has been reached for a BTree!";
    }

    _nodes.push_back(Node {});
    return _nodes.size() - 1;
  }

  //! Marks a node as 'deleted' for recycling in _newNode
  constexpr void _markNodeDeleted(unsigned nodeIndex) {
    _garbage.push_back(nodeIndex);
  }

  constexpr Node& _getNode(unsigned nodeIndex) {
    return _nodes.at(nodeIndex);
  }

  constexpr const Node& _getNode(unsigned nodeIndex) const PURITY_WEAK {
    return _nodes.at(nodeIndex);
  }

  //! Recursive search for an element in a subtree rooted at node
  constexpr Optional<unsigned> _search(unsigned nodeIndex, const KeyType& key) const PURITY_WEAK {
    auto node = _getNode(nodeIndex);

    auto keyLB = lowerBound<KeyType, LessThanComparator>(
      node.keys.begin(),
      node.keys.end(),
      key,
      _lt
    );

    // In case the lower bound is actually our sought key, return this node
    if(keyLB != node.keys.end() && _eq(*keyLB, key)) {
      return nodeIndex;
    }

    // If we haven't found the key and this node is a leaf, search fails
    if(node.isLeaf()) {
      return {};
    }

    /* Otherwise descend to the child at the same index as the lower bound in
     * keys (this also works in case LB is at keys.end() since children.size()
     * == keys.size() + 1) and look there
     */
    return _search(
      node.children.at(
        keyLB - node.keys.begin()
      ),
      key
    );
  }

  constexpr void _splitChild(unsigned nodeIndex, const unsigned& i) {
    // i is the child index in node's keys being split since that node is full
    auto& parent = _getNode(nodeIndex);

    // The node being split is afterwards considered the "left" node
    auto& left = _getNode(
      _getNode(nodeIndex).children.at(i)
    );

    // Allocate a new "right" node
    if(_nodes.size() == maxNodes) {
      throw "Inserting into full BTree";
    }

    auto rightIndex = _newNode();
    auto& right = _getNode(rightIndex);

    // Move keys
    right.keys = left.keys.splice(minDegree);

    // In case left is not a leaf, move the children too
    if(!left.isLeaf()) {
      right.children = left.children.splice(minDegree);
    }

    // Insert the original median key into the non-full node
    parent.keys.insertAt(
      parent.keys.begin() + i,
      left.keys.back()
    );

    // Have to remove it from left, too
    left.keys.pop_back();

    // And assign right as the child to the right of the inserted median key
    parent.children.insertAt(
      parent.children.begin() + i + 1,
      rightIndex
    );
  }

  constexpr void _insertNonFull(unsigned nodeIndex, const KeyType& key) {
    auto& node = _getNode(nodeIndex);

    if(node.isLeaf()) {
      auto keyLB = lowerBound<KeyType, LessThanComparator>(
        node.keys.begin(),
        node.keys.end(),
        key,
        _lt
      );

      if(keyLB != node.keys.end() && _eq(*keyLB, key)) {
        throw "Inserting an already-existent key!";
      }

      node.keys.insertAt(keyLB, key);
    } else {
      // Where to go?
      auto keyLB = lowerBound<KeyType, LessThanComparator>(
        node.keys.begin(),
        node.keys.end(),
        key,
        _lt
      );

      auto childPos = keyLB - node.keys.begin();

      // In case the purported child is full, split it!
      if(_getNode(node.children.at(childPos)).isFull()) {
        _splitChild(nodeIndex, childPos);

        /* Keys has an additional key from the split, check if index has to be
         * incremented
         */
        if(_lt(node.keys.at(childPos), key)) {
          ++childPos;
        }
      }

      // The target child cannot be full anymore, so we can call
      _insertNonFull(node.children.at(childPos), key);
    }
  }

  constexpr bool _isRootNode(unsigned nodeIndex) const PURITY_WEAK {
    return nodeIndex == _rootPtr;
  }

  //! Returns the smallest leaf node in the sub-tree rooted at nodePtr
  constexpr unsigned _smallestLeafNode(unsigned nodeIndex) const PURITY_WEAK {
    while(!_getNode(nodeIndex).isLeaf()) {
      nodeIndex = _getNode(nodeIndex).children.front();
    }

    return nodeIndex;
  }

  //! Returns the largest leaf node in the sub-tree rooted at nodePtr
  constexpr unsigned _largestLeafNode(unsigned nodeIndex) const PURITY_WEAK {
    while(!_getNode(nodeIndex).isLeaf()) {
      nodeIndex = _getNode(nodeIndex).children.back();
    }

    return nodeIndex;
  }

  //! Recursively deletes a key from a sub-tree rooted at node
  constexpr void _delete(unsigned nodeIndex, const KeyType& key) {
    auto& node = _getNode(nodeIndex);

    auto keyLB = lowerBound<KeyType, LessThanComparator>(
      node.keys.begin(),
      node.keys.end(),
      key,
      _lt
    );

    unsigned indexOfLB = keyLB - node.keys.begin();

    if(keyLB != node.keys.end() && _eq(*keyLB, key)) {
      // Key to remove is in this node's keys
      if(node.isLeaf()) {
        // Case 1
        node.keys.removeAt(keyLB);
      } else {
        // Case 2
        if(_getNode(node.children.at(indexOfLB)).keys.size() >= minDegree) {
          // Case 2a: Predecessor of key is maximum in subtree to the left

          // Predecessor key is largest key of largest leaf node in left subtree
          KeyType predecessor = _getNode(_largestLeafNode(
            node.children.at(indexOfLB)
          )).keys.back();

          // Replace the key to be deleted by its predecessor
          *keyLB = predecessor;

          // Recursively delete the predecessor
          _delete(node.children.at(indexOfLB), predecessor);
        } else if(_getNode(node.children.at(indexOfLB + 1)).keys.size() >= minDegree) {
          // Case 2b: Successor of key is minimum in subtree to the right

          // The successor key is the leftmost / smallest one
          KeyType successor = _getNode(_smallestLeafNode(
            node.children.at(indexOfLB + 1)
          )).keys.front();

          // Replace the key to be deleted by its successor
          *keyLB = successor;

          // Recursively delete the successor
          _delete(node.children.at(indexOfLB + 1), successor);
        } else {
          /* Case 2c: Merge the key to delete, all of the right child into the
           * left child. The current node loses both k and the pointer to the
           * right child
           */

          unsigned leftChildIndex = node.children.at(indexOfLB);
          auto& leftChild = _getNode(leftChildIndex);

          unsigned rightChildIndex = node.children.at(indexOfLB + 1);
          auto& rightChild = _getNode(rightChildIndex);

          // Add the key to the left child
          leftChild.keys.push_back(key);

          // Merge the right child into the left child
          leftChild.keys.copyIn(rightChild.keys);
          if(!leftChild.isLeaf()) {
            leftChild.children.copyIn(rightChild.children);
          }

          // Remove the key and child pointer to rightChild from left
          node.keys.removeAt(keyLB);
          node.children.removeAt(
            node.children.begin() + indexOfLB + 1
          );

          // Delete the right child
          _markNodeDeleted(rightChildIndex);

          // Delete the key recursively from the left child
          _delete(leftChildIndex, key);
        }
      }
    } else {
      /* Case 3: The key to delete is not in this node's keys and we have to
       * descend in the tree. Need to ensure that any node we descend to has
       * at least minDegree keys!
       */
      unsigned targetChildIndex = node.children.at(indexOfLB);
      auto& targetChild = _getNode(targetChildIndex);

      if(targetChild.keys.size() == minDegree - 1) {
        // Case 3a Move some keys around from left or right siblings

        if(
          indexOfLB != 0
          && _getNode(node.children.at(indexOfLB - 1)).keys.size() >= minDegree
        ) {
          unsigned leftSiblingIndex = node.children.at(indexOfLB - 1);
          auto& leftSibling = _getNode(leftSiblingIndex);

          // Move key at LB into targetChild
          targetChild.keys.insertAt(
            targetChild.keys.begin(),
            node.keys.at(indexOfLB - 1)
          );

          // Last key of left sibling replaces key at LB
          node.keys.at(indexOfLB - 1) = leftSibling.keys.back();
          leftSibling.keys.pop_back();

          // In case it is not a leaf, we move the child pointer too
          if(!targetChild.isLeaf()) {
            targetChild.children.insertAt(
              targetChild.children.begin(),
              leftSibling.children.back()
            );

            leftSibling.children.pop_back();
          }
        } else if(
          indexOfLB < node.keys.size()
          && _getNode(node.children.at(indexOfLB + 1)).keys.size() >= minDegree
        ) {
          unsigned rightSiblingIndex = node.children.at(indexOfLB + 1);
          auto& rightSibling = _getNode(rightSiblingIndex);

          // Move key at LB into targetChild
          targetChild.keys.push_back(
            node.keys.at(indexOfLB)
          );

          // First key of right sibling replaces key at LB
          *keyLB = rightSibling.keys.front();
          rightSibling.keys.removeAt(
            rightSibling.keys.begin()
          );

          // In case the target is not a leaf, we move the child pointer too
          if(!targetChild.isLeaf()) {
            targetChild.children.push_back(
              rightSibling.children.front()
            );

            rightSibling.children.removeAt(
              rightSibling.children.begin()
            );
          }
        } else {
          // Case 3b

          if(indexOfLB != 0) { // Merge with left sibling
            unsigned leftSiblingIndex = node.children.at(indexOfLB - 1);
            auto& leftSibling = _getNode(leftSiblingIndex);

            // Move key down to left sibling
            leftSibling.keys.push_back(
              node.keys.at(indexOfLB - 1)
            );
            node.keys.removeAt(
              node.keys.begin() + indexOfLB - 1
            );

            // Merge keys and children of targetChild into leftSibling
            leftSibling.keys.copyIn(targetChild.keys);
            if(!targetChild.isLeaf()) {
              leftSibling.children.copyIn(targetChild.children);
            }

            node.children.removeAt(
              node.children.begin() + indexOfLB
            );

            _markNodeDeleted(targetChildIndex);

            --indexOfLB;
          } else { // Merge with right sibling
            unsigned rightSiblingIndex = node.children.at(indexOfLB + 1);
            auto& rightSibling = _getNode(rightSiblingIndex);

            targetChild.keys.push_back(
              node.keys.at(indexOfLB)
            );
            node.keys.removeAt(
              node.keys.begin() + indexOfLB
            );

            targetChild.keys.copyIn(rightSibling.keys);
            if(!targetChild.isLeaf()) {
              targetChild.children.copyIn(rightSibling.children);
            }

            node.children.removeAt(
              node.children.begin() + indexOfLB + 1
            );

            _markNodeDeleted(rightSiblingIndex);
          }
        }
      }

      _delete(node.children.at(indexOfLB), key);
    }
  }

  //! Checks whether the node is a valid B-Tree node, and throws if anything is off
  constexpr void _validate(unsigned nodeIndex) const {
    // The node should not be in the garbage
    auto foundIter = _garbage.begin();

    while(foundIter != _garbage.end()) {
      if(*foundIter == nodeIndex) {
        break;
      }

      ++foundIter;
    }

    // Skip if the current node is in the garbage
    if(foundIter != _garbage.end()) {
      throw "An active node is marked as garbage!";
    }

    auto& node = _getNode(nodeIndex);

    // A non-root node has min. t-1 keys
    if(
      !_isRootNode(nodeIndex)
      && node.keys.size() < minDegree - 1
    ) {
      throw "Not every internal node has min. t-1 keys!";
    }

    // Every internal node with n keys has n+1 children
    if(
      !_isRootNode(nodeIndex)
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
            _getNode(node.children.at(i - 1)).keys.back(),
            node.keys.at(i - 1)
          ) || !_lt(
            node.keys.at(i - 1),
            _getNode(node.children.at(i)).keys.front()
          )
        ) {
          throw "Not all children's keys are bounded by the parent!";
        }
      }
    }
  }

public:
  constexpr BTree() : _rootPtr {0} {
    _rootPtr = _newNode();
  }

  /*! Add a key to the tree.
   *
   * Inserts a new key into the tree. This key must not already be in the tree.
   * Complexity is O(t log_t N), where t is the minimum degree of the tree and
   * N the number of contained elements.
   */
  constexpr void insert(const KeyType& key) {
    unsigned r = _rootPtr;

    if(_getNode(r).isFull()) { // Root is full, must be split
      unsigned s = _newNode();

      _rootPtr = s;

      _getNode(s).children.push_back(r);
      _splitChild(s, 0);
      _insertNonFull(s, key);
    } else {
      _insertNonFull(r, key);
    }

    ++_count;
  }

  /*! Check whether a key is stored in the tree.
   *
   * Check whether a key is stored in the tree. The complexity of this operation
   * is O(t log_t N), where t is the minimum degree of the tree and N the
   * number of contained elements.
   */
  constexpr bool contains(const KeyType& key) const PURITY_WEAK {
    auto nodeIndexOptional = _search(_rootPtr, key);

    return nodeIndexOptional.hasValue();
  }

  constexpr Optional<KeyType> getOption(const KeyType& key) const PURITY_WEAK {
    auto nodeIndexOptional = _search(_rootPtr, key);

    if(!nodeIndexOptional.hasValue()) {
      return {};
    }

    auto keyLB = lowerBound<KeyType, LessThanComparator>(
      _getNode(nodeIndexOptional.value()).keys.begin(),
      _getNode(nodeIndexOptional.value()).keys.end(),
      key,
      _lt
    );

    return *keyLB;
  }

  /*!
   * Deletes a key from the tree. This key must exist in the tree. The
   * complexity of this operation is O(t log_t N), where t is the minimum
   * degree of the tree and N the number of contained elements.
   */
  constexpr void remove(const KeyType& key) {
    _delete(_rootPtr, key);

    // In case the root node is keyless but has a child, shrink the tree
    if(_getNode(_rootPtr).keys.size() == 0 && !_getNode(_rootPtr).isLeaf()) {
      unsigned emptyRoot = _rootPtr;

      _rootPtr = _getNode(_rootPtr).children.front();

      _markNodeDeleted(emptyRoot);
    }

    --_count;
  }

  //! Dumps a graphViz representation of the B-Tree.
  std::string dumpGraphviz() const PURITY_WEAK {
    using namespace std::string_literals;

    std::stringstream graph;
    graph << "digraph g {\n"
      << "  node [shape=record, height=.1]\n\n";

    DynamicArray<unsigned, maxNodes> stack {_rootPtr};

    while(stack.size() > 0) {
      unsigned nodeIndex = stack.back();
      stack.pop_back();

      for(auto& childIndex : _getNode(nodeIndex).children) {
        stack.push_back(childIndex);
      }

      const auto& node = _getNode(nodeIndex);

      graph << "  node" << nodeIndex << "[label=\"";

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
          graph << "  \"node" << nodeIndex << "\":f" << i
            << " -> \"node" << node.children.at(i) << "\";\n";
        }
      }
    }

    graph << "}";

    return graph.str();
  }

  //! Validates the state of the tree by DFS traversal. Throws if anything is off.
  constexpr void validate() const {
    DynamicArray<unsigned, maxNodes> stack {_rootPtr};

    while(stack.size() > 0) {
      unsigned nodeIndex = stack.back();
      auto& node = _getNode(nodeIndex);
      stack.pop_back();

      for(auto& childIndex : node.children) {
        stack.push_back(childIndex);
      }

      _validate(nodeIndex);
    }
  }

  //! Returns the number of elements in the tree
  constexpr unsigned size() const PURITY_WEAK {
    return _count;
  }

  using ValueIteratorBase = std::iterator<
    std::bidirectional_iterator_tag, // iterator category
    KeyType,                         // value_type
    int,                             // difference_type
    const KeyType* const,            // pointer
    const KeyType&                   // reference
  >;

  class constIterator : public ValueIteratorBase {
  public:
    enum class InitializeAs : unsigned {
      Begin = 0,
      End = 1
    };

  private:
    const BTree& _baseRef;
    const unsigned _leftMostNode;
    const unsigned _rightMostNode;
    DynamicArray<unsigned, maxHeight + 1> _nodeStack;
    DynamicArray<unsigned, maxHeight + 1> _indexStack;

    constexpr const Node& _getCurrentNode() const PURITY_WEAK {
      return _baseRef._getNode(_nodeStack.back());
    }

    constexpr unsigned _currentNodeIndexLimit() const PURITY_WEAK {
      // For leaves, the past-the-end position is the size of keys
      if(_getCurrentNode().isLeaf()) {
        return _getCurrentNode().keys.size();
      }

      /* For internal nodes, the past-the-end position is the size of keys
       * plus the size of children + 1
       */
      return 2 * _getCurrentNode().keys.size() + 1;
    }

  public:
    constexpr constIterator(
      const BTree& tree,
      const InitializeAs& initDecision
    ) : _baseRef(tree),
        _leftMostNode(tree._smallestLeafNode(tree._rootPtr)),
        _rightMostNode(tree._largestLeafNode(tree._rootPtr)),
        _nodeStack {tree._rootPtr}
    {
      if(!static_cast<bool>(static_cast<unsigned>(initDecision))) {
        // 0 is Begin, so not-false is begin

        while(!_getCurrentNode().isLeaf()) {
          _indexStack.push_back(0);

          _nodeStack.push_back(
            _getCurrentNode().children.front()
          );
        }

        // At pos 0 of the leaf indices
        _indexStack.push_back(0);
      } else {
        // End constIterator initialization

        while(!_getCurrentNode().isLeaf()) {
          _indexStack.push_back(
            2 * _getCurrentNode().keys.size()
          );

          _nodeStack.push_back(
            _getCurrentNode().children.back()
          );
        }

        // past-the-end of indices
        _indexStack.push_back(_getCurrentNode().keys.size());
      }
    }

    constexpr constIterator(const constIterator& other)
      : _baseRef(other._baseRef),
        _leftMostNode(other._leftMostNode),
        _rightMostNode(other._rightMostNode),
        _nodeStack(other._nodeStack),
        _indexStack(other._indexStack)
    {}

    constexpr constIterator& operator = (const constIterator& other) {
      if(this->_baseRef != other._baseRef) {
        throw "Assigning BTree constIterator to another base instance!";
      }

      _leftMostNode = other._leftMostNode;
      _rightMostNode = other._rightMostNode;
      _nodeStack = other._nodeStack;
      _indexStack = other._indexStack;
    }

    constexpr constIterator& operator ++ () {
      auto indexLimit = _currentNodeIndexLimit();

      if(_indexStack.back() == indexLimit) { // We are already the end constIterator
        // Do nothing and return immediately
        return *this;
      }

      // In case we are a leaf, increment and re-check
      if(_getCurrentNode().isLeaf()) {
        ++_indexStack.back();

        /* If we hit the index limit for the node and we're not the rightmost
         * node, we have to go up the tree and to the right
         */
        if(
          _indexStack.back() == indexLimit
          && _nodeStack.back() != _rightMostNode
        ) {
          // Unwind the stack until we are at an incrementable position
          do {
            _indexStack.pop_back();
            _nodeStack.pop_back();
          } while(_indexStack.back() >= _currentNodeIndexLimit() - 1);

          // Increment here, now we are placed on a key at an internal node
          ++_indexStack.back();
        }

        return *this;
      }

      // We are on an internal node, incrementing puts us on a child
      ++_indexStack.back();
      _nodeStack.push_back(
        _getCurrentNode().children.at(
          _indexStack.back() / 2 // children are at even indices
        )
      );

      while(!_getCurrentNode().isLeaf()) {
        _indexStack.push_back(0);
        _nodeStack.push_back(
          _getCurrentNode().children.front()
        );
      }

      // Now we are a leaf, and placed on the first key
      _indexStack.push_back(0);

      return *this;
    }

    constexpr constIterator operator ++ (int) {
      constIterator retval = *this;
      ++(*this);
      return retval;
    }

    constexpr constIterator& operator -- () {
      if(_nodeStack.back() == _leftMostNode && _indexStack.back() == 0) {
        // We are the begin constIterator
        return *this;
      }

      // In case we are a leaf, decrement and re-check
      if(_getCurrentNode().isLeaf()) {
        // If we are at zero, we have to find a decrementable position
        if(_indexStack.back() == 0) {
          // Unwind the stack until we can decrement
          do {
            _indexStack.pop_back();
            _nodeStack.pop_back();
          } while(_indexStack.back() == 0);
        }

        // Decrement and return
        --_indexStack.back();
        return *this;
      }

      // We are an internal node, decrementing puts us on a child
      --_indexStack.back();
      _nodeStack.push_back(
        _getCurrentNode().children.at(
          _indexStack.back() / 2
        )
      );

      while(!_getCurrentNode().isLeaf()) {
        _indexStack.push_back(
          2 * _getCurrentNode().keys.size()
        );
        _nodeStack.push_back(
          _getCurrentNode().children.back()
        );
      }

      _indexStack.push_back(
         _getCurrentNode().keys.size() - 1
      );

      return *this;
    }

    constexpr constIterator operator -- (int) {
      constIterator retval = *this;
      --(*this);
      return retval;
    }

    constexpr bool operator == (const constIterator& other) const PURITY_WEAK {
      return (
        _nodeStack == other._nodeStack
        && _indexStack == other._indexStack
        && _leftMostNode == other._leftMostNode
        && _rightMostNode == other._rightMostNode
      );
    }

    constexpr bool operator != (const constIterator& other) const PURITY_WEAK {
      return !(
        *this == other
      );
    }

    constexpr typename ValueIteratorBase::reference operator *() const PURITY_WEAK {
      if(_getCurrentNode().isLeaf()) {
        return _getCurrentNode().keys.at(
          _indexStack.back()
        );
      }

      return _getCurrentNode().keys.at(
        _indexStack.back() / 2
      );
    }
  };

  //! Clears the tree. This operation is O(N)
  constexpr void clear() {
    // Refresh all nodes
    for(auto& node: _nodes) {
      node = Node {};
    }

    // Clear the garbage
    _garbage.clear();

    // Assign a new root
    _rootPtr = _newNode();
  }

  //! Alias for STL algorithm compatibility
  using const_iterator = constIterator;

  //! Returns a const iterator to the first key in the tree
  constexpr constIterator begin() const PURITY_WEAK {
    return constIterator(
      *this,
      constIterator::InitializeAs::Begin
    );
  }

  //! Returns a past-the-end const iterator
  constexpr constIterator end() const PURITY_WEAK {
    return constIterator(
      *this,
      constIterator::InitializeAs::End
    );
  }

  constexpr bool operator < (const BTree& other) const PURITY_WEAK {
    if(this->size() < other.size()) {
      return true;
    }

    if(this->size() > other.size()) {
      return false;
    }

    auto thisIter = this->begin();
    auto thisEnd = this->end();

    auto otherIter = other.begin();

    while(thisIter != thisEnd) {
      if(*thisIter < *otherIter) {
        return true;
      }

      if(*thisIter > *otherIter) {
        return false;
      }

      ++thisIter;
      ++otherIter;
    }

    return false;
  }

  constexpr bool operator > (const BTree& other) const PURITY_WEAK {
    return (other < *this);
  }

  constexpr bool operator == (const BTree& other) const PURITY_WEAK {
    if(this->size() != other.size()) {
      return false;
    }

    auto thisIter = this->begin();
    auto thisEnd = this->end();

    auto otherIter = other.begin();

    while(thisIter != thisEnd) {
      if(*thisIter != *otherIter) {
        return false;
      }

      ++thisIter;
      ++otherIter;
    }

    return true;
  }

  constexpr bool operator != (const BTree& other) const PURITY_WEAK {
    return !(*this == other);
  }
};

} // namespace temple

#endif
