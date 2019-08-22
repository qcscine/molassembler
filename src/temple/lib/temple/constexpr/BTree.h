/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief constexpr BTree
 *
 * Implements a constexpr BTree. Can be used for an ordered set-like container
 * with good complexity guarantees and turned into an associative container
 * with some comparator tweaking.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_BTREE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_BTREE_H

#include "temple/constexpr/Containers.h"
#include "temple/constexpr/DynamicArray.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Optional.h"

#include <sstream>

namespace temple {

namespace BTreeProperties {

/*!
 * Calculates the maximum number of nodes in a B-tree of a specific minDegree and
 * height.
 *
 * @complexity{@math{\Theta(1)}}
 */
PURITY_STRONG constexpr size_t maxNodesInTree(const size_t height, const size_t minDegree) {
  return static_cast<double>(
    Math::pow(2 * minDegree, static_cast<unsigned>(height + 1)) - 1
  ) / (
    2 * minDegree - 1
  );
}

/*!
 * Calculates the minimal height needed for a B-Tree of a specific minDegree to
 * be able to hold a certain number of values
 *
 * @complexity{@math{\Theta(1)}}
 */
PURITY_STRONG constexpr size_t minHeight(const size_t numValues, const size_t minDegree) {
  return Math::ceil(
    Math::log(
      static_cast<double>(numValues + 1),
      static_cast<double>(2 * minDegree)
    ) - 1
  );
}

/*!
 * Calculates the maximal height bound needed for a B-Tree of a specific
 * minDegree to be able to hold a certain number of values. The actual maximal
 * height may be lower.
 *
 * @complexity{@math{\Theta(1)}}
 */
PURITY_STRONG constexpr size_t maxHeightBound(const size_t numValues, const size_t minDegree) {
  return Math::floor(
    Math::log(
      static_cast<double>(numValues + 1) / 2,
      static_cast<double>(minDegree)
    )
  );
}

} // namespace BTreeImpl

/*! @brief A constexpr B-Tree
 *
 * This class is a B-Tree that stores values in its nodes.
 * (This can be turned into an associative container with a simple modification,
 * though: Set T to std::pair<KeyType, ValueType> and supply custom
 * LessThanComparator and EqualityComparators that merely compare using the
 * KeyType. Then return the ValueType from a lookup with a default-constructed
 * ValueType.)
 *
 * @tparam T The desired stored type
 * @tparam minDegree The minimum degree t of the B-Tree. Must be >= 2, since
 *   Nodes store a minimum of t-1 values. A Node of degree t can store up to 2t-1
 *   values and has up to 2t children.
 * @tparam numElements The intended maximimum amount of stored elements. The
 *   class must allocate all nodes that may possibly be stored at any time, so
 *   a tree height that allows the storage of at least numElements is chosen and
 *   then space is allocated for the case that this height is filled with nodes.
 *   This has the consequence that the instantiated tree can typically store
 *   quite a bit more elements than originally intended (see the static member
 *   maxElements). Play with the minimum order a bit if you need space-efficiency.
 * @tparam LessThanComparator A pure binary functor that takes two values and
 *   returns whether a is smaller than b. Must implement strict weak ordering.
 *   Defaults to std::less<T>.
 * @tparam EqualityComparator A pure binary functor that takes two values and
 *   returns whether they are equal. Defaults to std::equal_to<T>.
 *
 * @warning Since we are working in constexpr, we have to work around a large
 *   range of constraints, one of which is that we cannot dynamically allocate
 *   or delete objects. This is particularly bad for BTrees since the height
 *   of the tree can vary by a full level for the same contained information,
 *   so the total amount of allocated values may exceed the desired maximum
 *   number of elements by a factor proportional to the minimum degree. The
 *   amount of pre-allocated values is accessible via the maxElements static
 *   member.
 *
 * @parblock @note Nodes of the tree cannot be truly allocated or deleted
 * dynamically, so they are all pre-allocated on construction and merely marked
 * deleted.  Indices into the pre-allocated array of nodes are the 'pointers'
 * of this implementation.
 * @endparblock
 *
 * @parblock @note No differentiation is made in the type system for internal
 * or leaf nodes. Leaves just have no children.
 * @endparblock
 *
 * @parblock @note Regarding complexity annotations, @math{t} is used as the
 * symbol for the minimum degree of the tree, and @math{N} is the number of
 * contained elements.
 * @endparblock
 */
template<
  typename T,
  size_t minDegree,
  size_t numElements,
  class LessThanComparator = std::less<T>,
  class EqualityComparator = std::equal_to<T>
> class BTree {
private:
  // Forward-declare Node
  struct Node;

public:
//!@name Static properties
//!@{
  //! The height needed to be able to hold at least numElements elements
  static constexpr size_t maxHeight = BTreeProperties::maxHeightBound(numElements, minDegree);

  //! The maximum number of nodes that the tree can hold
  static constexpr size_t maxNodes = BTreeProperties::maxNodesInTree(maxHeight, minDegree);

  //! The maximum number of values that the tree can hold
  static constexpr size_t maxElements = maxNodes * (2 * minDegree - 1);
//!@}


public:
  constexpr BTree() : _rootPtr {0} {
    _rootPtr = _newNode();
  }

//!@name Modification
//!@{
  /*! @brief Clears the tree.
   *
   * @complexity{@math{O(N)}}
   */
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

  /*! @brief Add a value to the tree.
   *
   * Inserts a new value into the tree. This value must not already be in the
   * tree.
   *
   * @complexity{@math{\Theta(t \log_t N)}}
   * @pre The tree does not contain @p value
   */
  constexpr void insert(const T& value) {
    unsigned r = _rootPtr;

    if(_getNode(r).isFull()) { // Root is full, must be split
      unsigned s = _newNode();

      _rootPtr = s;

      _getNode(s).children.push_back(r);
      _splitChild(s, 0);
      _insertNonFull(s, value);
    } else {
      _insertNonFull(r, value);
    }

    ++_count;
  }

  /*! @brief Remove a value from the tree.
   *
   * Deletes a value from the tree.
   *
   * @complexity{@math{\Theta(t \log_t N)}}
   * @pre The value must exist in the tree.
   */
  constexpr void remove(const T& value) {
    _delete(_rootPtr, value);

    // In case the root node is valueless but has a child, shrink the tree
    if(_getNode(_rootPtr).values.size() == 0 && !_getNode(_rootPtr).isLeaf()) {
      unsigned emptyRoot = _rootPtr;

      _rootPtr = _getNode(_rootPtr).children.front();

      _markNodeDeleted(emptyRoot);
    }

    --_count;
  }

//!@}

//!@name Information
//!@{
  /*! @brief Check whether a value exists in the tree.
   *
   * Check whether a value is stored in the tree.
   *
   * @complexity{@math{\Theta(t \log_t N)}}
   */
  PURITY_WEAK constexpr bool contains(const T& value) const {
    auto nodeIndexOptional = _search(_rootPtr, value);

    return nodeIndexOptional.hasValue();
  }

  /*! @brief Checks if the tree contains T. If so, returns a const-ref to it
   *
   * @note This function may seem nonsensical, but it is an important building
   * block for maps in which what not all that is stored is compared against in
   * the search.
   *
   * @complexity{@math{\Theta(t \log_t N)}}
   */
  PURITY_WEAK constexpr Optional<const T&> getOption(const T& value) const {
    auto nodeIndexOptional = _search(_rootPtr, value);

    if(!nodeIndexOptional.hasValue()) {
      return {};
    }

    auto valueLowerBound = lowerBound<T, LessThanComparator>(
      _getNode(nodeIndexOptional.value()).values.begin(),
      _getNode(nodeIndexOptional.value()).values.end(),
      value,
      _lt
    );

    return Optional<const T&> {*valueLowerBound};
  }

  //! Dumps a graphViz representation of the B-Tree.
  PURITY_WEAK std::string dumpGraphviz() const {
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
        for(unsigned i = 0; i < node.values.size(); ++i) {
          graph << node.values.at(i);
          if(i != node.values.size() - 1) {
            graph << "|";
          }
        }
      } else {
        graph << "<f0>|";
        for(unsigned i = 1; i < node.children.size(); ++i) {
          graph << node.values.at(i-1) << "|<f" << i << ">";
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

  /*! @brief Validates the state of the tree by DFS traversal. Throws if
   *   anything is off.
   */
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

  /*! @brief Returns the number of elements in the tree
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr unsigned size() const {
    return _count;
  }
//!@}

//!@name Iterators
//!@{
  //! Nonmodifiable in-order iteration
  class const_iterator {
  public:
  //!@name Member types
  //!@{
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = T;
    using difference_type = int;
    using pointer = const T* const;
    using reference = const T&;

    enum class InitializeAs : unsigned {
      Begin = 0,
      End = 1
    };
  //!@}

  //!@name Constructor
  //!@{
    constexpr const_iterator(
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
        // End const_iterator initialization

        while(!_getCurrentNode().isLeaf()) {
          _indexStack.push_back(
            2 * _getCurrentNode().values.size()
          );

          _nodeStack.push_back(
            _getCurrentNode().children.back()
          );
        }

        // past-the-end of indices
        _indexStack.push_back(_getCurrentNode().values.size());
      }
    }
  //!@}

  //!@name Special member functions
  //!@{
    const_iterator() = delete;
    constexpr const_iterator(const const_iterator& other)
      : _baseRef(other._baseRef),
        _leftMostNode(other._leftMostNode),
        _rightMostNode(other._rightMostNode),
        _nodeStack(other._nodeStack),
        _indexStack(other._indexStack)
    {}
    constexpr const_iterator(const_iterator&& other) noexcept = default;
    constexpr const_iterator& operator = (const const_iterator& other) {
      if(this->_baseRef != other._baseRef) {
        throw "Assigning BTree const_iterator to another base instance!";
      }

      _leftMostNode = other._leftMostNode;
      _rightMostNode = other._rightMostNode;
      _nodeStack = other._nodeStack;
      _indexStack = other._indexStack;

      return *this;
    }
    constexpr const_iterator& operator = (const_iterator&& other) noexcept = default;
    ~const_iterator() = default;
  //!@}

    //! Prefix increment
    constexpr const_iterator& operator ++ () {
      auto indexLimit = _currentNodeIndexLimit();

      if(_indexStack.back() == indexLimit) { // We are already the end const_iterator
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

          // Increment here, now we are placed on a value at an internal node
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

      // Now we are a leaf, and placed on the first value
      _indexStack.push_back(0);

      return *this;
    }

    //! Postfix increment
    constexpr const_iterator operator ++ (int) {
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }

    //! Prefix decrement
    constexpr const_iterator& operator -- () {
      if(_nodeStack.back() == _leftMostNode && _indexStack.back() == 0) {
        // We are the begin const_iterator
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
          2 * _getCurrentNode().values.size()
        );
        _nodeStack.push_back(
          _getCurrentNode().children.back()
        );
      }

      _indexStack.push_back(
         _getCurrentNode().values.size() - 1
      );

      return *this;
    }

    //! Postfix decrement
    constexpr const_iterator operator -- (int) {
      const_iterator retval = *this;
      --(*this);
      return retval;
    }

    //! Compares on basis of positional equality
    PURITY_WEAK constexpr bool operator == (const const_iterator& other) const {
      return (
        _nodeStack == other._nodeStack
        && _indexStack == other._indexStack
        && _leftMostNode == other._leftMostNode
        && _rightMostNode == other._rightMostNode
      );
    }

    PURITY_WEAK constexpr bool operator != (const const_iterator& other) const {
      return !(
        *this == other
      );
    }

    //! Non-modifiable access
    PURITY_WEAK constexpr reference operator *() const {
      if(_getCurrentNode().isLeaf()) {
        return _getCurrentNode().values.at(
          _indexStack.back()
        );
      }

      return _getCurrentNode().values.at(
        _indexStack.back() / 2
      );
    }

  private:
  //!@name State
  //!@{
    const BTree& _baseRef;
    const unsigned _leftMostNode;
    const unsigned _rightMostNode;
    DynamicArray<unsigned, maxHeight + 1> _nodeStack;
    DynamicArray<unsigned, maxHeight + 1> _indexStack;
  //!@}

  //!@name Private member functions
  //!@{
    PURITY_WEAK constexpr const Node& _getCurrentNode() const {
      return _baseRef._getNode(_nodeStack.back());
    }

    PURITY_WEAK constexpr unsigned _currentNodeIndexLimit() const {
      // For leaves, the past-the-end position is the size of values
      if(_getCurrentNode().isLeaf()) {
        return _getCurrentNode().values.size();
      }

      /* For internal nodes, the past-the-end position is the size of values
       * plus the size of children + 1
       */
      return 2 * _getCurrentNode().values.size() + 1;
    }
  //@}
  };
  //! Returns a const iterator to the first value in the tree
  PURITY_WEAK constexpr const_iterator begin() const {
    return const_iterator(
      *this,
      const_iterator::InitializeAs::Begin
    );
  }

  //! Returns a past-the-end const iterator
  PURITY_WEAK constexpr const_iterator end() const {
    return const_iterator(
      *this,
      const_iterator::InitializeAs::End
    );
  }
//!@}

//!@name Operators
//!@{
  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator < (const BTree& other) const {
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

  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator > (const BTree& other) const {
    return (other < *this);
  }

  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator == (const BTree& other) const {
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

  //! Lexicographical comparison
  PURITY_WEAK constexpr bool operator != (const BTree& other) const {
    return !(*this == other);
  }
//!@}

private:
//!@name Private types
//!@{
  //! Type for Nodes
  struct Node {
  //!@name Static properties
  //!@{
    static constexpr unsigned minSize = minDegree - 1;
    static constexpr unsigned maxSize = 2 * minDegree - 1;
  //!@}

  //!@name State
  //!@{
    DynamicArray<T, maxSize> values;
    DynamicArray<unsigned, maxSize + 1> children;
  //!@}

    //! Default constructor
    constexpr Node() = default;

  //!@name Information
  //!@{
    //! Returns whether the node is a leaf node (i.e. has no children)
    PURITY_WEAK constexpr bool isLeaf() const {
      return children.size() == 0;
    }

    //! Returns whether the node is full
    PURITY_WEAK constexpr bool isFull() const {
      return values.size() == maxSize;
    }
  //!@}
  };
//!@}

//!@name Private state
//!@{
  //! 'Pointer' to root node
  unsigned _rootPtr;

  //! Array holding all tree nodes
  DynamicArray<Node, maxNodes> _nodes;

  //! Array holding 'pointers' to any 'deleted' tree nodes
  DynamicArray<unsigned, maxNodes> _garbage;

  //! Less-than comparator instance
  LessThanComparator _lt;

  //! Equality comparator instance
  EqualityComparator _eq;

  //! Number of contained elements
  unsigned _count = 0;
//!@}

//!@name Private member functions
//!@{
  //! 'Allocates' a new node and returns a 'pointer' to it
  constexpr unsigned _newNode() {
    // If there are nodes in the garbage, take those first
    if(_garbage.size() > 0) {
      unsigned newNodeIndex = _garbage.back();
      _garbage.pop_back();

      // Refresh the node
      _nodes.at(newNodeIndex) = Node {};

      return newNodeIndex;
    }

    if(_nodes.size() == maxNodes) {
      // We cannot get any new nodes! Nothing in the garbage
      throw "The maximum number of nodes has been reached for a BTree!";
    }

    // Default: Just expand the dynamic array with a fresh Node
    _nodes.push_back(Node {});
    return _nodes.size() - 1;
  }

  //! Marks a node as 'deleted' for recycling in _newNode
  constexpr void _markNodeDeleted(unsigned nodeIndex) {
    _garbage.push_back(nodeIndex);
  }

  //! Fetch a modifiable node by its 'pointer'
  constexpr Node& _getNode(unsigned nodeIndex) {
    return _nodes.at(nodeIndex);
  }

  //! Fetch an unmodifiable node by its 'pointer'
  PURITY_WEAK constexpr const Node& _getNode(unsigned nodeIndex) const {
    return _nodes.at(nodeIndex);
  }

  //! Recursive search for an element in a subtree rooted at node
  PURITY_WEAK constexpr Optional<unsigned> _search(unsigned nodeIndex, const T& value) const {
    auto node = _getNode(nodeIndex);

    auto valueLowerBound = lowerBound<T, LessThanComparator>(
      node.values.begin(),
      node.values.end(),
      value,
      _lt
    );

    // In case the lower bound is actually our sought value, return this node
    if(valueLowerBound != node.values.end() && _eq(*valueLowerBound, value)) {
      return Optional<unsigned> {nodeIndex};
    }

    // If we haven't found the value and this node is a leaf, search fails
    if(node.isLeaf()) {
      return Optional<unsigned> {};
    }

    /* Otherwise descend to the child at the same index as the lower bound in
     * values (this also works in case LB is at values.end() since children.size()
     * == values.size() + 1) and look there
     */
    return _search(
      node.children.at(
        valueLowerBound - node.values.begin()
      ),
      value
    );
  }

  /*!
   * @brief If a child we want to descend into during insertion is full, we
   *   have to split it in order to be certain we can insert into it if
   *   necessary.
   */
  constexpr void _splitChild(unsigned nodeIndex, const unsigned i) {
    // i is the child index in node's values being split since that node is full
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

    // Move values
    right.values = left.values.splice(minDegree);

    // In case left is not a leaf, move the children too
    if(!left.isLeaf()) {
      right.children = left.children.splice(minDegree);
    }

    // Insert the original median value into the non-full node
    parent.values.insertAt(
      parent.values.begin() + i,
      left.values.back()
    );

    // Have to remove it from left, too
    left.values.pop_back();

    // And assign right as the child to the right of the inserted median value
    parent.children.insertAt(
      parent.children.begin() + i + 1,
      rightIndex
    );
  }

  //! Insert case if the node we can insert into isn't full yet
  constexpr void _insertNonFull(unsigned nodeIndex, const T& value) {
    auto& node = _getNode(nodeIndex);

    if(node.isLeaf()) {
      auto valueLowerBound = lowerBound<T, LessThanComparator>(
        node.values.begin(),
        node.values.end(),
        value,
        _lt
      );

      if(valueLowerBound != node.values.end() && _eq(*valueLowerBound, value)) {
        throw "That key already exists in the tree!";
      }

      node.values.insertAt(valueLowerBound, value);
    } else {
      // Where to go?
      auto valueLowerBound = lowerBound<T, LessThanComparator>(
        node.values.begin(),
        node.values.end(),
        value,
        _lt
      );

      auto childPos = valueLowerBound - node.values.begin();

      // In case the purported child is full, split it!
      if(_getNode(node.children.at(childPos)).isFull()) {
        _splitChild(nodeIndex, childPos);

        /* values has an additional value from the split, check if index has to
         * be incremented
         */
        if(_lt(node.values.at(childPos), value)) {
          ++childPos;
        }
      }

      // The target child cannot be full anymore, so we can call
      _insertNonFull(node.children.at(childPos), value);
    }
  }

  //! Returns whether a node is the root
  PURITY_WEAK constexpr bool _isRootNode(unsigned nodeIndex) const {
    return nodeIndex == _rootPtr;
  }

  //! Returns the smallest leaf node in the sub-tree rooted at nodePtr
  PURITY_WEAK constexpr unsigned _smallestLeafNode(unsigned nodeIndex) const {
    while(!_getNode(nodeIndex).isLeaf()) {
      nodeIndex = _getNode(nodeIndex).children.front();
    }

    return nodeIndex;
  }

  //! Returns the largest leaf node in the sub-tree rooted at nodePtr
  PURITY_WEAK constexpr unsigned _largestLeafNode(unsigned nodeIndex) const {
    while(!_getNode(nodeIndex).isLeaf()) {
      nodeIndex = _getNode(nodeIndex).children.back();
    }

    return nodeIndex;
  }

  //! Recursively deletes a value from a sub-tree rooted at node
  constexpr void _delete(unsigned nodeIndex, const T& value) {
    auto& node = _getNode(nodeIndex);

    auto valueLowerBound = lowerBound<T, LessThanComparator>(
      node.values.begin(),
      node.values.end(),
      value,
      _lt
    );

    unsigned indexOfLB = valueLowerBound - node.values.begin();

    if(valueLowerBound != node.values.end() && _eq(*valueLowerBound, value)) {
      // value to remove is in this node's values
      if(node.isLeaf()) {
        // Case 1
        node.values.removeAt(valueLowerBound);
      } else {
        // Case 2
        if(_getNode(node.children.at(indexOfLB)).values.size() >= minDegree) {
          // Case 2a: Predecessor of value is maximum in subtree to the left

          // Predecessor value is largest value of largest leaf node in left subtree
          T predecessor = _getNode(_largestLeafNode(
            node.children.at(indexOfLB)
          )).values.back();

          // Replace the value to be deleted by its predecessor
          *valueLowerBound = predecessor;

          // Recursively delete the predecessor
          _delete(node.children.at(indexOfLB), predecessor);
        } else if(_getNode(node.children.at(indexOfLB + 1)).values.size() >= minDegree) {
          // Case 2b: Successor of value is minimum in subtree to the right

          // The successor value is the leftmost / smallest one
          T successor = _getNode(_smallestLeafNode(
            node.children.at(indexOfLB + 1)
          )).values.front();

          // Replace the value to be deleted by its successor
          *valueLowerBound = successor;

          // Recursively delete the successor
          _delete(node.children.at(indexOfLB + 1), successor);
        } else {
          /* Case 2c: Merge the value to delete, all of the right child into the
           * left child. The current node loses both k and the pointer to the
           * right child
           */

          unsigned leftChildIndex = node.children.at(indexOfLB);
          auto& leftChild = _getNode(leftChildIndex);

          unsigned rightChildIndex = node.children.at(indexOfLB + 1);
          auto& rightChild = _getNode(rightChildIndex);

          // Add the value to the left child
          leftChild.values.push_back(value);

          // Merge the right child into the left child
          leftChild.values.copyIn(rightChild.values);
          if(!leftChild.isLeaf()) {
            leftChild.children.copyIn(rightChild.children);
          }

          // Remove the value and child pointer to rightChild from left
          node.values.removeAt(valueLowerBound);
          node.children.removeAt(
            node.children.begin() + indexOfLB + 1
          );

          // Delete the right child
          _markNodeDeleted(rightChildIndex);

          // Delete the value recursively from the left child
          _delete(leftChildIndex, value);
        }
      }
    } else {
      /* Case 3: The value to delete is not in this node's values and we have to
       * descend in the tree. Need to ensure that any node we descend to has
       * at least minDegree values!
       */
      unsigned targetChildIndex = node.children.at(indexOfLB);
      auto& targetChild = _getNode(targetChildIndex);

      if(targetChild.values.size() == minDegree - 1) {
        // Case 3a Move some values around from left or right siblings

        if(
          indexOfLB != 0
          && _getNode(node.children.at(indexOfLB - 1)).values.size() >= minDegree
        ) {
          unsigned leftSiblingIndex = node.children.at(indexOfLB - 1);
          auto& leftSibling = _getNode(leftSiblingIndex);

          // Move value at LB into targetChild
          targetChild.values.insertAt(
            targetChild.values.begin(),
            node.values.at(indexOfLB - 1)
          );

          // Last value of left sibling replaces value at LB
          node.values.at(indexOfLB - 1) = leftSibling.values.back();
          leftSibling.values.pop_back();

          // In case it is not a leaf, we move the child pointer too
          if(!targetChild.isLeaf()) {
            targetChild.children.insertAt(
              targetChild.children.begin(),
              leftSibling.children.back()
            );

            leftSibling.children.pop_back();
          }
        } else if(
          indexOfLB < node.values.size()
          && _getNode(node.children.at(indexOfLB + 1)).values.size() >= minDegree
        ) {
          unsigned rightSiblingIndex = node.children.at(indexOfLB + 1);
          auto& rightSibling = _getNode(rightSiblingIndex);

          // Move value at LB into targetChild
          targetChild.values.push_back(
            node.values.at(indexOfLB)
          );

          // First value of right sibling replaces value at LB
          *valueLowerBound = rightSibling.values.front();
          rightSibling.values.removeAt(
            rightSibling.values.begin()
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

            // Move value down to left sibling
            leftSibling.values.push_back(
              node.values.at(indexOfLB - 1)
            );
            node.values.removeAt(
              node.values.begin() + indexOfLB - 1
            );

            // Merge values and children of targetChild into leftSibling
            leftSibling.values.copyIn(targetChild.values);
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

            targetChild.values.push_back(
              node.values.at(indexOfLB)
            );
            node.values.removeAt(
              node.values.begin() + indexOfLB
            );

            targetChild.values.copyIn(rightSibling.values);
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

      _delete(node.children.at(indexOfLB), value);
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

    // A non-root node has min. t-1 values
    if(
      !_isRootNode(nodeIndex)
      && node.values.size() < minDegree - 1
    ) {
      throw "Not every internal node has min. t-1 values!";
    }

    // Every internal node with n values has n+1 children
    if(
      !_isRootNode(nodeIndex)
      && !node.isLeaf()
      && node.values.size() != node.children.size() - 1
    ) {
      throw "Not every internal node with n values has n+1 children!";
    }

    // Every value list is ordered
    if(!isTotallyOrdered(node.values)) {
      throw "Not all value lists are totally ordered!";
    }

    /* Children 'between' values have values that are within the interval set by
     * the parents
     */
    if(!node.isLeaf()) {
      for(unsigned i = 1; i < node.children.size(); ++i) {
        if(
          !_lt(
            _getNode(node.children.at(i - 1)).values.back(),
            node.values.at(i - 1)
          ) || !_lt(
            node.values.at(i - 1),
            _getNode(node.children.at(i)).values.front()
          )
        ) {
          throw "Not all children's values are bounded by the parent!";
        }
      }
    }
  }
//!@}
};

} // namespace temple

#endif
