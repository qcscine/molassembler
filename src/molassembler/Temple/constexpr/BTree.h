/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief constexpr BTree
 *
 * Implements a constexpr BTree. Can be used for an ordered set-like container
 * with good complexity guarantees and turned into an associative container
 * with some comparator tweaking.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_BTREE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_BTREE_H

#include "molassembler/Temple/constexpr/Containers.h"
#include "molassembler/Temple/constexpr/DynamicArray.h"
#include "molassembler/Temple/constexpr/Math.h"
#include "molassembler/Temple/constexpr/Optional.h"

#include <sstream>

namespace Scine {
namespace Temple {
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
  constexpr BTree() : rootPtr_ {0} {
    rootPtr_ = newNode_();
  }

//!@name Modification
//!@{
  /*! @brief Clears the tree.
   *
   * @complexity{@math{O(N)}}
   */
  constexpr void clear() {
    // Refresh all nodes
    for(auto& node: nodes_) {
      node = Node {};
    }

    // Clear the garbage
    garbage_.clear();

    // Assign a new root
    rootPtr_ = newNode_();
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
    unsigned r = rootPtr_;

    if(getNode_(r).isFull()) { // Root is full, must be split
      unsigned s = newNode_();

      rootPtr_ = s;

      getNode_(s).children.push_back(r);
      splitChild_(s, 0);
      insertNonFull_(s, value);
    } else {
      insertNonFull_(r, value);
    }

    ++count_;
  }

  /*! @brief Remove a value from the tree.
   *
   * Deletes a value from the tree.
   *
   * @complexity{@math{\Theta(t \log_t N)}}
   * @pre The value must exist in the tree.
   */
  constexpr void remove(const T& value) {
    delete_(rootPtr_, value);

    // In case the root node is valueless but has a child, shrink the tree
    if(getNode_(rootPtr_).values.size() == 0 && !getNode_(rootPtr_).isLeaf()) {
      unsigned emptyRoot = rootPtr_;

      rootPtr_ = getNode_(rootPtr_).children.front();

      markNodeDeleted_(emptyRoot);
    }

    --count_;
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
    auto nodeIndexOptional = search_(rootPtr_, value);

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
    auto nodeIndexOptional = search_(rootPtr_, value);

    if(!nodeIndexOptional.hasValue()) {
      return {};
    }

    auto valueLowerBound = lowerBound<T, LessThanComparator>(
      getNode_(nodeIndexOptional.value()).values.begin(),
      getNode_(nodeIndexOptional.value()).values.end(),
      value,
      lt_
    );

    return Optional<const T&> {*valueLowerBound};
  }

  //! Dumps a graphViz representation of the B-Tree.
  PURITY_WEAK std::string dumpGraphviz() const {
    using namespace std::string_literals;

    std::stringstream graph;
    graph << "digraph g {\n"
      << "  node [shape=record, height=.1]\n\n";

    DynamicArray<unsigned, maxNodes> stack {rootPtr_};

    while(stack.size() > 0) {
      unsigned nodeIndex = stack.back();
      stack.pop_back();

      for(auto& childIndex : getNode_(nodeIndex).children) {
        stack.push_back(childIndex);
      }

      const auto& node = getNode_(nodeIndex);

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
    DynamicArray<unsigned, maxNodes> stack {rootPtr_};

    while(stack.size() > 0) {
      unsigned nodeIndex = stack.back();
      auto& node = getNode_(nodeIndex);
      stack.pop_back();

      for(auto& childIndex : node.children) {
        stack.push_back(childIndex);
      }

      validate_(nodeIndex);
    }
  }

  /*! @brief Returns the number of elements in the tree
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_WEAK constexpr unsigned size() const {
    return count_;
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
    ) : baseRef_(tree),
        leftMostNode_(tree.smallestLeafNode_(tree.rootPtr_)),
        rightMostNode_(tree.largestLeafNode_(tree.rootPtr_)),
        nodeStack_ {tree.rootPtr_}
    {
      if(!static_cast<bool>(static_cast<unsigned>(initDecision))) {
        // 0 is Begin, so not-false is begin

        while(!getCurrentNode_().isLeaf()) {
          indexStack_.push_back(0);

          nodeStack_.push_back(
            getCurrentNode_().children.front()
          );
        }

        // At pos 0 of the leaf indices
        indexStack_.push_back(0);
      } else {
        // End const_iterator initialization

        while(!getCurrentNode_().isLeaf()) {
          indexStack_.push_back(
            2 * getCurrentNode_().values.size()
          );

          nodeStack_.push_back(
            getCurrentNode_().children.back()
          );
        }

        // past-the-end of indices
        indexStack_.push_back(getCurrentNode_().values.size());
      }
    }
  //!@}

  //!@name Special member functions
  //!@{
    const_iterator() = delete;
    constexpr const_iterator(const const_iterator& other)
      : baseRef_(other.baseRef_),
        leftMostNode_(other.leftMostNode_),
        rightMostNode_(other.rightMostNode_),
        nodeStack_(other.nodeStack_),
        indexStack_(other.indexStack_)
    {}
    constexpr const_iterator(const_iterator&& other) noexcept = default;
    constexpr const_iterator& operator = (const const_iterator& other) {
      if(this->baseRef_ != other.baseRef_) {
        throw "Assigning BTree const_iterator to another base instance!";
      }

      leftMostNode_ = other.leftMostNode_;
      rightMostNode_ = other.rightMostNode_;
      nodeStack_ = other.nodeStack_;
      indexStack_ = other.indexStack_;

      return *this;
    }
    constexpr const_iterator& operator = (const_iterator&& other) noexcept = default;
    ~const_iterator() = default;
  //!@}

    //! Prefix increment
    constexpr const_iterator& operator ++ () {
      auto indexLimit = currentNodeIndexLimit_();

      if(indexStack_.back() == indexLimit) { // We are already the end const_iterator
        // Do nothing and return immediately
        return *this;
      }

      // In case we are a leaf, increment and re-check
      if(getCurrentNode_().isLeaf()) {
        ++indexStack_.back();

        /* If we hit the index limit for the node and we're not the rightmost
         * node, we have to go up the tree and to the right
         */
        if(
          indexStack_.back() == indexLimit
          && nodeStack_.back() != rightMostNode_
        ) {
          // Unwind the stack until we are at an incrementable position
          do {
            indexStack_.pop_back();
            nodeStack_.pop_back();
          } while(indexStack_.back() >= currentNodeIndexLimit_() - 1);

          // Increment here, now we are placed on a value at an internal node
          ++indexStack_.back();
        }

        return *this;
      }

      // We are on an internal node, incrementing puts us on a child
      ++indexStack_.back();
      nodeStack_.push_back(
        getCurrentNode_().children.at(
          indexStack_.back() / 2 // children are at even indices
        )
      );

      while(!getCurrentNode_().isLeaf()) {
        indexStack_.push_back(0);
        nodeStack_.push_back(
          getCurrentNode_().children.front()
        );
      }

      // Now we are a leaf, and placed on the first value
      indexStack_.push_back(0);

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
      if(nodeStack_.back() == leftMostNode_ && indexStack_.back() == 0) {
        // We are the begin const_iterator
        return *this;
      }

      // In case we are a leaf, decrement and re-check
      if(getCurrentNode_().isLeaf()) {
        // If we are at zero, we have to find a decrementable position
        if(indexStack_.back() == 0) {
          // Unwind the stack until we can decrement
          do {
            indexStack_.pop_back();
            nodeStack_.pop_back();
          } while(indexStack_.back() == 0);
        }

        // Decrement and return
        --indexStack_.back();
        return *this;
      }

      // We are an internal node, decrementing puts us on a child
      --indexStack_.back();
      nodeStack_.push_back(
        getCurrentNode_().children.at(
          indexStack_.back() / 2
        )
      );

      while(!getCurrentNode_().isLeaf()) {
        indexStack_.push_back(
          2 * getCurrentNode_().values.size()
        );
        nodeStack_.push_back(
          getCurrentNode_().children.back()
        );
      }

      indexStack_.push_back(
         getCurrentNode_().values.size() - 1
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
        nodeStack_ == other.nodeStack_
        && indexStack_ == other.indexStack_
        && leftMostNode_ == other.leftMostNode_
        && rightMostNode_ == other.rightMostNode_
      );
    }

    PURITY_WEAK constexpr bool operator != (const const_iterator& other) const {
      return !(
        *this == other
      );
    }

    //! Non-modifiable access
    PURITY_WEAK constexpr reference operator *() const {
      if(getCurrentNode_().isLeaf()) {
        return getCurrentNode_().values.at(
          indexStack_.back()
        );
      }

      return getCurrentNode_().values.at(
        indexStack_.back() / 2
      );
    }

  private:
  //!@name State
  //!@{
    const BTree& baseRef_;
    const unsigned leftMostNode_;
    const unsigned rightMostNode_;
    DynamicArray<unsigned, maxHeight + 1> nodeStack_;
    DynamicArray<unsigned, maxHeight + 1> indexStack_;
  //!@}

  //!@name Private member functions
  //!@{
    PURITY_WEAK constexpr const Node& getCurrentNode_() const {
      return baseRef_.getNode_(nodeStack_.back());
    }

    PURITY_WEAK constexpr unsigned currentNodeIndexLimit_() const {
      // For leaves, the past-the-end position is the size of values
      if(getCurrentNode_().isLeaf()) {
        return getCurrentNode_().values.size();
      }

      /* For internal nodes, the past-the-end position is the size of values
       * plus the size of children + 1
       */
      return 2 * getCurrentNode_().values.size() + 1;
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
  unsigned rootPtr_;

  //! Array holding all tree nodes
  DynamicArray<Node, maxNodes> nodes_;

  //! Array holding 'pointers' to any 'deleted' tree nodes
  DynamicArray<unsigned, maxNodes> garbage_;

  //! Less-than comparator instance
  LessThanComparator lt_;

  //! Equality comparator instance
  EqualityComparator eq_;

  //! Number of contained elements
  unsigned count_ = 0;
//!@}

//!@name Private member functions
//!@{
  //! 'Allocates' a new node and returns a 'pointer' to it
  constexpr unsigned newNode_() {
    // If there are nodes in the garbage, take those first
    if(garbage_.size() > 0) {
      unsigned newNodeIndex = garbage_.back();
      garbage_.pop_back();

      // Refresh the node
      nodes_.at(newNodeIndex) = Node {};

      return newNodeIndex;
    }

    if(nodes_.size() == maxNodes) {
      // We cannot get any new nodes! Nothing in the garbage
      throw "The maximum number of nodes has been reached for a BTree!";
    }

    // Default: Just expand the dynamic array with a fresh Node
    nodes_.push_back(Node {});
    return nodes_.size() - 1;
  }

  //! Marks a node as 'deleted' for recycling in newNode_
  constexpr void markNodeDeleted_(unsigned nodeIndex) {
    garbage_.push_back(nodeIndex);
  }

  //! Fetch a modifiable node by its 'pointer'
  constexpr Node& getNode_(unsigned nodeIndex) {
    return nodes_.at(nodeIndex);
  }

  //! Fetch an unmodifiable node by its 'pointer'
  PURITY_WEAK constexpr const Node& getNode_(unsigned nodeIndex) const {
    return nodes_.at(nodeIndex);
  }

  //! Recursive search for an element in a subtree rooted at node
  PURITY_WEAK constexpr Optional<unsigned> search_(unsigned nodeIndex, const T& value) const {
    auto node = getNode_(nodeIndex);

    auto valueLowerBound = lowerBound<T, LessThanComparator>(
      node.values.begin(),
      node.values.end(),
      value,
      lt_
    );

    // In case the lower bound is actually our sought value, return this node
    if(valueLowerBound != node.values.end() && eq_(*valueLowerBound, value)) {
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
    return search_(
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
  constexpr void splitChild_(unsigned nodeIndex, const unsigned i) {
    // i is the child index in node's values being split since that node is full
    auto& parent = getNode_(nodeIndex);

    // The node being split is afterwards considered the "left" node
    auto& left = getNode_(
      getNode_(nodeIndex).children.at(i)
    );

    // Allocate a new "right" node
    if(nodes_.size() == maxNodes) {
      throw "Inserting into full BTree";
    }

    auto rightIndex = newNode_();
    auto& right = getNode_(rightIndex);

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
  constexpr void insertNonFull_(unsigned nodeIndex, const T& value) {
    auto& node = getNode_(nodeIndex);

    if(node.isLeaf()) {
      auto valueLowerBound = lowerBound<T, LessThanComparator>(
        node.values.begin(),
        node.values.end(),
        value,
        lt_
      );

      if(valueLowerBound != node.values.end() && eq_(*valueLowerBound, value)) {
        throw "That key already exists in the tree!";
      }

      node.values.insertAt(valueLowerBound, value);
    } else {
      // Where to go?
      auto valueLowerBound = lowerBound<T, LessThanComparator>(
        node.values.begin(),
        node.values.end(),
        value,
        lt_
      );

      auto childPos = valueLowerBound - node.values.begin();

      // In case the purported child is full, split it!
      if(getNode_(node.children.at(childPos)).isFull()) {
        splitChild_(nodeIndex, childPos);

        /* values has an additional value from the split, check if index has to
         * be incremented
         */
        if(lt_(node.values.at(childPos), value)) {
          ++childPos;
        }
      }

      // The target child cannot be full anymore, so we can call
      insertNonFull_(node.children.at(childPos), value);
    }
  }

  //! Returns whether a node is the root
  PURITY_WEAK constexpr bool isRootNode_(unsigned nodeIndex) const {
    return nodeIndex == rootPtr_;
  }

  //! Returns the smallest leaf node in the sub-tree rooted at nodePtr
  PURITY_WEAK constexpr unsigned smallestLeafNode_(unsigned nodeIndex) const {
    while(!getNode_(nodeIndex).isLeaf()) {
      nodeIndex = getNode_(nodeIndex).children.front();
    }

    return nodeIndex;
  }

  //! Returns the largest leaf node in the sub-tree rooted at nodePtr
  PURITY_WEAK constexpr unsigned largestLeafNode_(unsigned nodeIndex) const {
    while(!getNode_(nodeIndex).isLeaf()) {
      nodeIndex = getNode_(nodeIndex).children.back();
    }

    return nodeIndex;
  }

  //! Recursively deletes a value from a sub-tree rooted at node
  constexpr void delete_(unsigned nodeIndex, const T& value) {
    auto& node = getNode_(nodeIndex);

    auto valueLowerBound = lowerBound<T, LessThanComparator>(
      node.values.begin(),
      node.values.end(),
      value,
      lt_
    );

    unsigned indexOfLB = valueLowerBound - node.values.begin();

    if(valueLowerBound != node.values.end() && eq_(*valueLowerBound, value)) {
      // value to remove is in this node's values
      if(node.isLeaf()) {
        // Case 1
        node.values.removeAt(valueLowerBound);
      } else {
        // Case 2
        if(getNode_(node.children.at(indexOfLB)).values.size() >= minDegree) {
          // Case 2a: Predecessor of value is maximum in subtree to the left

          // Predecessor value is largest value of largest leaf node in left subtree
          T predecessor = getNode_(largestLeafNode_(
            node.children.at(indexOfLB)
          )).values.back();

          // Replace the value to be deleted by its predecessor
          *valueLowerBound = predecessor;

          // Recursively delete the predecessor
          delete_(node.children.at(indexOfLB), predecessor);
        } else if(getNode_(node.children.at(indexOfLB + 1)).values.size() >= minDegree) {
          // Case 2b: Successor of value is minimum in subtree to the right

          // The successor value is the leftmost / smallest one
          T successor = getNode_(smallestLeafNode_(
            node.children.at(indexOfLB + 1)
          )).values.front();

          // Replace the value to be deleted by its successor
          *valueLowerBound = successor;

          // Recursively delete the successor
          delete_(node.children.at(indexOfLB + 1), successor);
        } else {
          /* Case 2c: Merge the value to delete, all of the right child into the
           * left child. The current node loses both k and the pointer to the
           * right child
           */

          unsigned leftChildIndex = node.children.at(indexOfLB);
          auto& leftChild = getNode_(leftChildIndex);

          unsigned rightChildIndex = node.children.at(indexOfLB + 1);
          auto& rightChild = getNode_(rightChildIndex);

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
          markNodeDeleted_(rightChildIndex);

          // Delete the value recursively from the left child
          delete_(leftChildIndex, value);
        }
      }
    } else {
      /* Case 3: The value to delete is not in this node's values and we have to
       * descend in the tree. Need to ensure that any node we descend to has
       * at least minDegree values!
       */
      unsigned targetChildIndex = node.children.at(indexOfLB);
      auto& targetChild = getNode_(targetChildIndex);

      if(targetChild.values.size() == minDegree - 1) {
        // Case 3a Move some values around from left or right siblings

        if(
          indexOfLB != 0
          && getNode_(node.children.at(indexOfLB - 1)).values.size() >= minDegree
        ) {
          unsigned leftSiblingIndex = node.children.at(indexOfLB - 1);
          auto& leftSibling = getNode_(leftSiblingIndex);

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
          && getNode_(node.children.at(indexOfLB + 1)).values.size() >= minDegree
        ) {
          unsigned rightSiblingIndex = node.children.at(indexOfLB + 1);
          auto& rightSibling = getNode_(rightSiblingIndex);

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
            auto& leftSibling = getNode_(leftSiblingIndex);

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

            markNodeDeleted_(targetChildIndex);

            --indexOfLB;
          } else { // Merge with right sibling
            unsigned rightSiblingIndex = node.children.at(indexOfLB + 1);
            auto& rightSibling = getNode_(rightSiblingIndex);

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

            markNodeDeleted_(rightSiblingIndex);
          }
        }
      }

      delete_(node.children.at(indexOfLB), value);
    }
  }

  //! Checks whether the node is a valid B-Tree node, and throws if anything is off
  constexpr void validate_(unsigned nodeIndex) const {
    // The node should not be in the garbage
    auto foundIter = garbage_.begin();

    while(foundIter != garbage_.end()) {
      if(*foundIter == nodeIndex) {
        break;
      }

      ++foundIter;
    }

    // Skip if the current node is in the garbage
    if(foundIter != garbage_.end()) {
      throw "An active node is marked as garbage!";
    }

    auto& node = getNode_(nodeIndex);

    // A non-root node has min. t-1 values
    if(
      !isRootNode_(nodeIndex)
      && node.values.size() < minDegree - 1
    ) {
      throw "Not every internal node has min. t-1 values!";
    }

    // Every internal node with n values has n+1 children
    if(
      !isRootNode_(nodeIndex)
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
          !lt_(
            getNode_(node.children.at(i - 1)).values.back(),
            node.values.at(i - 1)
          ) || !lt_(
            node.values.at(i - 1),
            getNode_(node.children.at(i)).values.front()
          )
        ) {
          throw "Not all children's values are bounded by the parent!";
        }
      }
    }
  }
//!@}
};

} // namespace Temple
} // namespace Scine

#endif
