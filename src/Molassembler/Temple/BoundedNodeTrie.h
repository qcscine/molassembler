/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Prefix trie of bounded values that tracks sub-tree fullness.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BOUNDED_NODE_TRIE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BOUNDED_NODE_TRIE_H

#include <memory>
#include <cassert>

#include "Molassembler/Temple/Functional.h"

#include "boost/dynamic_bitset.hpp"

namespace Scine {
namespace Molassembler {
namespace Temple {

/*!
 * @brief Data structure to store chains of discrete choices with an finite range
 *   of choices at each position.
 *
 * @tparam ChoiceIndex Value type in which the count of possible choices is
 *   representable. This needs to be an integer type.
 *
 * @note Complexiy annotations use @math{N} for the length of the set ChoiceList
 */
template<typename ChoiceIndex>
struct BoundedNodeTrie {
  static_assert(std::is_integral<ChoiceIndex>::value, "ChoiceIndex must be an integral type");

public:
//!@name Public types
//!@{
  //! Type used to represent choices made at each level of the tree
  using ChoiceList = std::vector<ChoiceIndex>;
  /*!
   * @brief Function signature needed to choose which child to descend to at a
   *   level in the tree
   */
  using ChoosingFunction = std::function<ChoiceIndex(const std::vector<ChoiceIndex>&, const boost::dynamic_bitset<>&)>;
//!@}

//!@name Constructors
//!@{
  BoundedNodeTrie() = default;
  /**
   * @brief Construct an empty trie, setting the upper exclusive bound on
   *   values at each level
   *
   * @param bounds The upper exclusive bound on values at each level.
   */
  BoundedNodeTrie(ChoiceList bounds) {
    setBounds(std::move(bounds));
  }
//!@}

//!@name Modification
//!@{
  /**
   * @brief Inserts a value list into the prefix trie if not already present
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @param values The value list to insert into the prefix trie
   *
   * @throws std::logic_error If no bounds are set.
   *
   * @returns Whether anything was inserted. I.e. returns true if @p values
   *   was not yet a member of the set.
   */
  bool insert(const ChoiceList& values) {
    // Ensure the bounds are set
    if(bounds_.empty()) {
      throw std::logic_error("No bounds are set!");
    }

    // Make sure we have a root to insert to
    if(!root_) {
      establishRoot_();
    }

    InsertResult result = root_->insert(values, bounds_, 0);

    if(result.insertedSomething) {
      ++size_;
    }

    return result.insertedSomething;
  }

  /**
   * @brief Create a new entry by choosing branch at each level of the trie
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @param chooseFunction Function that chooses the branch to descend to at
   *   each level. The function is supplied a list of viable choices at the
   *   depth (indicating which parts of the tree are full) and a bitset
   *   indicating which children exist.
   *
   * @throws std::logic_error If the trie is full, i.e. size() == capacity()
   *
   * @return The newly generated entry
   */
  ChoiceList generateNewEntry(const ChoosingFunction& chooseFunction) {
    if(bounds_.empty()) {
      throw std::logic_error("No bounds are set!");
    }

    if(!root_) {
      establishRoot_();
    }

    if(size_ == capacity_) {
      throw std::logic_error("Trie is at capacity!");
    }

    ChoiceList values;
    values.reserve(bounds_.size());

    InsertResult result = root_->generate(chooseFunction, values, bounds_, 0);

    if(!result.insertedSomething) {
      throw std::logic_error("Trie has failed to generate a new entry");
    }

    ++size_;

    return values;
  }

  /*! @brief Changes the bounds of the trie. Clears the trie.
   *
   * @complexity{@math{\Theta(N)}}
   */
  void setBounds(ChoiceList bounds) {
    bounds_ = std::move(bounds);

    // Make sure there is at least a single bound
    assert(!bounds_.empty());
    // For there to be a decision, there need to be at least two options
    assert(
      Temple::all_of(
        bounds_,
        [](const ChoiceIndex value) -> bool {
          return value > 1;
        }
      )
    );

    establishCapacity_();
    clear();
  }

  /*! @brief Removes all inserted lists from the trie
   *
   * @complexity{@math{\Theta(N)}}
   */
  void clear() {
    if(root_) {
      root_.reset();
    }

    size_ = 0;
  }
//!@}

//!@name Information
//!@{
  /**
   * @brief Checks whether a value list is present in the trie
   *
   * @complexity{@math{O(N)}}
   *
   * @note If you are planning to insert @p values if @p contains returns false,
   *   just use insert. That way you have only one trie traversal in all cases,
   *   and you get the information of whether it was already in the trie or not
   *   too.
   *
   * @param values The value list to check for
   * @return Whether the value list is in the tree
   */
  bool contains(const ChoiceList& values) const {
    if(!root_) {
      return false;
    }

    return root_->contains(values, 0);
  }

  //! Accessor for underlying bounds
  const ChoiceList& bounds() const {
    return bounds_;
  }

  /*! @brief Returns the number of value lists this trie contains
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned size() const {
    if(root_) {
      assert(root_->size() == size_);
    }

    return size_;
  }

  /*! @brief Returns the total number of value lists this trie might contain
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned capacity() const {
    return capacity_;
  }
//!@}

private:
//!@name Private types
//!@{
  struct InsertResult {
    bool subtreeIsFull;
    bool insertedSomething;
  };

  //! Abstract base class for nodes to allow differentiation
  struct AbstractNode {
    virtual ~AbstractNode() = default;

    virtual bool contains(const ChoiceList& values, unsigned depth) const = 0;

    virtual InsertResult insert(
      const ChoiceList& values,
      const ChoiceList& bounds,
      unsigned depth
    ) = 0;

    virtual InsertResult generate(
      const ChoosingFunction& choosingFunction,
      ChoiceList& values,
      const ChoiceList& bounds,
      unsigned depth
    ) = 0;

    virtual unsigned size() const = 0;
  };

  //! Owning pointer to base class
  using NodePtr = std::unique_ptr<AbstractNode>;

  //! Most common type of nodes have other AbstractNodes as children.
  class Node final : public AbstractNode {
  public:
    Node(const ChoiceIndex U) : children(U), fullChildren(U) {
      assert(U > 1);
    }

    bool contains(const ChoiceList& values, const unsigned depth) const final {
      const ChoiceIndex& branch = values.at(depth);
      const NodePtr& childPtr = children.at(branch);

      if(!childPtr) {
        return false;
      }

      return childPtr->contains(values, depth + 1);
    }

    InsertResult insert(const ChoiceList& values, const ChoiceList& bounds, const unsigned depth) final {
      const ChoiceIndex& branchChoice = values.at(depth);
      NodePtr& childPtr = children.at(branchChoice);

      InsertResult subtreeResult;

      bool insertedSomething = false;
      if(!childPtr) {
        if(depth + 1 == bounds.size() - 1) {
          childPtr = std::make_unique<Leaf>(bounds.at(depth + 1));
        } else {
          childPtr = std::make_unique<Node>(bounds.at(depth + 1));
        }

        insertedSomething = true;
      }

      subtreeResult = childPtr->insert(values, bounds, depth + 1);
      subtreeResult.insertedSomething |= insertedSomething;

      if(subtreeResult.subtreeIsFull) {
        fullChildren.set(branchChoice);
      }

      subtreeResult.subtreeIsFull = fullChildren.all();

      return subtreeResult;
    }

    InsertResult generate(
      const ChoosingFunction& choosingFunction,
      ChoiceList& values,
      const ChoiceList& bounds,
      const unsigned depth
    ) final {
      // Prepare parameters for the choosing function
      const unsigned N = children.size();
      boost::dynamic_bitset<> existingChildren(N);
      for(unsigned i = 0; i < N; ++i) {
        if(children[i]) {
          existingChildren.set(i);
        }
      }

      std::vector<ChoiceIndex> viableIndices;
      viableIndices.reserve(N);
      for(unsigned i = 0; i < N; ++i) {
        if(!fullChildren.test(i)) {
          viableIndices.push_back(i);
        }
      }

      const ChoiceIndex choice = choosingFunction(viableIndices, existingChildren);
      assert(choice < N);
      assert(Temple::makeContainsPredicate(viableIndices)(choice));

      values.push_back(choice);
      NodePtr& childPtr = children.at(choice);

      bool insertedSomething = false;
      if(!childPtr) {
        insertedSomething = true;
        if(depth + 1 == bounds.size() - 1) {
          childPtr = std::make_unique<Leaf>(bounds.at(depth + 1));
        } else {
          childPtr = std::make_unique<Node>(bounds.at(depth + 1));
        }
      }

      InsertResult subtreeResult = childPtr->generate(choosingFunction, values, bounds, depth + 1);

      subtreeResult.insertedSomething |= insertedSomething;

      if(subtreeResult.subtreeIsFull) {
        fullChildren.set(choice);
      }

      subtreeResult.subtreeIsFull = fullChildren.all();

      return subtreeResult;
    }

    unsigned size() const final {
      unsigned count = 0;

      for(const NodePtr& nodePtr : children) {
        if(nodePtr) {
          count += nodePtr->size();
        }
      }

      return count;
    }

  private:
    std::vector<NodePtr> children;
    boost::dynamic_bitset<> fullChildren;
  };

  //! Leaf node at bottom of the tree
  class Leaf final : public AbstractNode {
  public:
    Leaf(const ChoiceIndex U) : children(U) {
      assert(U > 1);
    }

    bool contains(const ChoiceList& values, const unsigned depth) const final {
      const ChoiceIndex& choice = values.at(depth);
      return children.test(choice);
    }

    InsertResult insert(const ChoiceList& values, const ChoiceList& /* bounds */, const unsigned depth) final {
      const ChoiceIndex& choice = values.at(depth);
      assert(choice < children.size());

      InsertResult result;
      result.insertedSomething = !children.test(choice);
      result.subtreeIsFull = children.all();

      children.set(choice);

      return result;
    }

    InsertResult generate(
      const ChoosingFunction& choosingFunction,
      ChoiceList& values,
      const ChoiceList& /* bounds */,
      const unsigned /* depth */
    ) final {
      std::vector<ChoiceIndex> viableIndices;
      const ChoiceIndex N = children.size();
      viableIndices.reserve(N);
      for(ChoiceIndex i = 0; i < N; ++i) {
        if(!children.test(i)) {
          viableIndices.push_back(i);
        }
      }

      const ChoiceIndex choice = choosingFunction(viableIndices, children);
      assert(choice < N);
      assert(Temple::makeContainsPredicate(viableIndices)(choice));

      InsertResult result;
      result.insertedSomething = !children.test(choice);

      values.push_back(choice);
      children.set(choice);

      result.subtreeIsFull = children.all();

      return result;
    }

    unsigned size() const final {
      return children.count();
    }

  private:
    boost::dynamic_bitset<> children;
  };
//!@}

//!@name Private member functions
//!@{
  //! Generate a root node
  void establishRoot_() {
    if(bounds_.size() == 1) {
      root_ = std::make_unique<Leaf>(bounds_.front());
    } else {
      root_ = std::make_unique<Node>(bounds_.front());
    }
  }

  void establishCapacity_() {
    capacity_ = Temple::accumulate(
      bounds_,
      1U,
      std::multiplies<>()
    );
  }
//!@}

//!@name Private members
//!@{
  ChoiceList bounds_;
  NodePtr root_;
  unsigned size_ = 0;
  unsigned capacity_ = 0;
//!@}
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
