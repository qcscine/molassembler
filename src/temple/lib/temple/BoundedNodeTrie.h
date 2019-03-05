/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Prefix trie of bounded values that tracks sub-tree fullness.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BOUNDED_NODE_TRIE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BOUNDED_NODE_TRIE_H

#include <memory>
#include <vector>
#include <cassert>

#include "temple/Functional.h"

#include "boost/dynamic_bitset.hpp"

namespace temple {

/*!
 * @brief Data structure to store chains of discrete choices with an finite range
 *   of choices at each position.
 *
 * @tparam ChoiceIndex Value type in which the count of possible choices is
 *   representable. This needs to be an integer type.
 */
template<typename ChoiceIndex>
struct BoundedNodeTrie {
  static_assert(std::is_integral<ChoiceIndex>::value, "ChoiceIndex must be an integral type");

public:
//!@name Public types
//!@{
  using ChoiceList = std::vector<ChoiceIndex>;
  using ChoosingFunction = std::function<ChoiceIndex(const std::vector<ChoiceIndex>&, const boost::dynamic_bitset<>&)>;
//!@}

//!@name Special member functions
//!@{
  BoundedNodeTrie() = default;
  BoundedNodeTrie(ChoiceList bounds) {
    setBounds(std::move(bounds));
  }
//!@}

//!@name Modification
//!@{
  /**
   * @brief Inserts a value list into the prefix trie if not already present
   *
   * @param values The value list to insert into the prefix trie
   *
   * This has complexity linear in the size of the boundaries with which the
   * trie was constructed.
   *
   * @returns Whether anything was inserted. I.e. returns true if @p values
   *   was not yet a member of the set.
   */
  bool insert(const ChoiceList& values) {
    // Make sure we have a root to insert to
    if(!_root) {
      _establishRoot();
    }

    InsertResult result = _root->insert(values, _bounds, 0);

    if(result.insertedSomething) {
      ++_size;
    }

    return result.insertedSomething;
  }

  ChoiceList generateNewEntry(const ChoosingFunction& chooseFunction) {
    if(!_root) {
      _establishRoot();
    }

    if(_size == _capacity) {
      throw std::logic_error("Trie is at capacity!");
    }

    ChoiceList values;
    values.reserve(_bounds.size());

    InsertResult result = _root->generate(chooseFunction, values, _bounds, 0);

    assert(result.insertedSomething);

    ++_size;

    return values;
  }

  //! Changes the bounds of the trie. Clears the trie.
  void setBounds(ChoiceList bounds) {
    _bounds = std::move(bounds);

    // Make sure there is at least a single bound
    assert(!_bounds.empty());
    // For there to be a decision, there need to be at least two options
    assert(
      temple::all_of(
        bounds,
        [](const ChoiceIndex value) -> bool {
          return value > 1;
        }
      )
    );

    _establishCapacity();
    clear();
  }

  //! Removes all inserted lists from the trie
  void clear() {
    if(_root) {
      _root.reset();
    }

    _size = 0;
  }
//!@}

//!@name Information
//!@{
  /**
   * @brief Checks whether a value list is present in the trie
   *
   * This has complexity at worst linear in the length of the boundaries with
   * which the trie was constructed.
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
    if(!_root) {
      return false;
    }

    return _root->contains(values, 0);
  }

  //! Accessor for underlying bounds
  const ChoiceList& bounds() const {
    return _bounds;
  }

  /*!
   * @brief Returns the number of value lists this trie contains
   * @note This function traverses the entire trie to calculate the count.
   */
  unsigned size() const {
    if(_root) {
      assert(_root->size() == _size);
    }

    return _size;
  }

  //! Returns the total number of value lists this trie might contain
  unsigned capacity() const {
    return _capacity;
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
      assert(temple::makeContainsPredicate(viableIndices)(choice));

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
      assert(temple::makeContainsPredicate(viableIndices)(choice));

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
  void _establishRoot() {
    if(_bounds.size() == 1) {
      _root = std::make_unique<Leaf>(_bounds.front());
    } else {
      _root = std::make_unique<Node>(_bounds.front());
    }
  }

  void _establishCapacity() {
    _capacity = temple::accumulate(
      _bounds,
      1u,
      std::multiplies<>()
    );
  }
//!@}

//!@name Private members
//!@{
  ChoiceList _bounds;
  NodePtr _root;
  unsigned _size = 0;
  unsigned _capacity;
//!@}
};

} // namespace temple

#endif
