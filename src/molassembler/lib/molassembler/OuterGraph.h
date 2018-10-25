// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H

#include "boost/optional/optional_fwd.hpp"
#include "Delib/ElementTypes.h"

#include "molassembler/Types.h"

#include <memory>

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

/*!@file
 *
 * @brief Interface class for the molecular graph
 *
 * An undirected graph consisting of atom vertices and bond edges. Vertices
 * store the atomic element type while vertices store a bond type that
 * distinguishes the bond orders one through six as well as a so-called eta
 * bond, which models connections between a central atom and a
 * haptically-bonded subset of atoms (i.e. a contiguous group of atoms all
 * bonded to a transition metal).
 *
 * The Graph class leaves a consumer a lot of freedom in the
 * specification of the molecule's graph, but does enforce some model
 * limitations.
 * - A molecule's graph must consist of a single connected
 *   component, meaning that there must be a path from any atom of the molecule
 *   to any other.
 * - Single atoms are not considered molecules, so a molecule must consist of
 *   at least two mutually bonded atoms.
 * - Removing atoms or bonds from diatomic molecules are disallowed operations.
 * - Disconnecting a molecule into two logical molecules by removing a
 *   particular bond or atom is also disallowed.
 */

// Forward-declarations
namespace Delib {
  class ElementTypeCollection;
} // namespace Delib

namespace molassembler {

// Forward-declare InnerGraph
class InnerGraph;
class Cycles;

class OuterGraph {
public:
//!@name Member types
//!@{
  template<typename T, bool isVertexInitialized>
  class InnerBasedIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = T;
    using reference = T;

    InnerBasedIterator(InnerBasedIterator&& other) noexcept;
    InnerBasedIterator& operator = (InnerBasedIterator&& other) noexcept;
    InnerBasedIterator(const InnerBasedIterator& other);
    InnerBasedIterator& operator = (const InnerBasedIterator& other);
    ~InnerBasedIterator();

    // Default constructor
    InnerBasedIterator();

    // Enable construction from Graph and bool if not vertex initialized
    template<bool Dependent = isVertexInitialized, std::enable_if_t<!Dependent, int>...>
    InnerBasedIterator(const InnerGraph& inner, bool begin);

    // Enable compound construction if vertex initialized
    template<bool Dependent = isVertexInitialized, std::enable_if_t<Dependent, int>...>
    InnerBasedIterator(AtomIndex a, const InnerGraph& inner, bool begin);

    InnerBasedIterator& operator ++ ();
    InnerBasedIterator operator ++ (int);
    value_type operator * () const;

    bool operator == (const InnerBasedIterator& other) const;
    bool operator != (const InnerBasedIterator& other) const;

  private:
    struct Impl;

#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
    std::experimental::propagate_const<
      std::unique_ptr<Impl>
    > _pImpl;
#else
    std::unique_ptr<Impl> _pImpl;
#endif
  };

  using AtomIterator = InnerBasedIterator<AtomIndex, false>;
  using BondIterator = InnerBasedIterator<BondIndex, false>;
  using AdjacencyIterator = InnerBasedIterator<AtomIndex, true>;
  using IncidentEdgesIterator = InnerBasedIterator<BondIndex, true>;

  template<typename Iter>
  using Range = std::pair<Iter, Iter>;
//!@}

//!@name Special member functions
//!@{
  OuterGraph(OuterGraph&& other) noexcept;
  OuterGraph& operator = (OuterGraph&& other) noexcept;
  OuterGraph(const OuterGraph& other);
  OuterGraph& operator = (const OuterGraph& other);
  ~OuterGraph();
//!@{

//!@name Constructors
//!@{
  OuterGraph();
  //! Wrapping constructor
  explicit OuterGraph(InnerGraph&& inner);
//!@}

//!@name Information
//!@{
  //! Returns whether two atoms are bonded
  bool adjacent(AtomIndex a, AtomIndex b) const;
  //! Optionally fetch the bond index of a possibly non-existent bond
  boost::optional<BondIndex> bond(AtomIndex a, AtomIndex b) const;
  //! Fetch the bond type at a particular bond
  BondType bondType(const BondIndex& edge) const;
  //! Returns whether an atom can be removed without disconnecting the graph
  bool canRemove(AtomIndex a) const;
  //! Returns whether a bond can be removed without disconnecting the graph
  bool canRemove(const BondIndex& edge) const;
  //! Fetch a reference to Cycles
  const Cycles& cycles() const;
  //! Returns the number of bonds incident upon an atom
  unsigned degree(AtomIndex a) const;
  //! Fetch an element collection of all atoms
  Delib::ElementTypeCollection elementCollection() const;
  //! Fetch the element type of an atom
  Delib::ElementType elementType(AtomIndex a) const;

  //! Number of atoms in the graph
  AtomIndex N() const;
  //! Number of bonds in the graph
  unsigned B() const;
//!@}

//!@name Ranges
//!@{
  //! A begin-end pair of iterators that yield the range of valid atom indices.
  Range<AtomIterator> atoms() const;
  //! A begin-end pair of iterators that yield the range of valid bond indices
  Range<BondIterator> bonds() const;
  /*!
   * @brief Fetch iterator pair yielding adjacents of an atom
   * @param a The atom whose adjacents are desired
   * @returns A begin-end pair of iterators that yield adjacent atoms of an atom
   */
  Range<AdjacencyIterator> adjacents(AtomIndex a) const;
  /*!
   * @brief Fetch iterator pair yielding bonds indices indicent to an atom
   * @param a The atom whose incident atoms are desired
   * @returns A begin-end pair of iterators that yield incident bond indices of
   *   an atom
   */
  Range<IncidentEdgesIterator> bonds(AtomIndex a) const;
//!@}

  /*! Access to library-internal graph representation class
   *
   * @warning This function is not intended for library consumers, merely used
   * for implementation purposes.
   */
  InnerGraph& inner() { return *_innerPtr; }

  /*! Const-access to library-internal graph representation class
   *
   * @warning This function is not intended for library consumers, merely used
   * for implementation purposes.
   */
  const InnerGraph& inner() const { return *_innerPtr; }

private:
#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
  std::experimental::propagate_const<
    std::unique_ptr<InnerGraph>
  > _innerPtr;

  mutable std::experimental::propagate_const<
    std::unique_ptr<Cycles>
  > _cachedCycles;
#else
  std::unique_ptr<InnerGraph> _innerPtr;
  mutable std::unique_ptr<Cycles> _cachedCycles;
#endif
};

} // namespace molassembler

#endif
