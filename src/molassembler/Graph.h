/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interface class for the molecular graph
 */

#ifndef INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H
#define INCLUDE_MOLASSEMBLER_OUTER_GRAPH_H

#include "boost/optional/optional_fwd.hpp"
#include "Utils/Geometry/ElementTypes.h"

#include "molassembler/Types.h"
#include "molassembler/IteratorRange.h"

#include <memory>
#include <vector>

namespace Scine {
// Forward-declarations
namespace Utils {
using ElementTypeCollection = std::vector<ElementType>;
class BondOrderCollection;
} // namespace Utils

namespace molassembler {

// Forward-declare PrivateGraph
class PrivateGraph;
class Cycles;

/**
 * @brief Represents the connectivity of atoms of a molecule
 *
 * An undirected graph consisting of atom vertices and bond edges. Vertices
 * store the atomic element type while vertices store a bond type that
 * distinguishes the bond orders one through six as well as a so-called eta
 * bond, which models connections between a central atom and a
 * haptically-bonded subset of atoms (i.e. a contiguous group of atoms all
 * bonded to a transition metal).
 *
 * The Graph class leaves a consumer a lot of freedom in the specification of
 * the molecule's graph, but does enforce some model limitations.
 * - A molecule's graph must consist of a single connected
 *   component, meaning that there must be a path from any atom of the molecule
 *   to any other.
 * - Single atoms are not considered molecules, so a molecule must consist of
 *   at least two mutually bonded atoms.
 * - Removing atoms or bonds from diatomic molecules are disallowed operations.
 * - Disconnecting a molecule into two logical molecules by removing a
 *   particular bond or atom is also disallowed.
 *
 * @note This class wraps PrivateGraph so that no Boost Graph types are exposed
 *   to library consumers.
 */
class MASM_EXPORT Graph {
public:
//!@name Member types
//!@{
  /*!
   * @brief Templated iterator facade based on PrivateGraph to provide iterative
   *   access to atom indices and edge indices
   *
   * @tparam T The type the iterator should yield
   * @tparam isVertexInitialized Whether this iterator accepts an atom index
   *   in its constructor.
   *
   * @note Although this is templated, there is no implementation in the header
   *   in order to hide the underlying Boost Graph types from library consumers.
   *   The required template specializations for the type defintions below are
   *   supplied in GraphIterators.cpp. Any other instantiations will fail.
   */
  template<typename T, bool isVertexInitialized>
  class InnerBasedIterator {
  public:
    static_assert(
      std::is_same<T, AtomIndex>::value || std::is_same<T, BondIndex>::value,
      "You may not instantiate this type for other Ts than AtomIndex or BondIndex"
    );

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

    /* We want to use this iterator facade in four ways:
     *
     * Iterate through all valid atom indices: (graph) -> AtomIndex
     * Iterate through all valid bond indices: (graph) -> BondIndex
     * Iterate through all atom indices adjacent to i: (graph, i) -> AtomIndex
     * Iterate through all bonds incident on i: (graph, i) -> BondIndex
     *
     * So, we template and enable-if constructors for the cases that the
     * template parameter isVertexInitialized (which denotes whether
     * additionally to the graph, a vertex descriptor is required in the
     * constructor):
     */

    /*!
     * @brief Construct an iterator from a graph and a boolean indicating
     *   begin/end
     *
     * @param inner The PrivateGraph (wrapper around BGL Types)
     * @param begin Whether this iterator denotes a begin or end iterator
     *
     * @note This constructor is enabled if the template parameter
     *   isVertexInitialized is false
     */
    template<bool Dependent = isVertexInitialized, std::enable_if_t<!Dependent, int>...>
    InnerBasedIterator(const PrivateGraph& inner, bool begin);

    /*!
     * @brief Construct an iterator from an atom index, a graph and a boolean
     *   indicating begin/end
     *
     * @param a An atom index around which adjacent vertices or incident edges
     *   are to be iterated over (depending on @p T)
     * @param inner The PrivateGraph (wrapper around BGL Types)
     * @param begin Whether this iterator denotes a begin or end iterator
     *
     * @note This constructor is enabled if the template parameter
     *   isVertexInitialized is true
     */
    template<bool Dependent = isVertexInitialized, std::enable_if_t<Dependent, int>...>
    InnerBasedIterator(AtomIndex a, const PrivateGraph& inner, bool begin);

    //! Prefix increment
    InnerBasedIterator& operator ++ ();
    //! Postfix increment
    InnerBasedIterator operator ++ (int);
    //! Dereference
    value_type operator * () const;

    //! Comparison operator
    bool operator == (const InnerBasedIterator& other) const;
    //! Inverts @see operator ==
    bool operator != (const InnerBasedIterator& other) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> pImpl_;
  };

  //! Iterator type yielding all valid atom indices
  using AtomIterator = InnerBasedIterator<AtomIndex, false>;
  //! Iterator type yielding all valid bond indices
  using BondIterator = InnerBasedIterator<BondIndex, false>;
  //! Iterator type yielding adjacent atoms to an atom
  using AdjacencyIterator = InnerBasedIterator<AtomIndex, true>;
  //! Iterator type yielding incident bonds to an atom
  using IncidentEdgesIterator = InnerBasedIterator<BondIndex, true>;
//!@}

//!@name Special member functions
//!@{
  Graph(Graph&& other) noexcept;
  Graph& operator = (Graph&& other) noexcept;
  Graph(const Graph& other);
  Graph& operator = (const Graph& other);
  ~Graph();
//!@{

//!@name Constructors
//!@{
  Graph();
  //! Wrapping constructor
  explicit Graph(PrivateGraph&& inner);
//!@}

//!@name Information
//!@{
  /*! @brief Returns whether two atoms are bonded
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool adjacent(AtomIndex a, AtomIndex b) const;
  /*! @brief Returns atoms matching an element type
   *
   * @complexity{@math{\Theta(N)}}
   */
  std::vector<AtomIndex> atomsOfElement(Utils::ElementType e) const;
  /*! @brief Optionally fetch the bond index of a possibly non-existent bond
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<BondIndex> bond(AtomIndex a, AtomIndex b) const;
  /*! @brief Generate a BondOrderCollection from the graph
   *
   * @complexity{@math{\Theta(B)}}
   */
  Utils::BondOrderCollection bondOrders() const;
  /*! @brief Fetch the bond type at a particular bond
   *
   * @complexity{@math{\Theta(1)}}
   */
  BondType bondType(const BondIndex& edge) const;
  /*! @brief Returns whether an atom can be removed without disconnecting the graph
   *
   * @complexity{@math{O(N)} worst case, if removal data is cached
   * @math{\Theta(1)}}
   *
   * @note This function is not thread-safe.
   */
  bool canRemove(AtomIndex a) const;
  /*! @brief Returns whether a bond can be removed without disconnecting the graph
   *
   * @complexity{@math{O(N)} worst case, if removal data is cached
   * @math{\Theta(1)}}
   *
   * @note This function is not thread-safe.
   */
  bool canRemove(const BondIndex& edge) const;
  /*! @brief Fetch a reference to Cycles
   *
   * @complexity{@math{O(B)} worst case where @math{B} is the number of bonds
   * in cycles, if cycles are cached @math{\Theta(1)}}
   *
   * @note This function is not thread-safe.
   */
  const Cycles& cycles() const;
  /*! @brief Returns the number of bonds incident upon an atom
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned degree(AtomIndex a) const;
  //! @brief Dumps a graphviz representation of the graph
  std::string dumpGraphviz() const;
  /*! @brief Fetch an element collection of all atoms
   *
   * @complexity{@math{\Theta(N)}}
   */
  Utils::ElementTypeCollection elementCollection() const;
  /*! @brief Fetch the element type of an atom
   *
   * @complexity{@math{\Theta(1)}}
   */
  Utils::ElementType elementType(AtomIndex a) const;

  /*! @brief Number of atoms in the graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  AtomIndex N() const;
  /*! @brief Number of bonds in the graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned B() const;

  /*! @brief Modular isomorphism comparison
   *
   * Returns None if the molecules are not isomorphic. Returns an index mapping
   * from this to other otherwise.
   */
  boost::optional<std::vector<AtomIndex>> modularIsomorphism(
    const Graph& other,
    AtomEnvironmentComponents components = AtomEnvironmentComponents::All
  ) const;

  /*! @brief Determine which vertices belong to which side of a bridge edge
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @note This function is not thread-safe.
   */
  std::pair<
    std::vector<AtomIndex>,
    std::vector<AtomIndex>
  > splitAlongBridge(BondIndex bridge) const;
//!@}

//!@name Operators
//!@{
  //! Full isomorphism comparison including element types and bond orders
  bool operator == (const Graph& other) const;
  inline bool operator != (const Graph& other) const {
    return !(*this == other);
  }
//!@}

//!@name Ranges
//!@{
  /*! @brief A begin-end pair of iterators that yield the range of valid atom
   *   indices.
   *
   * @complexity{@math{\Theta(N)}}
   */
  IteratorRange<AtomIterator> atoms() const;
  /*! @brief A begin-end pair of iterators that yield the range of valid bond
   *   indices
   *
   * @complexity{@math{\Theta(N)}}
   */
  IteratorRange<BondIterator> bonds() const;
  /*! @brief Fetch iterator pair yielding adjacents of an atom
   *
   * @param a The atom whose adjacents are desired
   *
   * @complexity{@math{\Theta(N)}}
   * @returns A begin-end pair of iterators that yield adjacent atoms of an atom
   */
  IteratorRange<AdjacencyIterator> adjacents(AtomIndex a) const;
  /*! @brief Fetch iterator pair yielding bonds indices indicent to an atom
   *
   * @param a The atom whose incident atoms are desired
   *
   * @complexity{@math{\Theta(N)}}
   * @returns A begin-end pair of iterators that yield incident bond indices of
   *   an atom
   */
  IteratorRange<IncidentEdgesIterator> bonds(AtomIndex a) const;
//!@}

  /*! @brief Access to library-internal graph representation class
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @warning This function is not intended for library consumers, merely used
   * for implementation purposes.
   */
  MASM_NO_EXPORT PrivateGraph& inner();

  /*! @brief Const-access to library-internal graph representation class
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @warning This function is not intended for library consumers, merely used
   * for implementation purposes.
   */
  MASM_NO_EXPORT const PrivateGraph& inner() const;

private:
  std::unique_ptr<PrivateGraph> innerPtr_;
};

} // namespace molassembler
} // namespace Scine

#endif
