/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Graph class to help in gradual discovery of ordering relations
 *
 * Implements a class that aids in the gradual discovery of ordering relations
 * for a set of values.
 *
 * @todo
 * - may be preferable to have OrderDiscoveryHelper emit pairs of Ts whose
 *   ordering relation is not yet known instead of iterating through the
 *   getUnorderedSets()
 *   particularly now since addAllFromOther and addRelationshipsFromOther
 *   exist, which can hypothetically create completely disjoint sets which can
 *   then get conflated in getSets and getUnorderedSets
 * - addTransferabilityEdges needs a better algorithm
 */

#ifndef INCLUDE_MOLASSEMBLER_ORDER_DISCORVERY_HELPER
#define INCLUDE_MOLASSEMBLER_ORDER_DISCORVERY_HELPER

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graphviz.hpp"
#include "molassembler/Temple/Adaptors/Transform.h"
#include "molassembler/Temple/Functional.h"

namespace Scine {

namespace molassembler {

/**
 * @brief Container aiding in gradual discovery of order
 * @tparam T The type whose order is to be discovered
 */
template<typename T>
class OrderDiscoveryHelper {
public:
  struct VertexData {
    T data;
  };

  using DependencyGraphType = boost::adjacency_list<
    /* OutEdgeListS = Type of Container for edges of a vertex
     * Options: vector, list, slist, set, multiset, unordered_set
     * Choice: setS
     */
    boost::setS,
    /* VertexListS = Type of Container for vertices
     * Options: vector, list, slist, set, multiset, unordered_set
     * Choice: vecS, removing vertices does not occur
     */
    boost::vecS,
    /* DirectedS = Is the graph directed or not?
     * Choice: Directed, do not need bidirectional
     */
    boost::directedS,
    /* VertexProperty = What information is stored about vertices?
     * Choice: Atom, containing an index and an element type
     */
    VertexData
  >;

  using VertexIndexType = typename DependencyGraphType::vertex_descriptor;

private:
  std::map<T, VertexIndexType> sourceMap_;
  DependencyGraphType graph_;

  /*!
   * @brief Get a grouped list of the ordered data sorted by out_degree ascending
   *
   * @complexity{@math{\Theta(N)}}
   */
  std::vector<
    std::vector<T>
  > getSetsByDegree_() const {
    /* Group all vertices's values by vertex out_degree
     *
     * The more out edges a vertex has, the "smaller" it is, since a less-than
     * relationship is represented by a directed edge from the smaller object
     * to the bigger object.
     */

    // Keep a mapping of out_degree to index in sets to avoid empty sets
    std::map<
      unsigned,
      std::vector<VertexIndexType>
    > degreeToSetMap;

    typename DependencyGraphType::vertex_iterator iter, end;
    std::tie(iter, end) = boost::vertices(graph_);

    while(iter != end) {
      auto outDegree = boost::out_degree(*iter, graph_);

      if(degreeToSetMap.count(outDegree) == 0) {
        degreeToSetMap[outDegree] = {*iter};
      } else {
        degreeToSetMap.at(outDegree).push_back(*iter);
      }

      ++iter;
    }

    /* The map is ordered by key (outDegree ASC), so the result of mapValues
     * should is ordered too, so we can just map the vertex descriptors for
     * every degree into the stored underlying data
     */
    return temple::map(
      degreeToSetMap,
      [&](const auto& mapIterPair) -> std::vector<T> {
        return temple::map_stl(
          mapIterPair.second,
          [&](const auto& index) -> T {
            return graph_[index].data;
          }
        );
      }
    );
  }

  VertexIndexType addItem_(const T& item) {
    auto newIndex = boost::add_vertex(graph_);

    graph_[newIndex].data = item;
    sourceMap_.emplace(
      item,
      newIndex
    );

    return newIndex;
  }


  class GraphvizWriter {
  private:
    const OrderDiscoveryHelper& baseRef_;

  public:
    GraphvizWriter(const OrderDiscoveryHelper& base) : baseRef_(base) {}

    void operator() (std::ostream& os) const {
      os << "graph [fontname = \"Arial\", layout = \"dot\"];\n"
        << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
        << "edge [fontname = \"Arial\"];\n";
    }

    void operator() (std::ostream& os, const VertexIndexType& vertexIndex) const {
      os << R"([label=")" << baseRef_.graph_[vertexIndex].data << R"("])";
    }

    void operator() (std::ostream& /* os */, const typename DependencyGraphType::edge_descriptor& /* edge */) const {
      // do nothing particular
    }
  };

public:
  OrderDiscoveryHelper() = default;

  /** @brief Arbitrary container constructor
   *
   * @tparam Container Container type
   * @param container Container instance
   *
   * @complexity{@math{\Theta(N)}}
   */
  template<typename Container>
  explicit OrderDiscoveryHelper(const Container& container) {
    setUnorderedValues(
      std::begin(container),
      std::end(container)
    );
  }

  /*! @brief Transfer relationships from another OrderDiscoveryHelper
   *
   * Adds any relationships from another OrderDiscoveryHelper that are not
   * yet present in this instance. No new vertices are added. Missing
   * transferability edges are added (if a < b && b < c, then a < c) If any
   * contradictory information is present, this function throws.
   *
   * @complexity{@math{\Theta(N^2)}}
   */
  void addRelationshipsFromOther(const OrderDiscoveryHelper& other) {
    // Build a mapping of vertices
    std::map<VertexIndexType, VertexIndexType> vertexMapping;

    for(VertexIndexType i = 0; i < boost::num_vertices(other.graph_); ++i) {
      // Does this vertex exist in this graph already?
      const auto& vertexData = other.graph_[i].data;
      if(sourceMap_.count(vertexData) > 0) {
        vertexMapping.emplace(
          i,
          sourceMap_.at(vertexData)
        );
      }
    }

    // Add non-contradicting edges from the other graph
    for(
      auto edgeIterPair = boost::edges(other.graph_);
      edgeIterPair.first != edgeIterPair.second;
      ++edgeIterPair.first
    ) {
      auto otherEdgeSource = boost::source(*edgeIterPair.first, other.graph_);
      auto otherEdgeTarget = boost::target(*edgeIterPair.first, other.graph_);

      if( // Both endpoints must have a counterpart in this graph
        vertexMapping.count(otherEdgeSource) > 0
        && vertexMapping.count(otherEdgeTarget) > 0
      ) {
        auto thisEdgeSource = vertexMapping.at(otherEdgeSource);
        auto thisEdgeTarget = vertexMapping.at(otherEdgeTarget);

        // Check if that edge already exists (or its inverse)
        auto thisGraphEdge = boost::edge(
          thisEdgeSource,
          thisEdgeTarget,
          graph_
        );

        // In case the graph edge is found, go to the next edge
        if(thisGraphEdge.second) {
          continue;
        }

        // If not, check if the inverse is present in the graph
        auto inverseGraphEdge = boost::edge(
          thisEdgeTarget,
          thisEdgeSource,
          graph_
        );

        if(inverseGraphEdge.second) {
          throw "Contradicting information in other OrderDiscoveryHelper graph!";
        }

        // Add the edge to this graph
        boost::add_edge(thisEdgeSource, thisEdgeTarget, graph_);

      }
    }

    addTransferabilityEdges();
  }

  /*! @brief Adds vertices and edges from another instance, then adds transferability edges
   *
   * Adds all information present in another OrderDiscoveryHelper that are not
   * yet present in this one. Any vertices not present in this graph are added,
   * plus any missing relationship edges. If any contradictory information is
   * present, this function throws. Missing transferability edges are added too
   * (if a < b && b < c, then a < c).
   *
   * @complexity{@math{O(N^2)}}
   */
  void addAllFromOther(const OrderDiscoveryHelper& other) {
    // Left is the other graph, right is this graph
    std::map<VertexIndexType, VertexIndexType> vertexMapping;

    // Add all new vertices from the other graph
    for(VertexIndexType i = 0; i < boost::num_vertices(other.graph_); ++i) {
      // Does this vertex exist in this graph already?
      const auto& vertexData = other.graph_[i].data;
      if(sourceMap_.count(vertexData) == 0) {
        auto newIndex = addItem_(other.graph_[i].data);
        vertexMapping.emplace(
          i,
          newIndex
        );
      } else {
        vertexMapping.emplace(
          i,
          sourceMap_.at(vertexData)
        );
      }
    }

    // Add non-contradicting edges from the other graph
    for(
      auto edgeIterPair = boost::edges(other.graph_);
      edgeIterPair.first != edgeIterPair.second;
      ++edgeIterPair.first
    ) {
      auto otherEdgeSource = boost::source(*edgeIterPair.first, other.graph_);
      auto otherEdgeTarget = boost::target(*edgeIterPair.first, other.graph_);

      auto thisEdgeSource = vertexMapping.at(otherEdgeSource);
      auto thisEdgeTarget = vertexMapping.at(otherEdgeTarget);

      // Check if that edge already exists (or its inverse)
      auto thisGraphEdge = boost::edge(
        thisEdgeSource,
        thisEdgeTarget,
        graph_
      );

      // In case the graph edge is found, go to the next edge
      if(thisGraphEdge.second) {
        continue;
      }

      // If not, check if the inverse is present in the graph
      auto inverseGraphEdge = boost::edge(
        thisEdgeTarget,
        thisEdgeSource,
        graph_
      );

      if(inverseGraphEdge.second) {
        throw "Contradicting information in other OrderDiscoveryHelper graph!";
      }

      // Add the edge to this graph
      boost::add_edge(thisEdgeSource, thisEdgeTarget, graph_);
    }

    addTransferabilityEdges();
  }

  /*! @brief Adds missing transferability edges
   *
   * Adds any missing transferability edges (if a < b && b < c, then a < c).
   *
   * @complexity{@math{O(N^2)}}
   *
   * This function is awful in terms of complexity, it's worst case of O(N²).
   * There must be better ways of doing this.
   */
  void addTransferabilityEdges() {
    for(VertexIndexType i = 0; i < boost::num_vertices(graph_); ++i) {
      std::vector<VertexIndexType> seeds;

      for(
        auto edgeIterPair = boost::out_edges(i, graph_);
        edgeIterPair.first != edgeIterPair.second;
        ++edgeIterPair.first
      ) {
        seeds.emplace_back(
          boost::target(*edgeIterPair.first, graph_)
        );
      }

      do {
        std::vector<VertexIndexType> newSeeds;

        for(const auto& seed: seeds) {
          if(!boost::edge(i, seed, graph_).second) {
            boost::add_edge(i, seed, graph_);
          }

          for(
            auto edgeIterPair = boost::out_edges(seed, graph_);
            edgeIterPair.first != edgeIterPair.second;
            ++edgeIterPair.first
          ) {
            newSeeds.emplace_back(
              boost::target(*edgeIterPair.first, graph_)
            );
          }
        }

        seeds = std::move(newSeeds);
      } while(!seeds.empty());
    }
  }

  /** @brief Range-abstracted @see setUnorderedValues
   *
   * @complexity{\math{\Theta(N)}}
   * @todo no checks against duplicate values
   */
  template<typename It>
  void setUnorderedValues(It iter, const It end) {
    if(boost::num_vertices(graph_) > 0) {
      graph_.clear();
    }

    sourceMap_.clear();

    for(/* */; iter != end; ++iter) {
      addItem_(*iter);
    }
  }

  /*!  @brief Container-abstracted @see setUnorderedValues
   *
   * @complexity{\math{\Theta(N)}}
   */
  template<typename Container>
  void setUnorderedValues(Container&& container) {
    static_assert(
      std::is_same<
        temple::traits::getValueType<Container>,
        T
      >::value,
      "Container value must match OrderDiscoveryHelper's template argument T"
    );

    setUnorderedValues(
      std::begin(container),
      std::end(container)
    );
  }

  /*! @brief Get a list of sets (in descending order) as currently discovered
   *
   * @complexity{@math{\Theta(N)}}
   */
  std::vector<
    std::vector<T>
  > getSets() const {
    // Keep only sets with at least one member
    return getSetsByDegree_();
  }

  /*! @brief Returns grouped data (in descending order) whose internal order is undecided
   *
   * @complexity{@math{\Theta(N)}}
   */
  std::vector<
    std::vector<T>
  > getUndecidedSets() const {
    // Keep only sets with more than one member
    return temple::copy_if(
      getSetsByDegree_(),
      [](const auto& set) -> bool {
        return set.size() > 1;
      }
    );
  }

  /*! @brief Returns whether the total order has been discovered or not
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool isTotallyOrdered() const {
    auto N = boost::num_vertices(graph_);

    return(boost::num_edges(graph_) == N * (N - 1) / 2);
  }

  /*! @brief Add a less-than relationship to the graph
   *
   * @complexity{@math{\Theta(1)}}
   */
  void addLessThanRelationship(
    const T& a,
    const T& b
  ) {
    assert(a != b);

    /* Less-than relationship is represented as directed edge from vertex
     * storing a to vertex storing b
     */
    boost::add_edge(
      sourceMap_.at(a),
      sourceMap_.at(b),
      graph_
    );
  }

  /*! @brief Dumps a graphviz string for visualization of relationships
   *
   * @complexity{@math{\Theta(N)}}
   */
  std::string dumpGraphviz() const {
    GraphvizWriter propertyWriter(*this);

    std::stringstream ss;

    boost::write_graphviz(
      ss,
      graph_,
      propertyWriter,
      propertyWriter,
      propertyWriter
    );

    return ss.str();
  }

  /*! @brief Check if a value is smaller than another
   *
   * This returns true if and only if both values are part of the graph and
   * if the less-than relationship is present in the order given, i.e. a < b.
   *
   * In all other cases (either of the values isn't in the graph, there is no
   * ordering relationship or the ordering relationship is the other way
   * around), this function returns false.
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool isSmaller(const T& a, const T& b) const {
    if(sourceMap_.count(a) == 0 || sourceMap_.count(b) == 0) {
      return false;
    }

    // is a < b?
    auto smallerEdge = boost::edge(sourceMap_.at(a), sourceMap_.at(b), graph_);
    return smallerEdge.second;
  }
};

} // namespace molassembler

} // namespace Scine

#endif
