/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"

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
  std::map<T, VertexIndexType> _sourceMap;
  DependencyGraphType _graph;

  //! Get a grouped list of the ordered data sorted by out_degree ascending
  std::vector<
    std::vector<T>
  > _getSetsByDegree() const {
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
    std::tie(iter, end) = boost::vertices(_graph);

    while(iter != end) {
      auto outDegree = boost::out_degree(*iter, _graph);

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
            return _graph[index].data;
          }
        );
      }
    );
  }

  VertexIndexType _addItem(const T& item) {
    auto newIndex = boost::add_vertex(_graph);

    _graph[newIndex].data = item;
    _sourceMap.emplace(
      item,
      newIndex
    );

    return newIndex;
  }


  class GraphvizWriter {
  private:
    const OrderDiscoveryHelper& _baseRef;

  public:
    GraphvizWriter(const OrderDiscoveryHelper& base) : _baseRef(base) {}

    void operator() (std::ostream& os) const {
      os << "graph [fontname = \"Arial\", layout = \"dot\"];\n"
        << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
        << "edge [fontname = \"Arial\"];\n";
    }

    void operator() (std::ostream& os, const VertexIndexType& vertexIndex) const {
      os << R"([label=")" << _baseRef._graph[vertexIndex].data << R"("])";
    }

    void operator() (std::ostream& /* os */, const typename DependencyGraphType::edge_descriptor& /* edge */) const {
      // do nothing particular
    }
  };

public:
  OrderDiscoveryHelper() = default;

  template<typename Container>
  explicit OrderDiscoveryHelper(const Container& container) {
    setUnorderedValues(
      std::begin(container),
      std::end(container)
    );
  }

  /*!
   * Adds any relationships from another OrderDiscoveryHelper that are not
   * yet present in this one. No new vertices are added. Missing transferability
   * edges are added (if a < b && b < c, then a < c) If any contradictory
   * information is present, this function throws.
   */
  void addRelationshipsFromOther(const OrderDiscoveryHelper& other) {
    // Build a mapping of vertices
    std::map<VertexIndexType, VertexIndexType> vertexMapping;

    for(VertexIndexType i = 0; i < boost::num_vertices(other._graph); ++i) {
      // Does this vertex exist in this graph already?
      const auto& vertexData = other._graph[i].data;
      if(_sourceMap.count(vertexData) > 0) {
        vertexMapping.emplace(
          i,
          _sourceMap.at(vertexData)
        );
      }
    }

    // Add non-contradicting edges from the other graph
    for(
      auto edgeIterPair = boost::edges(other._graph);
      edgeIterPair.first != edgeIterPair.second;
      ++edgeIterPair.first
    ) {
      auto otherEdgeSource = boost::source(*edgeIterPair.first, other._graph);
      auto otherEdgeTarget = boost::target(*edgeIterPair.first, other._graph);

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
          _graph
        );

        // In case the graph edge is found, go to the next edge
        if(thisGraphEdge.second) {
          continue;
        }

        // If not, check if the inverse is present in the graph
        auto inverseGraphEdge = boost::edge(
          thisEdgeTarget,
          thisEdgeSource,
          _graph
        );

        if(inverseGraphEdge.second) {
          throw "Contradicting information in other OrderDiscoveryHelper graph!";
        }

        // Add the edge to this graph
        boost::add_edge(thisEdgeSource, thisEdgeTarget, _graph);

      }
    }

    addTransferabilityEdges();
  }

  /*!
   * Adds all information present in another OrderDiscoveryHelper that are not
   * yet present in this one. Any vertices not present in this graph are added,
   * plus any missing relationship edges. If any contradictory information is
   * present, this function throws. Missing transferability edges are added too
   * (if a < b && b < c, then a < c).
   */
  void addAllFromOther(const OrderDiscoveryHelper& other) {
    // Left is the other graph, right is this graph
    std::map<VertexIndexType, VertexIndexType> vertexMapping;

    // Add all new vertices from the other graph
    for(VertexIndexType i = 0; i < boost::num_vertices(other._graph); ++i) {
      // Does this vertex exist in this graph already?
      const auto& vertexData = other._graph[i].data;
      if(_sourceMap.count(vertexData) == 0) {
        auto newIndex = _addItem(other._graph[i].data);
        vertexMapping.emplace(
          i,
          newIndex
        );
      } else {
        vertexMapping.emplace(
          i,
          _sourceMap.at(vertexData)
        );
      }
    }

    // Add non-contradicting edges from the other graph
    for(
      auto edgeIterPair = boost::edges(other._graph);
      edgeIterPair.first != edgeIterPair.second;
      ++edgeIterPair.first
    ) {
      auto otherEdgeSource = boost::source(*edgeIterPair.first, other._graph);
      auto otherEdgeTarget = boost::target(*edgeIterPair.first, other._graph);

      auto thisEdgeSource = vertexMapping.at(otherEdgeSource);
      auto thisEdgeTarget = vertexMapping.at(otherEdgeTarget);

      // Check if that edge already exists (or its inverse)
      auto thisGraphEdge = boost::edge(
        thisEdgeSource,
        thisEdgeTarget,
        _graph
      );

      // In case the graph edge is found, go to the next edge
      if(thisGraphEdge.second) {
        continue;
      }

      // If not, check if the inverse is present in the graph
      auto inverseGraphEdge = boost::edge(
        thisEdgeTarget,
        thisEdgeSource,
        _graph
      );

      if(inverseGraphEdge.second) {
        throw "Contradicting information in other OrderDiscoveryHelper graph!";
      }

      // Add the edge to this graph
      boost::add_edge(thisEdgeSource, thisEdgeTarget, _graph);
    }

    addTransferabilityEdges();
  }

  /*!
   * Adds any missing transferability edges (if a < b && b < c, then a < c).
   * This function is awful in terms of complexity, it's worst case of O(N²).
   * There must be better ways of doing this.
   */
  void addTransferabilityEdges() {
    for(VertexIndexType i = 0; i < boost::num_vertices(_graph); ++i) {
      std::vector<VertexIndexType> seeds;

      for(
        auto edgeIterPair = boost::out_edges(i, _graph);
        edgeIterPair.first != edgeIterPair.second;
        ++edgeIterPair.first
      ) {
        seeds.emplace_back(
          boost::target(*edgeIterPair.first, _graph)
        );
      }

      do {
        std::vector<VertexIndexType> newSeeds;

        for(const auto& seed: seeds) {
          if(!boost::edge(i, seed, _graph).second) {
            boost::add_edge(i, seed, _graph);
          }

          for(
            auto edgeIterPair = boost::out_edges(seed, _graph);
            edgeIterPair.first != edgeIterPair.second;
            ++edgeIterPair.first
          ) {
            newSeeds.emplace_back(
              boost::target(*edgeIterPair.first, _graph)
            );
          }
        }

        seeds = std::move(newSeeds);
      } while(!seeds.empty());
    }
  }

  //! Sets which as-yet unordered values are to be investigated
  void setUnorderedValues(const std::set<T>& unorderedValues) {
    if(boost::num_vertices(_graph) > 0) {
      _graph.clear();
    }

    _sourceMap.clear();

    for(const auto& value: unorderedValues) {
      _addItem(value);
    }
  }

  template<typename It>
  void setUnorderedValues(It iter, const It end) {
    if(boost::num_vertices(_graph) > 0) {
      _graph.clear();
    }

    _sourceMap.clear();

    for(/* */; iter != end; ++iter) {
      _addItem(*iter);
    }
  }

  //! @todo no checks against duplicate values
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

  //! Get a list of sets (in descending order) as currently discovered
  std::vector<
    std::vector<T>
  > getSets() const {
    // Keep only sets with at least one member
    return _getSetsByDegree();
  }

  //! Returns grouped data (in descending order) whose internal order is undecided
  std::vector<
    std::vector<T>
  > getUndecidedSets() const {
    // Keep only sets with more than one member
    return temple::copy_if(
      _getSetsByDegree(),
      [](const auto& set) -> bool {
        return set.size() > 1;
      }
    );
  }

  //! Returns whether the total order has been discovered or not
  bool isTotallyOrdered() const {
    auto N = boost::num_vertices(_graph);

    return(boost::num_edges(_graph) == N * (N - 1) / 2);
  }

  //! Add a less-than relationship to the graph
  void addLessThanRelationship(
    const T& a,
    const T& b
  ) {
    assert(a != b);

    /* Less-than relationship is represented as directed edge from vertex
     * storing a to vertex storing b
     */
    boost::add_edge(
      _sourceMap.at(a),
      _sourceMap.at(b),
      _graph
    );
  }

  /*!
   * Dumps a graphviz string that allows the visualization of the internal
   * graph structure.
   */
  std::string dumpGraphviz() const {
    GraphvizWriter propertyWriter(*this);

    std::stringstream ss;

    boost::write_graphviz(
      ss,
      _graph,
      propertyWriter,
      propertyWriter,
      propertyWriter
    );

    return ss.str();
  }

  /*!
   * This returns true if and only if both values are part of the graph and
   * if the less-than relationship is present in the order given, i.e. a < b.
   *
   * In all other cases (either of the values isn't in the graph, there is no
   * ordering relationship or the ordering relationship is the other way
   * around), this function returns false.
   */
  bool isSmaller(const T& a, const T& b) const {
    if(_sourceMap.count(a) == 0 || _sourceMap.count(b) == 0) {
      return false;
    }

    // is a < b?
    auto smallerEdge = boost::edge(_sourceMap.at(a), _sourceMap.at(b), _graph);
    return smallerEdge.second;
  }
};

} // namespace molassembler

} // namespace Scine

#endif
