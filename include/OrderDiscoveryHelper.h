#ifndef INCLUDE_MOLECULE_MANIP_ORDER_DISCORVERY_HELPER
#define INCLUDE_MOLECULE_MANIP_ORDER_DISCORVERY_HELPER

#include "boost/graph/graphviz.hpp"
#include "template_magic/Containers.h"

#include "common_typedefs.h"

/*! @file
 *
 * Implements a class that aids in the gradual discovery of ordering relations
 * for a set of values.
 */

namespace MoleculeManip {

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

  //! Get a list of sets by degree
  std::vector<
    std::vector<T>
  > _getSetsByDegree() const {
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

    /* The map should be ordered by key, meaning by outDegree, so the result of
     * mapValues should be ordered, too
     */
    return TemplateMagic::mapValues(
      degreeToSetMap,
      [&](const auto& indexSet) -> std::vector<T> {
        return TemplateMagic::map(
          indexSet,
          [&](const auto& index) -> T {
            return _graph[index].data;
          }
        );
      }
    );
  }

  class GraphvizWriter {
  private:
    const OrderDiscoveryHelper& _baseRef;

  public:
    GraphvizWriter(const OrderDiscoveryHelper& base) : _baseRef(base) {}

    void operator() (std::ostream& os) const {
      os << "graph [fontname = \"Arial\"];\n"
        << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
        << "edge [fontname = \"Arial\"];\n";
    }

    void operator() (std::ostream& os, const VertexIndexType& vertexIndex) const {
      os << R"([label=")" << _baseRef._graph[vertexIndex].data << R"("])";
    }

    void operator() (std::ostream&, const typename DependencyGraphType::edge_descriptor&) const {
      // do nothing particular
    }
  };

public:
  OrderDiscoveryHelper() = default;

  explicit OrderDiscoveryHelper(const std::set<T>& unorderedValues) {
    setUnorderedValues(unorderedValues);
  }

  void setUnorderedValues(const std::set<T>& unorderedValues) {
    if(boost::num_vertices(_graph) > 0) {
      _graph.clear();
    }

    _sourceMap.clear();

    for(const auto& value: unorderedValues) {
      auto newIndex = boost::add_vertex(_graph);
      _graph[newIndex].data = value;
      _sourceMap.emplace(
        value,
        newIndex
      );
    }
  }

  //! Get a list of sets (in ascending order) as currently discovered
  std::vector<
    std::vector<T>
  > getSets() const {
    // Keep only sets with at least one member
    return _getSetsByDegree();
  }

  /*! Returns a list of sets (in ascending order) whose internal order is yet
   * undecided
   */
  std::vector<
    std::vector<T>
  > getUndecidedSets() const {
    // Keep only sets with more than one member
    return TemplateMagic::moveIf(
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

    boost::add_edge(
      _sourceMap.at(a),
      _sourceMap.at(b),
      _graph
    );
  }

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
};

} // namespace MoleculeManip

#endif
