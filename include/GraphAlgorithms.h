#ifndef INCLUDE_GRAPH_ALGORITHMS_H
#define INCLUDE_GRAPH_ALGORITHMS_H

#include "common_typedefs.h"
#include <boost/graph/breadth_first_search.hpp>

#include <iostream>
#include "Delib/ElementInfo.h"

/*! @file
 *
 * Contains a number of graph-level algorithms where connectivity alone is
 * relevant.
 */

namespace MoleculeManip {

namespace GraphAlgorithms {

//! BGL BFS Visitors
namespace BFSVisitors {

/*!
 * Basis on which algorithms with limited depth BFS iteration can be built. Is 
 * inexpensive for small depths in sparse graphs, but probably very space 
 * inefficient for large depths in dense graphs, where pre-filtering nodes is
 * (maybe) preferable -> benchmark!
 *
 * In its current form entirely unused, but kept as reference.
 */
class LimitedDepthBFSVisitor : public boost::default_bfs_visitor {
private:
  std::map<AtomIndexType, unsigned> _depthMap;
  const unsigned _maxDepth;

public:
  class EarlyExit {};

  LimitedDepthBFSVisitor(
    const AtomIndexType& startingFrom,
    const unsigned& maxDepth
  ) : _depthMap({
        {startingFrom, 0}
      }),
      _maxDepth(maxDepth)
  {}

  template<typename Graph>
  void discover_vertex(const AtomIndexType& u, const Graph& g) {
    assert(_depthMap.count(u));
    /*std::cout << "Discovered vertex " << Delib::ElementInfo::symbol(
      g[u].elementType
    ) << u << std::endl;*/
    if(_depthMap.at(u) > _maxDepth) {
      throw EarlyExit();
    }
  }

  template<typename Graph>
  void tree_edge(const EdgeIndexType& e, const Graph& g) {
    _depthMap[boost::target(e, g)] = _depthMap[boost::source(e, g)] + 1;
  }
};

class SubstituentLinkSearcher : public boost::default_bfs_visitor {
private:
  // Const members
  const std::set<AtomIndexType> _soughtIndices;
  const AtomIndexType _source;

  // State
  std::map<AtomIndexType, unsigned> _depthMap;
  std::map<AtomIndexType, AtomIndexType> _parentMap;

  // Side-effect output
  std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >& _connectedPairs;

public:
  SubstituentLinkSearcher(
    const AtomIndexType& source,
    const std::set<AtomIndexType>& soughtIndices,
    std::set<
      std::pair<AtomIndexType, AtomIndexType>
    >& connectedPairsOutput
  ) : _soughtIndices(soughtIndices),
      _source(source),
      _depthMap({
        {source, 0}
      }),
      _connectedPairs(connectedPairsOutput)
  {
    for(const auto& index : soughtIndices) {
      _parentMap[index] = index;
    }
  }

  template<typename Graph>
  void tree_edge(const EdgeIndexType& e, const Graph& g) {
    const auto source = boost::source(e, g);
    const auto target = boost::target(e, g);

    _depthMap[target] = _depthMap[source] + 1;

    if(source == _source) {
      return;
    }

    bool sourceInMap = (_parentMap.count(source) == 1);
    bool targetInMap = (_parentMap.count(target) == 1);

    if(sourceInMap && !targetInMap) {
      _parentMap[target] = _parentMap[source];
    } else if(targetInMap && !sourceInMap) {
      _parentMap[source] = _parentMap[target];
    } 
  }

  template<typename Graph>
  void non_tree_edge(const EdgeIndexType& e, const Graph& g) {
    const auto source = boost::source(e, g);
    const auto target = boost::target(e, g);

    bool sourceInMap = _parentMap.count(source) == 1;
    bool targetInMap = _parentMap.count(target) == 1;

    if(
      sourceInMap 
      && targetInMap 
      && _depthMap[target] == _depthMap[source] + 1
      && _parentMap[source] != _parentMap[target]
    ) {
      _connectedPairs.emplace(
        std::min(
          _parentMap[source],
          _parentMap[target]
        ),
        std::max(
          _parentMap[source],
          _parentMap[target]
        )
      );
    }
  }
};

/*!
 * BFS Visitor for better understanding the sequence of events in BGL BFS
 * traversal. Merely prints out a message to std::cout for each event.
 */
class ShowAllEventsBFSVisitor : public boost::default_bfs_visitor {
public:
  template<typename Graph>
  void discover_vertex(const AtomIndexType& u, const Graph& g) {
    std::cout << "Discovered vertex " << Delib::ElementInfo::symbol(
      g[u].elementType
    ) << u << std::endl;
  }

  template<typename Graph>
  void examine_vertex(const AtomIndexType& u, const Graph& g) {
    std::cout << "Examining vertex " << Delib::ElementInfo::symbol(
      g[u].elementType
    ) << u << std::endl;
  }

  template<typename Graph>
  void finish_vertex(const AtomIndexType& u, const Graph& g) {
    std::cout << "Finishing vertex " << Delib::ElementInfo::symbol(
      g[u].elementType
    ) << u << std::endl;
  }

  template<typename Graph>
  void examine_edge(const EdgeIndexType& e, const Graph& g) {
    std::cout << "Examining: " << boost::source(e, g) << " -> " << boost::target(e, g) << std::endl;
  }

  template<typename Graph>
  void tree_edge(const EdgeIndexType& e, const Graph& g) {
    std::cout << "Tree: " << boost::source(e, g) << " -> " << boost::target(e, g) << std::endl;
  }

  template<typename Graph>
  void non_tree_edge(const EdgeIndexType& e, const Graph& g) {
    std::cout << "Non-tree: " << boost::source(e, g) << " -> " << boost::target(e, g) << std::endl;
  }

  template<typename Graph>
  void gray_target(const EdgeIndexType& e, const Graph& g) {
    std::cout << boost::source(e, g) << " -> " << boost::target(e, g) 
      << " (gray)" << std::endl;
  }

  template<typename Graph>
  void black_target(const EdgeIndexType& e, const Graph& g) {
    std::cout << boost::source(e, g) << " -> " << boost::target(e, g) 
      << " (black)" << std::endl;
  }
};

} // namespace BFSVisitors

/*! Shorthand function for using the SubstituentLinkSearcher. Determines whether
 * substituents at a central atom are linked or not at the grpah level.
 */
std::set<
  std::pair<AtomIndexType, AtomIndexType>
> findSubstituentLinks(
  const GraphType& graph,
  const AtomIndexType& source,
  const std::set<AtomIndexType>& activeSubstituents
);

/*! 
 * Returns the number of connected components of the graph. This is a central
 * property as the library enforces this number to be always one for any given
 * Molecule. The data representation of a Molecule should not contain two 
 * disconnected graphs!
 */
unsigned numConnectedComponents(const GraphType& graph);

//! Data class to return removal safety information on the graph
struct RemovalSafetyData {
  std::set<EdgeIndexType> bridges;
  std::set<AtomIndexType> articulationVertices;
};

/*! 
 * Determines which vertices of the graph are articulation points and which
 * edges are bridges. Removing an articulation vertex would incur the removal of
 * a bridge, whose removal always disconnects a graph. Therefore, neither can
 * be removed without prior graph changes.
 */
RemovalSafetyData getRemovalSafetyData(const GraphType& graph);

} // namespace GraphAlgorithms

} // namespace MoleculeManip

#endif
