#ifndef INCLUDE_GRAPH_ALGORITHMS_H
#define INCLUDE_GRAPH_ALGORITHMS_H

#include "common_typedefs.h"
#include <boost/graph/breadth_first_search.hpp>

#include <iostream>
#include "Delib/ElementInfo.h"
#include "temple/Containers.h"
#include "temple/UnorderedSets.h"
#include "temple/TinySet.h"

/*! @file
 *
 * Contains a number of graph-level algorithms where connectivity alone is
 * relevant.
 */

namespace molassembler {

// Forward-declare CycleData
class CycleData;

//! Core graph-level algorithms (not requiring stereocenter information)
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

struct LinkInformation {
  //! An (asc) ordered pair of the substituent indices that are linked
  std::pair<AtomIndexType, AtomIndexType> indexPair;

  /*! The in-order atom sequence of the cycle atom indices.
   *
   * NOTE: The front and back indices are repeated.
   */
  std::vector<AtomIndexType> cycleSequence;

  /* TODO this operator isn't quite correct, in principle, the cycle sequence
   * is not the completely reduced form / has degrees of freedom, and a
   * lexicographical comparison isn't correct for it
   */
  bool operator == (const LinkInformation& other) const {
    return (
      indexPair == other.indexPair
      && cycleSequence == other.cycleSequence
    );
  }
};

std::vector<LinkInformation> substituentLinks(
  const GraphType& graph,
  const CycleData& cycleData,
  const AtomIndexType& source,
  const std::vector<AtomIndexType>& activeAdjacents
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

} // namespace molassembler

#endif
