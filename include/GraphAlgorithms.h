#ifndef INCLUDE_GRAPH_ALGORITHMS_H
#define INCLUDE_GRAPH_ALGORITHMS_H

#include "common_typedefs.h"
#include "Tree.h"
#include <boost/graph/breadth_first_search.hpp>

#include <iostream>
#include "Delib/ElementInfo.h"

/* TODO
 * - Refactor to BGL
 *
 * NOTES
 * - Yes, it has to be a deque. A queue, although seemingly the minimal 
 *   functionality needed here for the FIFO structure, does not offer traversal
 *   without removal, which we need here.
 *
 */

namespace MoleculeManip {

namespace GraphAlgorithms {

//! BGL BFS Visitors
namespace BFSVisitors {

/* Basis on which algorithms with limited depth BFS iteration can be built. Is 
 * inexpensive for small depths in sparse graphs, but probably very space 
 * inefficient for large depths in dense graphs, where pre-filtering nodes is
 * (maybe) preferable.
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

class TreeGenerator : public boost::default_bfs_visitor {
public:
  using NodeType = Tree::Node<AtomIndexType>;  
  using KeyToNodeListMap = std::map<
    AtomIndexType,
    std::vector<
      std::weak_ptr<NodeType>
    >
  >;

private:
/* Internal state */
  // Map to keep track of depth from starting vertex
  std::map<AtomIndexType, unsigned> _depthMap;
  // Maximal depth as passed in constructor
  const unsigned _maxDepth;

  // Pointer to Tree root on heap (where we generate it)
  std::shared_ptr<NodeType> _rootPtr;
  // Mapping of atomic index to nodes in the Tree currently being generated
  KeyToNodeListMap _existingNodePtrMap;
  // List of edges to add in current depth
  std::vector<
    std::pair<AtomIndexType, AtomIndexType>
  > _edgesToAddBeforeDescent;


public:
  struct EarlyExit {};

  TreeGenerator(
    const AtomIndexType& startingFrom,
    const unsigned& maxDepth,
    std::shared_ptr<NodeType>& rootPtr
  ) : _depthMap({
        {startingFrom, 0}
      }),
      _maxDepth(maxDepth),
      _rootPtr(rootPtr),
      _existingNodePtrMap({
        {startingFrom, {_rootPtr}}
      })
  {}

  template<typename Graph>
  void discover_vertex(const AtomIndexType& u, const Graph& g) {
    assert(_depthMap.count(u));
    /*std::cout << "Discovered vertex " << Delib::ElementInfo::symbol(
      g[u].elementType
    ) << u << std::endl;*/
    if(_depthMap.at(u) > _maxDepth && _maxDepth != 0) {
      throw EarlyExit();
    }
  }

  template<typename Graph>
  void finish_vertex(const AtomIndexType& u, const Graph& g) {
    _edgesToAddBeforeDescent.erase(
      std::remove_if(
        _edgesToAddBeforeDescent.begin(),
        _edgesToAddBeforeDescent.end(),
        [&](const auto& edgePair) {
          bool canBeAdded = (
            _existingNodePtrMap.count(edgePair.first) > 0
            && _existingNodePtrMap.count(edgePair.second) > 0
          );

          // For every node with the target index
          if(canBeAdded) {
            for(const auto& parentWeakPtr : _existingNodePtrMap.at(edgePair.second)) {
              auto newNode = std::make_shared<NodeType>(edgePair.first);

              // Get a shared pointer from the weak one
              if(auto parentPtr = parentWeakPtr.lock()) {
                /* if the depth of the target node in the generated tree is
                 * less or equal to the depth of the source vertex in the graph
                 */
                if(parentPtr->depth() <= _depthMap[edgePair.first]) {
                  // Add the new child
                  parentPtr -> addChild(newNode);

                  // Add it to existingNodePtrMap
                  if(_existingNodePtrMap.count(edgePair.first) == 0) {
                    _existingNodePtrMap[edgePair.first] = std::vector<
                      std::weak_ptr<NodeType>
                    >({newNode});
                  } else {
                    _existingNodePtrMap[edgePair.first].push_back(newNode);
                  }
                }
              }
            }
          }

          return canBeAdded;
        }
      ),
      _edgesToAddBeforeDescent.end()
    );
  }

  template<typename Graph>
  void tree_edge(const EdgeIndexType& e, const Graph& g) {
    _depthMap[boost::target(e, g)] = _depthMap[boost::source(e, g)] + 1;
  }

  template<typename Graph>
  void non_tree_edge(const EdgeIndexType& e, const Graph& g) {
    auto source = boost::source(e, g), // where we are now
         target = boost::target(e, g); // maybe points to lower depth vertex

    /* Add the source as child of target, provided target is in the tree and
     * is at a lower depth than the source
     */
    if(
      _existingNodePtrMap.count(target)
      && _depthMap[target] < _depthMap[source]
    ) {
      // For every node with the target index
      for(const auto& parentWeakPtr : _existingNodePtrMap.at(target)) {
        auto newNode = std::make_shared<NodeType>(source);

        // Get a shared pointer from the weak one
        if(auto parentPtr = parentWeakPtr.lock()) {
          // Add the new child
          parentPtr -> addChild(newNode);
        }

        // Add it to existingNodePtrMap
        if(_existingNodePtrMap.count(source) == 0) {
          _existingNodePtrMap[source] = std::vector<
            std::weak_ptr<NodeType>
          >({newNode});
        } else {
          _existingNodePtrMap[source].push_back(newNode);
        }
      }
    } else if(_depthMap[target] == _depthMap[source]) { 
      // Need to take care of this special case for correct cycle decomposition
      // Don't know what exactly to do yet
      /*std::cout << source << " & " << target << " are on the same depth!" 
        << std::endl;
      std::cout << std::boolalpha << "Target exists already: " 
        << _existingNodePtrMap.count(target) << std::endl;*/
      _edgesToAddBeforeDescent.emplace_back(source, target);
    }
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

std::shared_ptr<
  BFSVisitors::TreeGenerator::NodeType
> makeTree(
  const GraphType& graph,
  const AtomIndexType& source,
  const unsigned& maxDepth
);

std::shared_ptr<
  BFSVisitors::TreeGenerator::NodeType
> makeTree(
  const GraphType& graph,
  const AtomIndexType& source
);

std::shared_ptr<
  BFSVisitors::TreeGenerator::NodeType
> makeTree(
  const GraphType& graph
);

std::set<
  std::pair<AtomIndexType, AtomIndexType>
> findSubstituentLinks(
  const GraphType& graph,
  const AtomIndexType& source,
  const std::set<AtomIndexType>& activeSubstituents
);

unsigned numConnectedComponents(const GraphType& graph);

} // namespace GraphAlgorithms

} // namespace MoleculeManip

#endif
