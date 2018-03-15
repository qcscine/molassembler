/* From GraphAlgorithms
 * - Not really useful since it will return only a single cycle the selected
 *   starting atom is in. Just the size information isn't really helpful during
 *   bounds collection, it would be much better if you could get all cycle 
 *   information en bloc -> URFs
 */

class SmallestContainingCycle : public boost::default_bfs_visitor {
private:
  // Side effect output
  boost::optional<
    std::set<AtomIndexType>
  >& _smallestCycle;

  // Const input parameter
  const unsigned _maxDepth;

  // State
  std::map<AtomIndexType, unsigned> _depthMap;
  std::vector<
    std::vector<AtomIndexType>
  > _sequences;

public:
  class EarlyExit {};

  SmallestContainingCycle(
    boost::optional<
      std::set<AtomIndexType>
    >& smallestCycle,
    const AtomIndexType& startingFrom,
    const unsigned& maxDepth,
    const unsigned& N
  ) : _smallestCycle(smallestCycle),
      _maxDepth(maxDepth),
      _depthMap({
        {startingFrom, 0}
      }),
      _sequences(N)
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
    // Maintain depth map
    _depthMap[boost::target(e, g)] = _depthMap[boost::source(e, g)] + 1;

    // Copy sequences from source to target and add source to it
    _sequences[boost::target(e, g)] = _sequences[boost::source(e, g)];
    _sequences[boost::target(e, g)].emplace_back(boost::source(e, g));
  }

  template<typename Graph>
  void non_tree_edge(const EdgeIndexType& e, const Graph& g) {
    _smallestCycle = std::set<AtomIndexType> {};
    for(const auto& index: _sequences[boost::source(e, g)]) {
      _smallestCycle->insert(index);
    }
    for(const auto& index: _sequences[boost::target(e, g)]) {
      _smallestCycle->insert(index);
    }

    throw EarlyExit();
  }
};

