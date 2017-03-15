/* These are remnants of the old tree generation algorithm
 * and the cycle detection algorithm. Use them for inspiration if needed when
 * cycle detection becomes relevant again.
 */

struct MakeTreeReturnType {
  // A pointer to the root node of the resulting tree
  std::shared_ptr<NodeType> rootPtr;

  // A vector (why optionals?) of pointers to the nodes
  std::vector<
    boost::optional<
      std::shared_ptr<NodeType>
    >
  > nodes;

  // A vector of pointers to duplicate nodes
  std::vector<
    std::shared_ptr<NodeType> 
  > duplicateNodes;
};

/* TODO 
 * - This algorithm is in dire need of a refactor. It's reaaaally hard to work
 *   with and has led to some abstruse APIs in Tree.h which make no sense
 */
std::shared_ptr<NodeType> makeNodeRecursive(
  AdjacencyList& adjacencies, 
  MakeTreeReturnType& workStruct,
  const AtomIndexType& index,
  boost::optional<
    std::shared_ptr<NodeType>
  > parentPtrOption 
) {
  auto nodePtr = std::make_shared<NodeType>(index);

  if(workStruct.nodes[index]) { // if this index has been registered before
    // do not add any children to it, stop here
    workStruct.duplicateNodes.emplace_back(nodePtr);
  } else {
    /* record that we have registered this node by copying it's pointer into
     * workStruct.nodes
     */
    workStruct.nodes[index] = nodePtr;

    // as long as there are adjacencies for this index in the AdjacencyList
    while(adjacencies.getNumAdjacencies(index) != 0) {
      // pick the first adjacency
      // TODO this is 
      auto adjacency = adjacencies.getAdjacencies(index).front();

      /* Since adjacencies are bi-directional (A has an entry for B and B has
       * an entry for A), we must exclude this node's parent from recursion
       */
      if(!( 
          parentPtrOption
          && (parentPtrOption.value()) -> key == adjacency 
      )) {
        // remove adjacency (this removes both entries)
        adjacencies.removeAdjacency(index, adjacency);

        // add children to the current node by recursion
        nodePtr -> addChild(
          makeNodeRecursive( 
            adjacencies, // same adjacencyList
            workStruct, // same struct
            adjacency, // the new index is the current adjacency index
            nodePtr // a copy of the shared_ptr is passed to optional ctor
          )
        );
      }
    }
  }

  return nodePtr;
}

/* detectCycles is deprecated along with makeTree
 * makeTree, with it's underlying recursive algorithm, does not create acyclic
 * trees like I want, so it's been superseded with a new algorithm, but I don't
 * want to delete this algorithm below outright, since it well need adaptation.
 * It's not currently in use since the aromatic ring detection algorithm has 
 * gone.
 */
/* Prerequisite: AdjacencyList has a single connected component. */
std::vector<
  std::vector<
    AtomIndexType
  >
> detectCycles(const AdjacencyList& adjacencies) {
  /* Rough outline of algorithm:
   *
   * Create a top-down tree from the AdjacencyList using a DFS-like algorithm.
   * Detect ring closure atoms (by detecting vertices that have been visited
   * before). The stored key in a node is the atom's index. Store node pointers
   * to repeated nodes.
   *
   * Then traverse the tree upwards from repeated nodes simultaneously until we hit 
   * a common ancestor. The common path is then the cycle.
   *
   * e.g. 
   *  
   *    0
   *   / \
   *  1 – 2
   *
   * is expanded into a tree, starting from 0, as:
   *
   * root/top ->   0 – 1 – 2 – 0   <- leaf/bottom
   *               ^           ^
   * repeated:     a           b
   *
   * then, as long as one of both node pointers has a parentOption and the sets 
   * of node keys of both pointers' traversal up the tree have an empty 
   * intersection, we move the pointers up the tree. If along the way we 
   * re-encounter the duplicate key, that pointer's traversal path is the cycle
   * set.
   *
   * step: 0 – 1 – 2 – 0    aSet = {}, bSet = {2}
   *       ^       ^
   *       a       b
   *
   * step: 0 – 1 – 2 – 0    aSet = {}, bSet = {2, 1}
   *       ^   ^
   *       a   b
   *
   * step: 0 – 1 – 2 – 0    aSet = {}, bSet = {2, 1, 0}
   *       ^^                                        ^
   *       ab                                        REPEAT FOUND
   *
   * -> Due to found repeating key, {2, 1, 0} is the cycle.
   */
  auto treeStruct = makeTree(adjacencies);

  std::vector<
    std::vector<
      AtomIndexType
    >
  > cycles;

  for(auto& candidateNodePtr : treeStruct.duplicateNodes) {
    AtomIndexType currentIndex = candidateNodePtr -> key;
    // get matching ptr from visited
    assert(treeStruct.nodes[candidateNodePtr -> key]);
    auto matchingPtr = treeStruct.nodes[candidateNodePtr -> key].value();

    // traverse the tree up to common ancestor

    std::pair<
      std::vector<AtomIndexType>,
      std::vector<AtomIndexType>
    > backtrackingPaths = { {candidateNodePtr -> key}, {matchingPtr -> key} };

    std::pair<
      std::set<AtomIndexType>,
      std::set<AtomIndexType>
    > backtrackingPathSets = {};

    std::vector<AtomIndexType> intersection;
    std::set_intersection(
      backtrackingPathSets.first.begin(),
      backtrackingPathSets.first.end(),
      backtrackingPathSets.second.begin(),
      backtrackingPathSets.second.end(),
      std::back_inserter(intersection)
    );

    boost::optional<
      std::vector<AtomIndexType>
    > optionFoundWhileBacktracking;

    //std::cout << "Constructing path for candidate: " << candidateNodePtr -> key << std::endl;

    while(
      (
        !(candidateNodePtr -> isRoot())
        || !(matchingPtr -> isRoot())
      ) && intersection.size() == 0
    ) {
      /*std::cout << "Step" << std::endl
        << candidateNodePtr << std::endl
        << matchingPtr << std::endl;*/
      /* Special case:
       * Entire ring is in single branch. In this case, if one of the
       * backtracking pointers encounters the current index as a key, then that
       * backtracking chain is a ring!
       */
      // move up the tree one step for both pointers if you can
      //if(candidateNodePtr -> parentOption) {
      if(candidateNodePtr = (candidateNodePtr -> parentWeakPtr).lock()) {
        if(candidateNodePtr -> key == currentIndex) {
          optionFoundWhileBacktracking = backtrackingPaths.first;
          break;
        }
        backtrackingPaths.first.push_back(candidateNodePtr -> key);
        backtrackingPathSets.first.insert(candidateNodePtr -> key);
      }
      if(matchingPtr = (matchingPtr -> parentWeakPtr).lock()) {
        if(matchingPtr -> key == currentIndex) {
          optionFoundWhileBacktracking = backtrackingPaths.second;
          break;
        }
        backtrackingPaths.second.push_back(matchingPtr -> key);
        backtrackingPathSets.second.insert(matchingPtr -> key);
      }

      // recompute the intersection
      intersection.clear();
      std::set_intersection(
        backtrackingPathSets.first.begin(),
        backtrackingPathSets.first.end(),
        backtrackingPathSets.second.begin(),
        backtrackingPathSets.second.end(),
        std::back_inserter(intersection)
      );
    }
    
    if(optionFoundWhileBacktracking) {
      cycles.push_back(optionFoundWhileBacktracking.value());
    } else {
      // the intersection is the "pivot" between both backtracking paths
      /*for(const auto& index: intersection) {
        std::cout << index << ", ";
      }
      std::cout << std::endl;*/
      assert(intersection.size() == 1);

      /* ring chain is then:
       * backtrackingPaths.first (forwards) up to but not including intersection[0]
       * + backtrackingPaths.second (backwards) from intersection[0] up to but not
       *    including backtrackingPaths.second[0] (this would be duplicate with
       *    backtrackingPaths.first[0], the common ring index)
       */

      std::vector<AtomIndexType> cycleChain;

      std::copy(
        backtrackingPaths.first.begin(),
        std::find( // position of the pivot
          backtrackingPaths.first.begin(),
          backtrackingPaths.first.end(),
          intersection[0]
        ),
        std::back_inserter(cycleChain)
      );

      std::copy(
        std::find( // position of pivot in reverse iteration
          backtrackingPaths.second.rbegin(),
          backtrackingPaths.second.rend(),
          intersection[0]
        ),
        backtrackingPaths.second.rend() - 1, // beginning + 1 (skip duplicate)
        std::back_inserter(cycleChain)
      );

      // add to cycles
      cycles.push_back(cycleChain);
    }
  }

  return cycles;
}
