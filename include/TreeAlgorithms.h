#ifndef INCLUDE_TREE_ALGORITHMS_H
#define INCLUDE_TREE_ALGORITHMS_H

#include "Tree.h"
#include "AdjacencyList.h"

namespace BasicTree {

using namespace MoleculeManip; 

using NodeType = BasicTree::Node<AtomIndexType>;

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
    while(adjacencies.getAdjacencies(index).size() != 0) {
      // pick the first adjacency
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

MakeTreeReturnType makeTree(
  const AdjacencyList& adjacencies
) {
  // initialize return struct
  MakeTreeReturnType workStruct; 
  workStruct.nodes.resize(adjacencies.size(), boost::none);

  // copy AdjacencyList
  AdjacencyList adjacencyCopy = adjacencies;

  // enter recursive call
  workStruct.rootPtr = makeNodeRecursive(
    adjacencyCopy,
    workStruct,
    0,
    boost::none
  );

  return workStruct;
}

}

#endif
