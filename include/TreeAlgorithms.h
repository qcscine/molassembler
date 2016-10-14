#ifndef INCLUDE_TREE_ALGORITHMS_H
#define INCLUDE_TREE_ALGORITHMS_H

#include "Tree.h"
#include "AdjacencyList.h"

namespace BasicTree {

using namespace MoleculeManip; 

using NodeType = BasicTree::Node<AtomIndexType>;

struct MakeTreeReturnType {
  std::shared_ptr<NodeType> rootPtr;
  std::vector<
    std::experimental::optional<
      std::shared_ptr<NodeType>
    >
  > nodes;
  std::vector<
    std::shared_ptr<NodeType> 
  > duplicateNodes;

  MakeTreeReturnType(const AdjacencyList& adjacencies) {
    nodes = std::vector<
      std::experimental::optional<
        std::shared_ptr<NodeType>
      >
    >(adjacencies.size(), std::experimental::nullopt);
  }
};

std::shared_ptr<NodeType> makeNodeRecursive(
  const AdjacencyList& adjacencies,
  MakeTreeReturnType& workStruct,
  const AtomIndexType& index,
  std::experimental::optional<
    std::shared_ptr<NodeType>
  > parentPtrOption 
) {
  auto node = std::make_shared<NodeType>(index);
  node -> parentOption = parentPtrOption;

  if(workStruct.nodes[index]) {
    // do not add any children
    workStruct.duplicateNodes.emplace_back(node);
  } else {
    workStruct.nodes[index] = node;
    for(const auto& adjacency : adjacencies.getAdjacencies(index)) {
      if(!( // exclude parent value from children nodes
          node -> parentOption
          && node -> parentOption.value() -> key == adjacency // was index before
      )) {
        node -> addChild(
          makeNodeRecursive( // recursive call
            adjacencies,
            workStruct,
            adjacency, 
            node  // a copy of the shared_ptr node is passed to optional ctor
          )
        );
      }
    }
  }

  return node;
}

MakeTreeReturnType makeTree(
  const AdjacencyList& adjacencies
) {
  MakeTreeReturnType workStruct(adjacencies); // passed to reveal size

  workStruct.rootPtr = makeNodeRecursive(
    adjacencies,
    workStruct,
    0,
    std::experimental::nullopt
  );

  return workStruct;
}


 /* const AtomIndexType& index,
  std::experimental::optional<
    std::shared_ptr<NodeType> 
  >& parentPtrOption
) {

}*/

}

#endif
