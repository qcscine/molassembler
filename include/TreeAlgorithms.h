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

  // Constructor with AdjacencyList
  MakeTreeReturnType(const AdjacencyList& adjacencies) {
    nodes = std::vector<
      boost::optional<
        std::shared_ptr<NodeType>
      >
    >(
      adjacencies.size(),
      boost::none
    ); 
  }
};

std::shared_ptr<NodeType> makeNodeRecursive(
  AdjacencyList& adjacencies, // this is modified in the recursive calls!
  MakeTreeReturnType& workStruct,
  const AtomIndexType& index,
  boost::optional<
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
    //for(const auto& adjacency : adjacencies.getAdjacencies(index)) {
    while(adjacencies.getAdjacencies(index).size() != 0) {
      auto adjacency = adjacencies.getAdjacencies(index)[0];
      if(!( // exclude parent value from children nodes
          node -> parentOption
          && node -> parentOption.value() -> key == adjacency // was index before
      )) {
          // remove adjacency
          adjacencies.removeAdjacency(index, adjacency);
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
  AdjacencyList adjacencyCopy = adjacencies;

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
