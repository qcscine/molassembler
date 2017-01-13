#include "AdjacencyListAlgorithms.h"

#include <iostream>

namespace MoleculeManip {

namespace AdjacencyListAlgorithms {

using NodeType = Tree::Node<AtomIndexType>;

std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies,
  const AtomIndexType& startingFrom
) {
  std::shared_ptr<NodeType> rootPtr = std::make_shared<NodeType>(startingFrom);

  std::map<
    AtomIndexType,
    std::vector<
      std::weak_ptr<NodeType>
    >
  > existingNodePtrMap {
    {startingFrom, {rootPtr}}
  };

  std::vector<AtomIndexType> onSameDepth;
  unsigned depthLevel = 0;
  
  auto indexVisitor = [&](const auto& atomIndex, const unsigned& depth) -> bool {
    // skip initial visit to root key
    if(atomIndex == startingFrom) return true;

    if(depth != depthLevel) {
      onSameDepth.clear();
      depthLevel = depth;
    }
    onSameDepth.push_back(atomIndex);

    std::vector<AtomIndexType> addAsChildToAll;

    for(const auto& adjacent: adjacencies.getAdjacencies(atomIndex)) {
      if(existingNodePtrMap.count(adjacent) == 1) {
        for(auto& existingWeakPtr : existingNodePtrMap.at(adjacent)) {
          auto newNode = std::make_shared<NodeType>(atomIndex);

          if(auto parentPtr = existingWeakPtr.lock()) {
            parentPtr -> addChild(newNode);
          }

          // add to existing NodePtrMap
          if(existingNodePtrMap.count(atomIndex) == 0) {
            existingNodePtrMap[atomIndex] = std::vector<
              std::weak_ptr<NodeType>
            >({newNode});
          } else {
            existingNodePtrMap[atomIndex].push_back(newNode);
          }
        }

        if(TemplateMagic::makeContainsPredicate(onSameDepth)(adjacent)) {
          addAsChildToAll.push_back(adjacent);
        }
      }
    }

    for(const auto& index : addAsChildToAll) {
      // except! nodes that have index as parent key
      for(auto& weakPtr : existingNodePtrMap.at(atomIndex)) {
        if(auto nodePtr = weakPtr.lock()) {
          if(auto parentPtr = nodePtr -> parentWeakPtr.lock()) {
            if(parentPtr -> key != index) {
              auto newNode = nodePtr -> addChild(index);
              existingNodePtrMap.at(index).push_back(newNode);
            } 
          } else {
            auto newNode = nodePtr -> addChild(index);
            existingNodePtrMap.at(index).push_back(newNode);
          }
        } 
      }
    }
    
    return true;
  };

  BFSVisit(
    adjacencies,
    startingFrom,
    indexVisitor,
    0
  );

  return rootPtr;
}

std::shared_ptr<NodeType> makeTree(
  const AdjacencyList& adjacencies
) {
  return makeTree(
    adjacencies,
    0
  );
}

std::vector<unsigned> _connectedComponents(const AdjacencyList& adjacencies) {
  if(adjacencies.size() == 0) return {};

  std::vector<unsigned> visited (
    adjacencies.size(),
    0
  );
  std::vector<unsigned>::iterator visIter = visited.begin();
  std::deque<AtomIndexType> toVisit = {0};
  unsigned counter = 1;
  AtomIndexType current;

  while(toVisit.size() > 0) {

    // take the first element
    current = toVisit[0];
    toVisit.pop_front();

    // mark it
    visited[current] = counter;

    /* add it's connections to toVisit if they are unvisited and not in the 
     *  deque already
     */
    for(const auto& connectedIndex : adjacencies[current]) {
      if(
        visited[connectedIndex] == 0 
        && std::find(
          toVisit.begin(),
          toVisit.end(),
          connectedIndex
        ) == toVisit.end()
      ) {
        toVisit.push_back(connectedIndex);
      }
    }

    // if toVisit is empty, increment the counter and add the next 
    // non-visited element to toVisit
    if(
      toVisit.size() == 0 &&
      visIter != visited.end()
    ) {
      // move visIter to the next unvisited element
      while(
        visIter != visited.end() 
        && *visIter != 0
      ) {
        visIter++;
      }
      // if we're not at the end, then
      if(visIter != visited.end()) {
        // increment the counter
        counter += 1;
        // and add it's index to toVisit
        toVisit.push_back(
          visIter - visited.begin()
        );
      }
    }
  }

  // all atoms must be visited
  assert(std::accumulate(
    visited.begin(),
    visited.end(),
    true,
    [](const bool& carry, const unsigned& element) {
      return (
        carry 
        && element != 0
      );
    }
  ));

  return visited;
}

unsigned numConnectedComponents(const AdjacencyList& adjacencies) {
  auto visited = _connectedComponents(adjacencies);
  if(visited.size() == 0) return 0;
  else return *max_element(visited.begin(), visited.end());
}

std::vector<
  std::vector<AtomIndexType>
> connectedComponentGroups(const AdjacencyList& adjacencies) {
  auto visited = _connectedComponents(adjacencies);

  // guard against 0-length
  if(visited.size() == 0) return std::vector<
    std::vector<AtomIndexType>
  >();

  unsigned num_groups = *max_element(visited.begin(), visited.end());
  std::vector<
    std::vector<AtomIndexType>
  > groups (num_groups);

  // go through visited
  // it could be e.g. [1, 2, 1, 3, 4, 3]
  for(auto it = visited.begin(); it != visited.end(); it++) {
    // make space if the group number is bigger
    if(groups.size() < *it) {
      groups.resize(*it);
    }

    groups[(*it) - 1].push_back( 
      it - visited.begin() 
    );
  }

  return groups;
} 

} // eo namespace AdjacencyListAlgorithms

} // eo namespace MoleculeManip
