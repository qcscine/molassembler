#ifndef INCLUDE_BASIC_TREE_H
#define INCLUDE_BASIC_TREE_H

#include <experimental/optional>
#include <memory>
#include <vector>
#include <iostream>
#include <deque>
#include <sstream>
#include <map>


namespace BasicTree {

template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::map<T1, T2>& map) {
  for(const auto& mapping : map) {
    os << mapping.first << " => " << mapping.second << std::endl;
  }

  return os;
}

template<typename T>
struct Node {
/* Public members */
  std::experimental::optional<
    std::shared_ptr<
      Node
    >
  > parentOption;
  std::vector<
    std::shared_ptr<
      Node
    >
  > children;

  T key;

/* Public member functions */
  /* Constructors */
  Node(const T& passKey) : key(passKey) {};

  void addChild(const std::shared_ptr<Node>& nodePtr) {
    children.push_back(nodePtr);
  }

  std::shared_ptr<Node>& addChild(const T& key) {
    children.emplace_back(key);
    return children.back();
  }

  /* Information */
  bool isRoot() const {
    return !parentOption;
  }

  bool isLeaf() const {
    return children.size() == 0;
  }
};

template<typename T>
std::ostream& operator << (
    std::ostream& os,
    const std::shared_ptr<
      Node<T> 
    >& rootPtr
) {
  std::vector<
    std::string
  > horizontal;

  std::map<T, unsigned> rowMap;

  auto printRows = [&os](const std::vector<std::string>& rowVector) {
    for(const auto& row : rowVector) {
      os << row << std::endl;
    }
  };

  std::deque<
    std::shared_ptr<
      Node<T>
    >
  > nodesToVisit = {rootPtr};

  unsigned initialSpace = 0;
  bool firstNode = true;
  bool firstRow = true;

  while(nodesToVisit.size() != 0) {
    //printRows(horizontal);

    auto current = nodesToVisit[0];
    nodesToVisit.pop_front();

    // add all children of current to nodesToVisit
    std::copy(
      current -> children.begin(),
      current -> children.end(),
      std::back_inserter(nodesToVisit)
    );

    if(firstNode) {
      unsigned rowNumber = 0;
      for(const auto& childrenPtr : current -> children) {
        std::stringstream row;
        if(firstRow) row << current -> key << "-" << childrenPtr -> key;
        else row << " `" << childrenPtr -> key;

        horizontal.push_back(row.str());
        rowMap[childrenPtr -> key] = rowNumber;

        rowNumber += 1;
        firstRow = false;
      }
      initialSpace += 3;

      firstNode = false;
    } else {
      firstRow = true;
      // get row this child key is in
      /*std::cout << "accessing index " << current -> key << " in map: " 
        << std::endl << rowMap << std::endl;*/
      auto row = rowMap.at(current -> key);
      auto rowNumber = row;
      for(const auto& childrenPtr : current -> children) {
        if(firstRow) {
          std::stringstream rowSS(
            horizontal[row],
            std::ios_base::in | std::ios_base::out | std::ios_base::ate
          );
          rowSS << "-" << childrenPtr -> key;
          rowMap[childrenPtr -> key] = rowNumber;
          horizontal[row] = rowSS.str();

          firstRow = false;
          rowNumber += 1;
        } else {
          std::stringstream rowSS;
          rowSS << std::string(initialSpace, ' ') << "`" << childrenPtr -> key;
          rowMap[childrenPtr -> key] = rowNumber;
          horizontal.insert(
            horizontal.begin() + rowNumber,
            rowSS.str()
          );
          // increase all unvisited map keys > rowNumber by one
          for(const auto& mapping : rowMap) {
            if(mapping.second >= rowNumber) {
              rowMap.at(mapping.first) += 1;
            }
          }
        }
      }
      initialSpace += 2;
    }
  }

  printRows(horizontal);

  return os;
}

} // eo namespace

#endif
