#ifndef STDLIB_TYPE_ALGORITHMS_H
#define STDLIB_TYPE_ALGORITHMS_H

#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>

template<typename T>
std::ostream& operator << (std::ostream& os, const std::set<T>& rhs) {
  os << "set{";
  bool first = true;
  for(const auto& element: rhs) {
    if(first) {
      first = false;
      os << element;
    } else os << ", " << element;
  }
  os << "}";
  return os;
}

template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1, T2>& rhs) {
  os << "(" << rhs.first << ", " << rhs.second << ")";
  return os;
}

template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& rhs) {
  os << "vector{";
  bool first = true;
  for(const auto& element: rhs) {
    if(first) {
      first = false;
      os << element;
    } else os << ", " << element;
  }
  os << "}";
  return os;
}

namespace StdlibTypeAlgorithms {

template<typename T>
void mergeOverlappingSetsInplace(
  std::vector<
    std::set<T>
  >& sets
) {
  std::vector<T> intersection;
  for(unsigned i = 0; i < sets.size() - 1 && sets.size() > 1; i++) {
    for(unsigned j = i + 1; j < sets.size() && sets.size() > 1; j++) {
      intersection.clear();
      std::set_intersection(
        sets[i].begin(),
        sets[i].end(),
        sets[j].begin(),
        sets[j].end(),
        std::back_inserter(intersection)
      );
      if(intersection.size() > 0) {
        sets[i].insert(
          sets[j].begin(),
          sets[j].end()
        );
        sets.erase(sets.begin() + j);
        j = i; // redo entire row
      }
    }
  }
}

template<typename T>
std::vector<
  std::set<T>
> mergeOverlappingSets(
  const std::vector<
    std::set<T>
  >& immutableSets
) {
  auto setsCopy = immutableSets;
  mergeOverlappingSetsInplace(setsCopy);
  return setsCopy;
}


template<typename T>
std::vector<
  std::set<T>
> makeIndividualSets(
  const std::set<
    std::pair<
      T,
      T
    >
  >& pairsSet
) {

  std::vector<
    std::set<T>
  > sets;

  for(const auto& pair : pairsSet) {
    /* find a set in the vector containing one of the elements of this pair
     * and add the non-contained to the set
     */
    bool addedToExisting = false;
    std::vector<unsigned> countOverlap;
    for(auto& set : sets) {
      countOverlap.push_back(
        set.count(pair.first)
        + set.count(pair.second)
      );
    }

    unsigned numSetOverlaps = std::accumulate(
      countOverlap.begin(),
      countOverlap.end(),
      0u,
      [](const unsigned& carry, const unsigned& current) {
        if(current > 0) return carry + 1;
        else return carry;
      }
    );
    if(numSetOverlaps == 0) {
      // add a new set
      sets.emplace_back(std::set<T>({
        pair.first,
        pair.second
      }));
    } else {
      // add it to one set it has overlaps with
      for(unsigned i = 0; i < countOverlap.size(); i++) {
        if(countOverlap[i] > 0) {
          sets[i].insert(pair.first);
          sets[i].insert(pair.second);
          break;
        }
      }
    }

    if(!addedToExisting) {
      sets.emplace_back(std::set<T>({
        pair.first,
        pair.second
      }));
    }
  }

  // merge sets with overlap
  mergeOverlappingSetsInplace(sets);
  return sets;
}

template<typename T>
std::vector<T> copyMerge(
  const std::vector<T>& a,
  const std::vector<T>& b
) {
  std::vector<T> returnVector(a.size()+b.size());
  returnVector.insert(returnVector.end(), a.begin(), a.end());
  returnVector.insert(returnVector.end(), b.begin(), b.end());
  return returnVector;
}

template<typename T1, typename T2, typename ReturnType>
ReturnType minMaxAdaptor(
  const std::function<
    ReturnType(const T1&, const T2&)
  >& function,
  const T1& a,
  const T2& b
) {
  return function(
    std::min(a, b),
    std::max(a, b)
  );
}

template<typename T>
std::function<
  typename std::enable_if_t<std::is_function<T>::value, T>
> makeFunction(T *t) {
    return { t };
}

} // eo namespace

#endif
