/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/McSplit.h"
#include "Molassembler/Molecule/AtomEnvironmentHash.h"
#include "Molassembler/Temple/Permutations.h"
#include "Molassembler/Temple/Functional.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <set>

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {
namespace McSplit {

LabeledGraph::LabeledGraph(
  const PrivateGraph& g,
  const bool labelEdges
) {
  n = g.V();
  adjmat = AdjacencyMatrix(n, std::vector<unsigned int>(n, 0U));
  label = std::vector<unsigned int>(n, 0U);

  std::vector<unsigned> degrees(n);
  for(const auto i : g.vertices()) {
    degrees[i] = g.degree(i);
  }
  permutation.resize(n);
  std::iota(std::begin(permutation), std::end(permutation), 0);
  std::stable_sort(
    std::begin(permutation),
    std::end(permutation),
    [&](int a, int b) { return degrees[a] > degrees[b]; }
  );

  auto inverse = Temple::Permutation::from(permutation).inverse();

  for(const PrivateGraph::Edge edge : g.edges()) {
    add_edge(
      inverse.at(g.source(edge)),
      inverse.at(g.target(edge)),
      labelEdges ? static_cast<unsigned>(g.bondType(edge)) : 0
    );
  }

  // Vertex labels are the element types only
  for(int i = 0; i < n; ++i) {
    label[i] = static_cast<unsigned>(g.elementType(permutation[i]));
  }
}

//! Calculate the upper bound on matchable vertices
inline int calculateBound(const std::vector<Bidomain>& domains) {
  int bound = 0;
  for (const Bidomain& bd : domains) {
    bound += std::min(bd.left_len, bd.right_len);
  }
  return bound;
}

//! Find the minimum value in part of a vector
inline int findMinValue(const std::vector<int>& arr, int start_idx, int len) {
  // NOTE: len is not always equal to the size of the array
  int min_v = std::numeric_limits<int>::max();
  for (int i=0; i<len; i++) {
    if (arr[start_idx + i] < min_v) {
      min_v = arr[start_idx + i];
    }
  }
  return min_v;
}

/* Returns the index of the smallest value in arr that is >w.
 *
 * Assumptions:
 * - such a value exists
 * - arr contains no duplicates
 * - arr has no values==INT_MAX
 */
inline int indexOfNextSmallest(
  const std::vector<int>& arr,
  const int start_idx,
  const int len,
  const int w
) {
  int idx = -1;
  int smallest = std::numeric_limits<int>::max();
  for (int i=0; i<len; i++) {
    if (arr[start_idx + i]>w && arr[start_idx + i]<smallest) {
      smallest = arr[start_idx + i];
      idx = i;
    }
  }
  return idx;
}

inline void removeVtxFromLeftDomain(
  std::vector<int>& left,
  Bidomain& bd,
  const int v
) {
  int i = 0;
  while(left[bd.l + i] != v) {
    i++;
  }
  std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
  bd.left_len--;
}

inline void removeBidomain(std::vector<Bidomain>& domains, const int idx) {
  domains[idx] = std::move(domains[domains.size()-1]);
  domains.pop_back();
}


// Returns length of left half of array
inline int partition(
  std::vector<int>& all_vv,
  const int start,
  const int len,
  const std::vector<unsigned>& adjrow
) {
  int i=0;
  for (int j=0; j < len; j++) {
    if (adjrow[all_vv[start+j]] > 0) {
      std::swap(all_vv[start+i], all_vv[start+j]);
      i++;
    }
  }
  return i;
}


bool checkSolution(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  const AllVertexMappings::Mapping& solution
) {
  std::vector<bool> used_left(g0.n, false);
  std::vector<bool> used_right(g1.n, false);
  const auto end = std::end(solution.left);
  for(auto i = std::begin(solution.left); i != end; ++i) {
    const unsigned v = i->first;
    const unsigned w = i->second;
    if(used_left[v] || used_right[w]) {
      return false;
    }
    used_left[v] = true;
    used_right[w] = true;
    if(g0.label[v] != g1.label[w]) {
      return false;
    }

    auto j = i;
    for(++j; j != end; ++j) {
      const unsigned x = j->first;
      const unsigned y = j->second;
      if(g0.adjmat[v][x] != g1.adjmat[w][y]) {
        return false;
      }
    }
  }

  return true;
}

int selectBidomain(
  const std::vector<Bidomain>& domains,
  const std::vector<int>& left,
  int current_matching_size,
  const Arguments& arguments
) {
  // Select the bidomain with the smallest max(leftsize, rightsize), breaking
  // ties on the smallest vertex index in the left set
  int min_size = std::numeric_limits<int>::max();
  int min_tie_breaker = std::numeric_limits<int>::max();
  int best = -1;
  const int numDomains = domains.size();
  for (int i=0; i<numDomains; i++) {
    const Bidomain& bd = domains[i];
    if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) {
      continue;
    }
    const int len = arguments.heuristic == Heuristic::min_max ?
      std::max(bd.left_len, bd.right_len) :
      bd.left_len * bd.right_len;
    if (len < min_size) {
      min_size = len;
      min_tie_breaker = findMinValue(left, bd.l, bd.left_len);
      best = i;
    } else if (len == min_size) {
      int tie_breaker = findMinValue(left, bd.l, bd.left_len);
      if (tie_breaker < min_tie_breaker) {
        min_tie_breaker = tie_breaker;
        best = i;
      }
    }
  }
  return best;
}

// Returns length of left half of array
// multiway is for directed and/or labelled graphs
std::vector<Bidomain> filterDomains(
  const std::vector<Bidomain>& d,
  std::vector<int>& left,
  std::vector<int>& right,
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  const int v,
  const int w,
  const bool multiway
) {
  std::vector<Bidomain> new_d;
  new_d.reserve(d.size());
  for (const Bidomain& old_bd : d) {
    int l = old_bd.l;
    int r = old_bd.r;
    /* After these two partitions, left_len and right_len are the lengths of the
     * arrays of vertices with edges from v or w (int the directed case, edges
     * either from or to v or w)
     */
    const int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
    const int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
    const int left_len_noedge = old_bd.left_len - left_len;
    const int right_len_noedge = old_bd.right_len - right_len;
    if ((left_len_noedge != 0) && (right_len_noedge != 0)) {
      new_d.emplace_back(l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent);
    }
    if (multiway && (left_len != 0) && (right_len != 0)) {
      const auto& adjrow_v = g0.adjmat[v];
      const auto& adjrow_w = g1.adjmat[w];
      auto l_begin = std::begin(left) + l;
      auto r_begin = std::begin(right) + r;
      std::sort(
        l_begin,
        l_begin + left_len,
        [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; }
      );
      std::sort(
        r_begin,
        r_begin + right_len,
        [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; }
      );
      const int l_top = l + left_len;
      const int r_top = r + right_len;
      while (l < l_top && r < r_top) {
        unsigned left_label = adjrow_v[left[l]];
        unsigned right_label = adjrow_w[right[r]];
        if (left_label < right_label) {
          l++;
        } else if (left_label > right_label) {
          r++;
        } else {
          const int lmin = l;
          const int rmin = r;
          do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
          do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
          new_d.emplace_back(lmin, rmin, l-lmin, r-rmin, true);
        }
      }
    } else if ((left_len != 0) && (right_len != 0)) {
      new_d.emplace_back(l, r, left_len, right_len, true);
    }
  }
  return new_d;
}

bool AllVertexMappings::isHydrogenPermutation(
  const Mapping& a,
  const Mapping& b,
  const LabeledGraph& g1
) {
  /* Left maps are ordered in their first index, so we can
   * sequentially compare elements of the mapping.
   *
   * The idea here is that if the mapping is sequence identical
   * regarding non-hydrogen elements, then it must be a permutation
   * thereof. This is much cheaper than using std::is_permutation.
   */
  return std::equal(
    std::begin(a.left),
    std::end(a.left),
    std::begin(b.left),
    std::end(b.left),
    [&](const auto& firstMap, const auto& secondMap) -> bool {
      constexpr auto hLabel = static_cast<unsigned>(Utils::ElementType::H);
      // Ensure left-sequence identical
      if(firstMap.first != secondMap.first) {
        return false;
      }

      /* Do not compare mapped indices if the elements of the mapped vertices
       * are hydrogen
       */
      if(
        g1.label[firstMap.second] == hLabel
        && g1.label[secondMap.second] == hLabel
      ) {
        return true;
      }

      // Compare target vertices of the mappings
      return (firstMap.second == secondMap.second);
    }
  );
}

void solve(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  AllVertexMappings& incumbent,
  AllVertexMappings::Mapping& current,
  std::vector<Bidomain>& domains,
  std::vector<int>& left,
  std::vector<int>& right,
  const unsigned matching_size_goal,
  const Arguments& arguments
) {
  if (current.size() > incumbent.size) {
    incumbent.mappings = {current};
    incumbent.size = current.size();
  } else if(current.size() == incumbent.size) {
    if(
      !Temple::any_of(
        incumbent.mappings,
        [&](const auto& mapping) -> bool {
          return AllVertexMappings::isHydrogenPermutation(mapping, current, g1);
        }
      )
    ) {
      incumbent.mappings.push_back(current);
    }
  }

  unsigned bound = current.size() + calculateBound(domains);
  if (bound < incumbent.size || bound < matching_size_goal) {
    return;
  }

  if (arguments.big_first && incumbent.size == matching_size_goal) {
    return;
  }

  int bd_idx = selectBidomain(domains, left, current.size(), arguments);
  // In the MCCS case, there may be nothing to branch on
  if (bd_idx == -1) {
    return;
  }
  Bidomain& bd = domains[bd_idx];

  int v = findMinValue(left, bd.l, bd.left_len);
  removeVtxFromLeftDomain(left, domains[bd_idx], v);

  // Try assigning v to each vertex w in the color class beginning at bd.r
  int w = -1;
  bd.right_len--;
  for (int i=0; i<=bd.right_len; i++) {
    const int idx = indexOfNextSmallest(right, bd.r, bd.right_len+1, w);
    w = right[bd.r + idx];

    // swap w to the end of its colour class
    right[bd.r + idx] = right[bd.r + bd.right_len];
    right[bd.r + bd.right_len] = w;

    auto new_domains = filterDomains(
      domains, left, right, g0, g1, v, w,
      arguments.edge_labelled
    );
    current.insert(AllVertexMappings::Mapping::value_type(v, w));
    solve(g0, g1, incumbent, current, new_domains, left, right, matching_size_goal, arguments);
    current.left.erase(v);
  }
  bd.right_len++;
  if (bd.left_len == 0) {
    removeBidomain(domains, bd_idx);
  }
  solve(g0, g1, incumbent, current, domains, left, right, matching_size_goal, arguments);
}

AllVertexMappings mcs(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  const Arguments& arguments
) {
  // Buffers of vertices for each partition
  std::vector<int> left;
  std::vector<int> right;

  std::vector<Bidomain> domains;

  std::set<unsigned> left_labels {std::begin(g0.label), std::end(g0.label)};
  std::set<unsigned> right_labels {std::begin(g1.label), std::end(g1.label)};
  std::vector<unsigned> common_labels;
  std::set_intersection(std::begin(left_labels),
      std::end(left_labels),
      std::begin(right_labels),
      std::end(right_labels),
      std::back_inserter(common_labels));

  // Create a bidomain for each label that appears in both graphs
  for (const unsigned label : common_labels) {
    const int start_l = left.size();
    const int start_r = right.size();

    for (int i=0; i<g0.n; i++) {
      if (g0.label[i]==label) {
        left.push_back(i);
      }
    }
    for (int i=0; i<g1.n; i++) {
      if (g1.label[i]==label) {
        right.push_back(i);
      }
    }

    const int left_len = left.size() - start_l;
    const int right_len = right.size() - start_r;
    domains.emplace_back(start_l, start_r, left_len, right_len, false);
  }

  AllVertexMappings incumbent;

  if (arguments.big_first) {
    // McSplit ↓
    for (int k=0; k < g0.n; k++) {
      const unsigned goal = g0.n - k;
      auto left_copy = left;
      auto right_copy = right;
      auto domains_copy = domains;
      AllVertexMappings::Mapping current;
      solve(g0, g1, incumbent, current, domains_copy, left_copy, right_copy, goal, arguments);
      if (incumbent.size == goal) {
        break;
      }
    }
  } else {
    // Regular McSplit
    AllVertexMappings::Mapping current;
    solve(g0, g1, incumbent, current, domains, left, right, 1, arguments);
  }

  return incumbent;
}

AllVertexMappings mcs(
  const PrivateGraph& g0,
  const PrivateGraph& g1,
  const bool connected,
  const bool labelEdges
) {
  Arguments arguments;
  arguments.edge_labelled = labelEdges;
  arguments.connected = connected;

  // Reorder the graphs by vertex degree and label edges
  const LabeledGraph l0 {g0, labelEdges};
  const LabeledGraph l1 {g1, labelEdges};

  const auto allMappings = mcs(l0, l1, arguments);
  assert(
    Temple::all_of(
      allMappings.mappings,
      [&](const auto& mapping) { return checkSolution(l0, l1, mapping); }
    )
  );

  // Resolve reordering back to original graph vertices
  AllVertexMappings resolved;
  resolved.size = allMappings.size;
  resolved.mappings = Temple::map(
    allMappings.mappings,
    [&](const auto& mapping) {
    AllVertexMappings::Mapping resolvedMapping;
      for(const auto& p : mapping.left) {
        resolvedMapping.insert(
          AllVertexMappings::Mapping::value_type(
            l0.permutation.at(p.first),
            l1.permutation.at(p.second)
          )
        );
      }
      return resolvedMapping;
    }
  );
  return resolved;
}

} // namespace McSplit
} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine
