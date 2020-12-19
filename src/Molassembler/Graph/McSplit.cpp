/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/McSplit.h"
#include "Molassembler/Molecule/AtomEnvironmentHash.h"
#include "Molassembler/Temple/Permutations.h"

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

  auto inverse = Temple::make_permutation(permutation).inverse();

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

bool check_sol(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  const std::vector<VtxPair>& solution
) {
  std::vector<bool> used_left(g0.n, false);
  std::vector<bool> used_right(g1.n, false);
  for (unsigned i=0; i<solution.size(); i++) {
    const VtxPair p0 = solution[i];
    if (used_left[p0.first] || used_right[p0.second]) {
      return false;
    }
    used_left[p0.first] = true;
    used_right[p0.second] = true;
    if (g0.label[p0.first] != g1.label[p0.second]) {
      return false;
    }
    for (unsigned j=i+1; j<solution.size(); j++) {
      VtxPair p1 = solution[j];
      if (g0.adjmat[p0.first][p1.first] != g1.adjmat[p0.second][p1.second]) {
        return false;
      }
    }
  }
  return true;
}

int select_bidomain(
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
      min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
      best = i;
    } else if (len == min_size) {
      int tie_breaker = find_min_value(left, bd.l, bd.left_len);
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
std::vector<Bidomain> filter_domains(
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

void solve(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  std::vector<VtxPair>& incumbent,
  std::vector<VtxPair>& current,
  std::vector<Bidomain>& domains,
  std::vector<int>& left,
  std::vector<int>& right,
  const unsigned matching_size_goal,
  const Arguments& arguments
) {
  if (current.size() > incumbent.size()) {
    incumbent = current;
  }

  unsigned bound = current.size() + calc_bound(domains);
  if (bound <= incumbent.size() || bound < matching_size_goal) {
    return;
  }

  if (arguments.big_first && incumbent.size()==matching_size_goal) {
    return;
  }

  int bd_idx = select_bidomain(domains, left, current.size(), arguments);
  // In the MCCS case, there may be nothing to branch on
  if (bd_idx == -1) {
    return;
  }
  Bidomain &bd = domains[bd_idx];

  int v = find_min_value(left, bd.l, bd.left_len);
  remove_vtx_from_left_domain(left, domains[bd_idx], v);

  // Try assigning v to each vertex w in the color class beginning at bd.r
  int w = -1;
  bd.right_len--;
  for (int i=0; i<=bd.right_len; i++) {
    const int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
    w = right[bd.r + idx];

    // swap w to the end of its colour class
    right[bd.r + idx] = right[bd.r + bd.right_len];
    right[bd.r + bd.right_len] = w;

    auto new_domains = filter_domains(
      domains, left, right, g0, g1, v, w,
      arguments.edge_labelled
    );
    current.emplace_back(v, w);
    solve(g0, g1, incumbent, current, new_domains, left, right, matching_size_goal, arguments);
    current.pop_back();
  }
  bd.right_len++;
  if (bd.left_len == 0) {
    remove_bidomain(domains, bd_idx);
  }
  solve(g0, g1, incumbent, current, domains, left, right, matching_size_goal, arguments);
}

std::vector<VtxPair> mcs(
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

  std::vector<VtxPair> incumbent;

  if (arguments.big_first) {
    for (int k=0; k < g0.n; k++) {
      const unsigned goal = g0.n - k;
      auto left_copy = left;
      auto right_copy = right;
      auto domains_copy = domains;
      std::vector<VtxPair> current;
      solve(g0, g1, incumbent, current, domains_copy, left_copy, right_copy, goal, arguments);
      if (incumbent.size() == goal/* || abort_due_to_timeout */) {
        break;
      }
    }
  } else {
    std::vector<VtxPair> current;
    solve(g0, g1, incumbent, current, domains, left, right, 1, arguments);
  }

  return incumbent;
}

std::vector<VtxPair> mcs(
  const PrivateGraph& g0,
  const PrivateGraph& g1,
  const bool connected,
  const bool labelEdges
) {
  Arguments arguments;
  arguments.edge_labelled = labelEdges;
  arguments.connected = connected;

  // Reorder the graphs by vertex degree and label edges
  LabeledGraph l0 {g0, labelEdges};
  LabeledGraph l1 {g1, labelEdges};

  auto mapping = mcs(l0, l1, arguments);
  assert(check_sol(l0, l1, mapping));

  // Resolve reordering back to original graph vertices
  std::vector<VtxPair> resolved;
  resolved.reserve(mapping.size());
  for(const VtxPair& p : mapping) {
    resolved.emplace_back(
      l0.permutation.at(p.first),
      l1.permutation.at(p.second)
    );
  }
  return resolved;
}

} // namespace McSplit
} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine
