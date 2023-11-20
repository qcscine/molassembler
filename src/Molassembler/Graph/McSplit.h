/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_MCSPLIT_H
#define INCLUDE_MOLASSEMBLER_GRAPH_MCSPLIT_H

#include "Molassembler/Graph/PrivateGraph.h"
#include "boost/bimap.hpp"
#include <vector>

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {
namespace McSplit {

/* Implementation adapted from C++ implementation by James Trimble
 * from github.com/jamestrimble/ijcai2017-partitioning-common-subgraph
 * licensed under MIT License, reproduced below:
 *
 * Copyright (c) 2020 James Trimble
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Changes:
 * - Cleanup of unnecessary struct keywords
 * - Refactoring of global variables
 * - Addition of const all over the place to better indicate variables
 * - Alteration of return type and result tracking to yield all subgraphs of
 *   the maximum size
 */

enum class Heuristic { min_max, min_product };

struct Arguments {
  // Whether the maximum common subgraph must be a single connected component
  bool connected = false;
  // Whether the edges of the graphs are labelled
  bool edge_labelled = true;
  // Whether the vertices of the graphs are labelled
  bool vertex_labelled = true;
  // First try to find an induced subgraph isomorphism, then decrement the target size
  bool big_first = false;
  // Don't quite know what difference this makes yet
  Heuristic heuristic = Heuristic::min_max;
};

struct LabeledGraph {
  using AdjacencyMatrix = std::vector<std::vector<unsigned int>>;

  LabeledGraph(const PrivateGraph& g, bool labelEdges);

  inline void add_edge(unsigned v, unsigned w, unsigned edge_label) {
    adjmat[v][w] = 1 + edge_label;
    adjmat[w][v] = 1 + edge_label;
  }

  // Size
  int n;
  // Symmetric adjacency matrix with possibly edge labeled entries
  AdjacencyMatrix adjmat;
  // Vertex labeling
  std::vector<unsigned> label;
  // Vertices in labeled graph to source vertex indices permutation
  std::vector<unsigned> permutation;
};

using VtxPair = std::pair<int, int>;

struct Bidomain {
  Bidomain(int pl, int pr, int pleft_len, int pright_len, bool adjacent):
    l(pl),
    r(pr),
    left_len(pleft_len),
    right_len(pright_len),
    is_adjacent(adjacent) {};

  // Starting indices of left and right sets
  int l, r;
  int left_len, right_len;
  bool is_adjacent;
};

//! Verify a solution of the maximum common subgraph algorithm
bool checkSolution(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  const std::vector<VtxPair>& solution
);

int selectBidomain(
  const std::vector<Bidomain>& domains,
  const std::vector<int>& left,
  int current_matching_size,
  const Arguments& arguments
);

// multiway is for directed and/or labelled graphs
std::vector<Bidomain> filterDomains(
  const std::vector<Bidomain>& d,
  std::vector<int>& left,
  std::vector<int>& right,
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  int v,
  int w,
  bool multiway
);

struct AllVertexMappings {
  using Mapping = boost::bimap<AtomIndex, AtomIndex>;

  unsigned size = 0;
  std::vector<Mapping> mappings;

  static bool isHydrogenPermutation(
    const Mapping& a,
    const Mapping& b,
    const LabeledGraph& g1
  );
};

void solve(
  const PrivateGraph& g0,
  const PrivateGraph& g1,
  AllVertexMappings& incumbent,
  std::vector<VtxPair>& current,
  std::vector<Bidomain>& domains,
  std::vector<int>& left,
  std::vector<int>& right,
  unsigned matching_size_goal,
  const Arguments& arguments
);

AllVertexMappings mcs(
  const LabeledGraph& g0,
  const LabeledGraph& g1,
  const Arguments& arguments
);

/**
 * @brief Maximum common subgraph
 *
 * @param connected Whether the maximum common subgraph must be connected
 * @param labelEdges Whether to match bond types of edges. Element types are
 *   always matched on vertices.
 *
 * @return Vertex pairs in the common subgraph from @p g0 to @g1
 */
AllVertexMappings mcs(
  const PrivateGraph& g0,
  const PrivateGraph& g1,
  bool connected,
  bool labelEdges
);

} // namespace McSplit
} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine

#endif
