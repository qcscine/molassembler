/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 *
 * Canonicalization through nauty implementation.
 *
 * Note that it was attempted to use nauty to do graph isomorphisms by
 * canonicalizing both and then identity comparing, but the incurred graph
 * copies are just too expensive to compete with boost::isomorphism.
 */

#include "molassembler/Graph/Canonicalization.h"

#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Graph.h"

#include "temple/Adaptors/SequentialPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"

extern "C" {
#include "nauty/nausparse.h"

/*!
 * @brief Establish canonical labeling for a sparse graph representation
 *
 * @complexity{Generally sub-exponential in the number of vertices: @math{c^N}
 * with @math{c} some fixed constant.}
 *
 * @post lab contains an (inverse) index permutation yielding canonical graph
 *   labeling according to the coloring provided.
 */
void molassembler_nauty_canonicalize(int nv, size_t nde, size_t* v, int* d, int* e, size_t vlen, size_t dlen, size_t elen, int* lab, int* ptn) {
  DEFAULTOPTIONS_SPARSEGRAPH(options);
  // We want to specify a non-uniform vertex coloring
  options.defaultptn = false;
  // We want to use distance vertex invariants
  options.invarproc = distances_sg;
  /* Get a canonical graph (required for postcondition that lab contains
   * canonical labeling) despite us having no interest in the canonical graph
   * itself
   */
  options.getcanon = true;

  statsblk stats;

  SG_DECL(canong);

  DYNALLSTAT(int, orbits, orbits_sz);
  DYNALLOC1(int, orbits, orbits_sz, nv, "malloc");
  sparsegraph source {nde, v, nv, d, e, nullptr, vlen, dlen, elen, 0};

  int m = SETWORDSNEEDED(nv);
  nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

  sparsenauty(&source, lab, ptn, orbits, &options, &stats, &canong);

  /* After execution, vertex indices should be mapped pursuant to lab:
   *
   *   [i] 0 1 2 3 = new_idx
   *   lab 4 3 1 2 = old_idx
   */
  SG_FREE(canong);
  DYNFREE(orbits, orbits_sz);
}

} // end extern "C"

namespace Scine {
namespace molassembler {

namespace detail {

struct NautySparseGraph {
  // Number of vertices
  int nv;
  // Number of directed edges (all undirected edges count as 2)
  std::size_t nde;
  /* Adjacency-list structure:
   * - For each vertex i, d[i] is the out-degree of the vertex
   * - For each vertex i, v[i] is an index into the array e so that
   *   { e[v[i]], e[v[i] + 1], ..., e[v[i] + d[i] - 1] }
   *   are the vertices to which i is adjacent. This list does not need to be
   *   sorted.
   */
  std::vector<std::size_t> v;
  std::vector<int> d, e;
  /* Coloring and partioning:
   * - lab is a sorted list of vertex indices
   * - Each corresponding index in ptn indicates whether the current group of
   *   identically colored vertices (cell) is closed by its vertex in lab.
   */
  std::vector<int> lab, ptn;

  NautySparseGraph(
    const PrivateGraph& inner,
    const std::vector<hashes::WideHashType>& hashes
  ) {
    /* Sparse graph data structure:
     * - n is the number of vertices
     * - 'setword' = unsigned integer type of either 32 / 64 bits, depending on
     *   WORDSIZE compile-time parameter (link to nautyW / nautyL), this can be
     *   re-checked in C at runtime
     * - 'set' is a subset of vertices, that is represented by an array of m
     *   'setwords' so that WORDSIZE \cdot m >= n
     * - sparsegraph is a structure with the following fields
     *   - int nv: number of vertices
     *   - size_t nde: number of directed edges (all undirected edges count as 2)
     *   - size_t* v: pointer to array of length at least nv
     *   - int* d: pointer to array of length at least nv
     *   - int* e: pointer to array of length at least nde
     *   - sg_weight* w: not used, should be a nullptr
     *   - size_t vlen, dlen, elen, wlen: actual lengths of the vectors v, d, e
     *     and w. Since w is a nullptr, wlen should be zero
     *
     *   - For each vertex i, d[i] is the out-degree of the vertex.
     *   - v[i] is an index into the array e so that e[v[i]], e[v[i]+1], ...,
     *     e[v[i]+d[i]-1] are the vertices to which the vertex is joined. This
     *     list does not need to be sorted.
     *
     * Our strategy here is to allocate all the data we need to call nauty
     * in the C++ section and just pass pointers to it along to the C section so
     * as little as possible manual allocation happens there. Plus, we get to
     * keep the data afterwards.
     */

    const AtomIndex N = inner.N();
    // std::size_t can exceed int
    if(N > static_cast<AtomIndex>(std::numeric_limits<int>::max())) {
      throw std::domain_error("Graph size exceeds canonical labeling algorithm size limits");
    }

    if(hashes.size() != N) {
      throw std::invalid_argument("Supplied hashes do not match number of vertices");
    }

    nv = inner.N();
    nde = 2 * inner.B();
    v.reserve(nv);
    d.reserve(nv);
    e.reserve(nde);

    // Construct the adjacency list representation within a 'sparsegraph'
    for(AtomIndex i : boost::make_iterator_range(inner.vertices())) {
      d.push_back(inner.degree(i));

      // Create an adjacency list for i in v
      v.push_back(e.size());
      for(AtomIndex j : boost::make_iterator_range(inner.adjacents(i))) {
        e.push_back(j);
      }
    }

    /* A coloring is specified by a pair of integer arrays, called lab and ptn:
     * - lab contains a list of vertices in some order
     * - ptn indicates whether the group of identically colored vertices is
     *   over(0) or continues(1). The ordering of colors is encapsulated by the
     *   order of groups in lab!
     *
     * E.g.: lab: 4 1 3 2
     *       ptn: 0 1 0 0
     *       -> {4}, {1, 3}, {2}
     *
     * Our coloring is the hashing. Sort vertices by their hash and then make
     * groups according to the partitioning ptn.
     */

    // We use the hashes to order our vertices
    lab = temple::iota<int>(nv);
    std::sort(
      std::begin(lab),
      std::end(lab),
      [&hashes](const int a, const int b) -> bool {
        return hashes.at(a) < hashes.at(b);
      }
    );

    // And then generate the partition from checking adjacent hash equality
    ptn = temple::map(
      temple::adaptors::sequentialPairs(lab),
      [&hashes](const int a, const int b) -> int {
        return static_cast<int>(
          hashes.at(a) == hashes.at(b)
        );
      }
    );

    /* Add final element (since taking sequential pairs yields an element fewer),
     * which is always the end of a cell (group of vertices with identical
     * coloring)
     */
    ptn.push_back(0);
  }
};

} // namespace detail

std::vector<int> canonicalAutomorphism(
  const PrivateGraph& inner,
  const std::vector<hashes::WideHashType>& hashes
) {
  detail::NautySparseGraph nautyGraph(inner, hashes);

  /* Call the C function with addresses to the start of the underlying arrays
   * and their sizes. The sizes of lab and ptn are known since they have to be
   * nv elements long.
   */
  molassembler_nauty_canonicalize(
    nautyGraph.nv,
    nautyGraph.nde,
    &nautyGraph.v[0],
    &nautyGraph.d[0],
    &nautyGraph.e[0],
    nautyGraph.v.size(),
    nautyGraph.d.size(),
    nautyGraph.e.size(),
    &nautyGraph.lab[0],
    &nautyGraph.ptn[0]
  );

  // lab is now the (inverse) permutation we need to apply for a canonical graph
  return nautyGraph.lab;
}

} // namespace molassembler
} // namespace Scine
