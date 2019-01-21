/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Graph/Canonicalization.h"

#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/OuterGraph.h"

#include "temple/Adaptors/SequentialPairs.h"
#include "temple/Functional.h"

extern "C" {
#include "nauty/traces.h"

void traces_canonicalize(int nv, size_t nde, size_t* v, int* d, int* e, size_t vlen, size_t dlen, size_t elen, int* lab, int* ptn) {
  DEFAULTOPTIONS_TRACES(options);
  // We want to specify a non-uniform vertex coloring
  options.defaultptn = false;

  TracesStats stats;

  DYNALLSTAT(int, orbits, orbits_sz);
  DYNALLOC1(int, orbits, orbits_sz, nv, "malloc");
  sparsegraph source {nde, v, nv, d, e, nullptr, vlen, dlen, elen, 0};

  Traces(&source, lab, ptn, orbits, &options, &stats, nullptr);

  /* After execution, vertex indices should be mapped pursuant to lab:
   *
   *   [i] 0 1 2 3 = new_idx
   *   lab 4 3 1 2 = old_idx
   */
  DYNFREE(orbits, orbits_sz);
}

} // end extern "C"

namespace Scine {
namespace molassembler {

std::vector<int> canonicalAutomorphism(
  const Molecule& mol,
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
   * Our strategy here is to allocate all the data we need to call Traces
   * in the C++ section and just pass pointers to it along to the C section so
   * as little as possible manual allocation happens there.
   */

  AtomIndex N = mol.graph().N();
  // std::size_t can exceed int
  if(N > std::numeric_limits<int>::max()) {
    throw std::domain_error("Graph size exceeds canonical labeling algorithm size limits");
  }

  if(hashes.size() != N) {
    throw std::invalid_argument("Supplied hashes do not match number of vertices");
  }

  int nv = mol.graph().N();
  std::size_t nde = 2 * mol.graph().B();
  std::vector<std::size_t> v(nv);
  std::vector<int> d(nv), e(nde);

  const InnerGraph& inner = mol.graph().inner();

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
  auto lab = temple::iota<int>(nv);
  std::sort(
    std::begin(lab),
    std::end(lab),
    [&hashes](const int a, const int b) -> bool {
      return hashes.at(a) < hashes.at(b);
    }
  );

  // And then generate the partition from checking adjacent hash equality
  auto ptn = temple::map(
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

  /* Call the C function with addresses to the start of the underlying arrays
   * and their sizes. The sizes of lab and ptn are known since they have to be
   * nv elements long.
   */
  traces_canonicalize(nv, nde, &v[0], &d[0], &e[0], v.size(), d.size(), e.size(), &lab[0], &ptn[0]);

  // lab now contains our result
  return lab;
}

} // namespace molassembler
} // namespace Scine
