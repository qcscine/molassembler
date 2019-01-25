/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Graph/Canonicalization.h"

#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/OuterGraph.h"

#include "temple/Adaptors/SequentialPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"

extern "C" {
#include "nauty/nausparse.h"

/*!
 * @brief Establish canonical labeling for a sparse graph representation
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

/*!
 * @brief Check whether the two groups of components of sparsegraphs could
 *   possibly be isomorphic and establish canonical labeling for both
 *
 * @post Establishes canonical labeling for each graph in @p flab and @p slab.
 *   The proposed isomorphism is, for each i in [0, N), flab[i] <-> slab[i].
 *
 * @returns The result of a basic graph comparison (not isomorphism) of
 *   canonical graphs without considering each graph's labeling or color cell
 *   values (which information is not present in a sparsegraph).
 *
 * @warning Heed the function name's warning and do not consider graphs that
 *   yield true here as isomorphic without comparing back to hash values.
 */
bool molassembler_nauty_possibly_isomorphic(
  int fnv,
  size_t fnde,
  size_t* fv,
  int* fd,
  int* fe,
  size_t fvlen,
  size_t fdlen,
  size_t felen,
  int* flab,
  int* fptn,

  int snv,
  size_t snde,
  size_t* sv,
  int* sd,
  int* se,
  size_t svlen,
  size_t sdlen,
  size_t selen,
  int* slab,
  int* sptn
) {
  DEFAULTOPTIONS_SPARSEGRAPH(options);
  // We want to specify non-uniform vertex colors
  options.defaultptn = false;
  // We want canonical graphs to be generated
  options.getcanon = true;
  // We want to use distance vertex invariants
  options.invarproc = distances_sg;

  statsblk stats;

  DYNALLSTAT(int, orbits, orbits_sz);
  DYNALLOC1(int, orbits, orbits_sz, fnv, "malloc");
  sparsegraph first_graph {fnde, fv, fnv, fd, fe, nullptr, fvlen, fdlen, felen, 0};
  sparsegraph second_graph {snde, sv, snv, sd, se, nullptr, svlen, sdlen, selen, 0};

  int m = SETWORDSNEEDED(fnv);
  nauty_check(WORDSIZE, m, fnv, NAUTYVERSIONID);

  SG_DECL(first_canon);
  SG_DECL(second_canon);

  sparsenauty(&first_graph, flab, fptn, orbits, &options, &stats, &first_canon);
  sparsenauty(&second_graph, slab, sptn, orbits, &options, &stats, &second_canon);

  /* lab, ptn is not taken into account here for direct graph comparison!
   * also, information of cell color values is lost already beforehand: methane
   * and silane are identically colored graphs!
   */
  bool possibly_isomorphic = aresame_sg(&first_canon, &second_canon);

  // Clean up
  SG_FREE(first_canon);
  SG_FREE(second_canon);
  DYNFREE(orbits, orbits_sz);

  return possibly_isomorphic;
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
    const InnerGraph& inner,
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
    if(N > std::numeric_limits<int>::max()) {
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
  const InnerGraph& inner,
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

bool isomorphic(
  const InnerGraph& firstGraph,
  const std::vector<hashes::WideHashType>& firstHashes,
  const InnerGraph& secondGraph,
  const std::vector<hashes::WideHashType>& secondHashes
) {
  assert(firstGraph.N() == secondGraph.N());
  detail::NautySparseGraph first(firstGraph, firstHashes), second(secondGraph, secondHashes);

  /* Generates a canonical labeling for both graphs and makes a very basic
   * direct graph comparison that does not include coloring or color cell
   * values. E.g. Silane and Methane are yielded as possibly isomorphic here.
   *
   * However, if this returns false, then we can early-exit, since if the direct
   * graph comparison post canonical labeling fails, then the graphs truly
   * cannot be isomorphic anymore, even when checking coloring.
   */
  if(
    !molassembler_nauty_possibly_isomorphic(
      first.nv,
      first.nde,
      &first.v[0],
      &first.d[0],
      &first.e[0],
      first.v.size(),
      first.d.size(),
      first.e.size(),
      &first.lab[0],
      &first.ptn[0],
      second.nv,
      second.nde,
      &second.v[0],
      &second.d[0],
      &second.e[0],
      second.v.size(),
      second.d.size(),
      second.e.size(),
      &second.lab[0],
      &second.ptn[0]
    )
  ) {
    return false;
  }

  /* Use the proposed isomorphism in first.lab and second.lab to compare vertex
   * hash values to figure out definitively whether the graphs are isomorphic.
   */
  return temple::all_of(
    temple::adaptors::range(firstGraph.N()),
    [&](const AtomIndex i) -> bool {
      return firstHashes.at(
        first.lab.at(i)
      ) == secondHashes.at(
        second.lab.at(i)
      );
    }
  );
}

} // namespace molassembler
} // namespace Scine
