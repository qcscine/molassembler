#ifndef INCLUDE_MOLASSEMBLER_SUBGRAPHS_H
#define INCLUDE_MOLASSEMBLER_SUBGRAPHS_H

#include "boost/bimap.hpp"
#include "chemical_symmetries/Names.h"
#include "temple/constexpr/UpperTriangularMatrix.h"

#include "molassembler/Types.h"

#include <memory>

/*!@file
 *
 * Enables common subgraph computations
 */

namespace molassembler {

// Forward-declare Molecule
class Molecule;

namespace subgraphs {

//! Maybe consider using boost::bimap for this instead?
using IndexMap = boost::bimap<AtomIndex, AtomIndex>;

enum class VertexStrictness : unsigned {
  ElementType,
  LowEffortTransitionToLargerSymmetry,
  SameSymmetry
  // What else?
};

enum class EdgeStrictness : unsigned {
  BondType,
  EZIdentical
};

/*! Mappings for the maximum common subgraph between two molecules
 *
 * Finds an index mapping from a to b representing the maximum common subgraph
 * (if present).
 *
 * \warning For subgraph comparison, only element and bond types are considered
 * Stereocenters and Stereopermutations are not graph-local properties suitable
 * to a common substructure matching.
 */
std::vector<IndexMap> maximum(
  const Molecule& a,
  const Molecule& b,
  VertexStrictness vertexStrictness = VertexStrictness::ElementType,
  EdgeStrictness edgeStrictness = EdgeStrictness::BondType
);

namespace detail {

constexpr std::size_t upperTrigSize = Symmetry::nSymmetries * (Symmetry::nSymmetries - 1) / 2;

extern std::unique_ptr<
  temple::UpperTriangularMatrix<bool, upperTrigSize>
> lowEffortMappings;

void generateLowEffortMappings();

} // namespace detail

} // namespace subgraphs

} // namespace molassembler

#endif
