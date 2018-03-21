#ifndef INCLUDE_GEOMETRY_MODELING_MODEL_H
#define INCLUDE_GEOMETRY_MODELING_MODEL_H

#include "chemical_symmetries/Symmetries.h"

#include "common_typedefs.h"
#include "boost/variant.hpp"

/*! @file
 *
 * Declarations for the general interface with which a number of classes can
 * determine the local geometry that a specific arrangement of atoms should
 * have.
 */

namespace molassembler {

namespace LocalGeometry {

/* Typedefs */

using ElementAndBondPair = std::pair<
  Delib::ElementType,
  BondType
>;

using LigandType = std::tuple<
  unsigned, // L
  unsigned, // X
  std::vector<
    ElementAndBondPair
  >
>;

// Mapping of bond type to a floating-point weight
extern const std::map<BondType, double> bondWeights;

/* Models */
boost::optional<Symmetry::Name> vsepr(
  const Delib::ElementType& centerAtomType,
  const unsigned& nSites,
  const std::vector<LigandType>& ligands,
  const int& formalCharge
);

boost::optional<Symmetry::Name> firstOfSize(const unsigned& size);

} // namespace LocalGeometry

} // namespace molassembler

#endif
