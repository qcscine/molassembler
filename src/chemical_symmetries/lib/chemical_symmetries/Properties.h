#ifndef INCLUDE_SYMMETRY_INFORMATION_PROPERTIES_H
#define INCLUDE_SYMMETRY_INFORMATION_PROPERTIES_H

#include "ConstexprProperties.h"
#include "DynamicProperties.h"

/*! @file
 *
 * Contains the interface for property generation and access.
 */

namespace Symmetry {

/* Derived stored constexpr data */
//! The smallest angle between ligands in any symmetry
constexpr double smallestAngle __attribute__ ((unused))
= temple::TupleType::unpackToFunction<
  data::allSymmetryDataTypes,
  constexprProperties::minAngleFunctor
>();

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
//! All 0, +1 symmetry transition mappings
extern const temple::UpperTriangularMatrix<
  temple::Optional<constexprProperties::MappingsReturnType>,
  nSymmetries * (nSymmetries - 1) / 2
> allMappings;
#endif

/* Dynamic access to constexpr data */
//! Cache for on-the-fly generated mappings
extern temple::MinimalCache<
  std::tuple<Symmetry::Name, Symmetry::Name, boost::optional<unsigned>>,
  properties::SymmetryTransitionGroup
> mappingsCache;

//! Cached access to mappings. Populates the cache from constexpr if generated.
const boost::optional<const properties::SymmetryTransitionGroup&> getMapping(
  const Symmetry::Name a,
  const Symmetry::Name b,
  const boost::optional<unsigned>& removedIndexOption = boost::none
);

#ifdef USE_CONSTEXPR_HAS_MULTIPLE_UNLINKED_ASSIGNMENTS
extern const temple::Array<
  temple::DynamicArray<bool, constexprProperties::maxSymmetrySize>,
  nSymmetries
> allHasMultipleUnlinkedAssignments;
#endif

extern temple::MinimalCache<
  Symmetry::Name,
  std::vector<bool>
> hasMultipleUnlinkedCache;

bool hasMultipleUnlinkedAssignments(
  const Symmetry::Name symmetryName,
  unsigned nIdenticalLigands
);

} // namespace Symmetry

#endif
