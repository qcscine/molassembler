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
= ConstexprMagic::TupleType::unpackToFunction<
  data::allSymmetryDataTypes,
  constexprProperties::minAngleFunctor
>();

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
//! All 0, +1 symmetry transition mappings, if calculated at compile-time
extern const ConstexprMagic::UpperTriangularMatrix<
  ConstexprMagic::Optional<constexprProperties::MappingsReturnType>,
  nSymmetries * (nSymmetries - 1) / 2
> allMappings;
#endif

/* Dynamic access to constexpr data */
//! Cache for on-the-fly generated mappings
extern TemplateMagic::MinimalCache<
  std::pair<Symmetry::Name, Symmetry::Name>,
  properties::SymmetryTransitionGroup
> mappingsCache;

//! Dynamic access to constexpr mappings
const boost::optional<const properties::SymmetryTransitionGroup&> getMapping(
  const Symmetry::Name& a,
  const Symmetry::Name& b
);

} // namespace Symmetry

#endif
