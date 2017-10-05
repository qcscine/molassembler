#ifndef INCLUDE_SYMMETRY_INFORMATION_PROPERTIES_H
#define INCLUDE_SYMMETRY_INFORMATION_PROPERTIES_H

#include "ConstexprProperties.h"
#include "DynamicProperties.h"

/*! @file
 *
 * Contains the interface for property generation and access.
 */

namespace Symmetry {

/* Since the pointer-to-function of the instantiated function template is
 * identical for all symmetry data types (when using the optional-extended
 * version to avoid instantiating symmetryTransitionMappings with non-adjacent
 * symmetries), we can generate an upper triangular matrix of function pointers!
 */

template<typename SymmetrySource, typename SymmetryTarget>
struct mappingCalculationFunctionPointerFunctor {
  static constexpr auto value() {
    // Return merely the function, not the evaluated result
    return constexprProperties::calculateMapping<SymmetrySource, SymmetryTarget>;
  }
};

/* Make function pointers to symmetryMapping for all possible combinations of
 * symmetry types
 */
constexpr auto allMappingFunctions __attribute__ ((unused)) 
= ConstexprMagic::makeUpperTriangularMatrix(
  ConstexprMagic::TupleType::mapAllPairs<
    data::allSymmetryDataTypes,
    mappingCalculationFunctionPointerFunctor
  >()
);

/* Derived stored constexpr data */
//! The smallest angle between ligands in any symmetry
constexpr double smallestAngle __attribute__ ((unused)) 
= ConstexprMagic::TupleType::unpackToFunction<
  data::allSymmetryDataTypes,
  constexprProperties::minAngleFunctor
>();

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
extern const ConstexprMagic::UpperTriangularMatrix<
  ConstexprMagic::Optional<constexprProperties::MappingsReturnType>,
  nSymmetries * (nSymmetries - 1) / 2
> allMappings;
#endif

/* Dynamic access to constexpr data */
//! Cache for on-the-fly generated mappings
extern TemplateMagic::MinimalCache<
  std::pair<Symmetry::Name, Symmetry::Name>,
  ConstexprMagic::Optional<constexprProperties::MappingsReturnType>
> mappingsCache;

//! Dynamic access to constexpr mappings
const ConstexprMagic::Optional<constexprProperties::MappingsReturnType>& getMapping(
  const Symmetry::Name& a,
  const Symmetry::Name& b
);

} // namespace Symmetry

#endif
