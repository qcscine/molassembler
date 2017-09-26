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
 *
 * Is proven by the following static_assert (commented out for performance)
 */
/*static_assert(
  std::is_same<
    decltype(calculateMapping<data::Linear, data::Bent>),
    decltype(calculateMapping<data::Linear, data::TrigonalPlanar>)
  >::value,
  "pointer-to-function types are not identical"
);*/

template<typename SymmetrySource, typename SymmetryTarget>
struct mappingCalculationFunctionPointerFunctor {
  static constexpr auto value() {
    // Return merely the function, not the evaluated result
    return constexprProperties::calculateMapping<SymmetrySource, SymmetryTarget>;
  }
};

constexpr auto allMappingFunctions
= ConstexprMagic::makeUpperTriangularMatrix(
  ConstexprMagic::TupleType::mapAllPairs<
    data::limitedSymmetryDataTypes,
    mappingCalculationFunctionPointerFunctor
  >()
);

template<typename SymmetrySource, typename SymmetryTarget>
struct mappingCalculationFunctor {
  static constexpr ConstexprMagic::Optional<constexprProperties::MappingsReturnType> value() {
    // Return the evaluated result
    return allMappingFunctions.at(
      static_cast<unsigned>(SymmetrySource::name),
      static_cast<unsigned>(SymmetryTarget::name)
    )();
  }
};

/* Derived stored constexpr data */
constexpr double smallestAngle __attribute__ ((unused)) 
= ConstexprMagic::TupleType::unpackToFunction<
  data::allSymmetryDataTypes,
  constexprProperties::minAngleFunctor
>();

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
constexpr auto allMappings = ConstexprMagic::makeUpperTriangularMatrix(
  ConstexprMagic::TupleType::mapAllPairs<
    data::limitedSymmetryDataTypes,
    mappingCalculationFunctor
  >()
);
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
