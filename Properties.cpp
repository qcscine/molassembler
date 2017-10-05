#include "Properties.h"
#include "constexpr_magic/ToSTL.h"

namespace Symmetry {

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
/* Since the pointer-to-function of the instantiated function template is
 * identical for all symmetry data types (when using the optional-extended
 * version to avoid instantiating symmetryTransitionMappings with non-adjacent
 * symmetries), we can generate an upper triangular matrix of function pointers!
 */

template<typename SymmetrySource, typename SymmetryTarget>
struct mappingCalculationFunctor {
  static constexpr auto value() {
    // Return merely the function, not the evaluated result
    return constexprProperties::calculateMapping<SymmetrySource, SymmetryTarget>();
  }
};

/* Make function pointers to symmetryMapping for all possible combinations of
 * symmetry types
 */
constexpr auto precalculatedMappings = ConstexprMagic::makeUpperTriangularMatrix(
  ConstexprMagic::TupleType::mapAllPairs<
    data::allSymmetryDataTypes,
    mappingCalculationFunctor
  >()
);

const ConstexprMagic::UpperTriangularMatrix<
  ConstexprMagic::Optional<constexprProperties::MappingsReturnType>,
  nSymmetries * (nSymmetries - 1) / 2
> allMappings = precalculatedMappings;
#endif

TemplateMagic::MinimalCache<
  std::pair<Symmetry::Name, Symmetry::Name>,
  properties::SymmetryTransitionGroup
> mappingsCache;

const boost::optional<const properties::SymmetryTransitionGroup&> getMapping(
  const Symmetry::Name& a,
  const Symmetry::Name& b
) {
  assert(a != b);

  auto indexPair = std::make_pair(a, b);

  if(mappingsCache.has(indexPair)) {
    return mappingsCache.get(indexPair);
  }

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
  /* Is the desired mapping in the generated list of mappings?
   * It can be that allMappings contains only a limited set of mappings!
   *
   * WARNING: this assumes that the enum containing the names of symmetries has
   * the same order as the tuple specifying all (or some) symmetry data types
   * used at compile-time to generate allMappings
   */
  if(
    allMappings.N <= std::max(
      static_cast<unsigned>(a),
      static_cast<unsigned>(b)
    )
  ) {
    auto& constexprOption = allMappings.at(
      static_cast<unsigned>(a),
      static_cast<unsigned>(b)
    );

    if(constexprOption.hasValue()) {
      const auto& constexprMappings = constexprOption.value();

      properties::SymmetryTransitionGroup STLResult;
      STLResult.indexMappings = TemplateMagic::map(
        ConstexprMagic::toSTL(constexprMappings.mappings),
        [&](const auto& indexList) -> std::vector<unsigned> {
          return {
            indexList.begin(),
            indexList.end()
          };
        }
      );

      STLResult.angularDistortion = constexprMappings.angularDistortion;
      STLResult.chiralDistortion = constexprMappings.chiralDistortion;

      mappingsCache.add(
        indexPair,
        STLResult
      );
    }
  }
#else
  if(!mappingsCache.has(indexPair)
    && (std::set<int> {0, 1}).count(
      static_cast<int>(Symmetry::size(a))
      - static_cast<int>(Symmetry::size(b))
    ) == 1
  ) {
    mappingsCache.add(
      indexPair,
      properties::symmetryTransitionMappings(a, b)
    );
  }
#endif
  return mappingsCache.getOption(indexPair);
}

} // namespace Symmetry
