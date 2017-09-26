#include "ConstexprProperties.h"
#include "Properties.h"

namespace Symmetry {

namespace constexprProperties {

TemplateMagic::MinimalCache<
  std::pair<Symmetry::Name, Symmetry::Name>,
  ConstexprMagic::Optional<MappingsReturnType>
> mappingsCache;

const ConstexprMagic::Optional<MappingsReturnType>& getMapping(
  const Symmetry::Name& a,
  const Symmetry::Name& b
) {
  assert(a != b);

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
    return allMappings.at(
      static_cast<unsigned>(a),
      static_cast<unsigned>(b)
    );
  }
#endif

  /* Okay, so the mapping is definitely not in the set generated at
   * compile-time, so we have to fetch it from the cache or generate it
   */
  auto indexPair = std::make_pair(a, b);

  if(mappingsCache.has(indexPair)) {
    return mappingsCache.get(indexPair);
  }

  // NOTE: Compiling this expression below leads to the severe cost
  mappingsCache.add(
    indexPair,
    allMappingFunctions.at(
      static_cast<unsigned>(a),
      static_cast<unsigned>(b)
    )()
  );

  return mappingsCache.get(indexPair);
}

} // namespace constexprProperties

} // namespace Symmetry
