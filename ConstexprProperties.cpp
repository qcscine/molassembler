#include "ConstexprProperties.h"

namespace Symmetry {

namespace constexprProperties {

const ConstexprMagic::Optional<MappingsReturnType>& getMapping(
  const Symmetry::Name& a,
  const Symmetry::Name& b
) {
  assert(a != b);

  return allMappings.at(
    static_cast<unsigned>(a),
    static_cast<unsigned>(b)
  );
}

} // namespace constexprProperties

} // namespace Symmetry
