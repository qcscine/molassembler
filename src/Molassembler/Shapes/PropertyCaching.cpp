/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Shapes/PropertyCaching.h"

#include "Molassembler/Temple/constexpr/ToStl.h"
#include "Molassembler/Temple/constexpr/TupleTypePairs.h"
#include "Molassembler/Temple/Functional.h"

namespace Scine {
namespace Molassembler {
namespace Shapes {

constexpr Temple::Array<std::pair<double, double>, nShapes> symmetryAngleBounds = Temple::Tuples::map<
  Data::allShapeDataTypes,
  ConstexprProperties::AngleBoundsFunctor
>();

double minimumAngle(const Shape shape) {
  return symmetryAngleBounds.at(
    static_cast<
      std::underlying_type_t<Shape>
    >(shape)
  ).first;
}

double maximumAngle(const Shape shape) {
  return symmetryAngleBounds.at(
    static_cast<
      std::underlying_type_t<Shape>
    >(shape)
  ).second;
}


boost::optional<const Properties::ShapeTransitionGroup&> getMapping(
  const Shape a,
  const Shape b,
  const boost::optional<Vertex>& removedIndexOption
) {
  using Key = std::tuple<Shape, Shape, boost::optional<unsigned>>;
  static Temple::MinimalCache<Key, Properties::ShapeTransitionGroup> mappingsCache;

  const Key key {a, b, removedIndexOption};

  if(mappingsCache.has(key)) {
    return mappingsCache.getOption(key);
  }

  int sizeDiff = static_cast<int>(Shapes::size(b)) - static_cast<int>(Shapes::size(a));

  if(sizeDiff == 1 || sizeDiff == 0) {
    mappingsCache.add(
      key,
      Properties::selectBestTransitionMappings(
        Properties::shapeTransitionMappings(a, b)
      )
    );
  } else if(sizeDiff == -1 && removedIndexOption) {
    // Deletion case (always dynamic)
    mappingsCache.add(
      key,
      Properties::selectBestTransitionMappings(
        Properties::ligandLossTransitionMappings(a, b, removedIndexOption.value())
      )
    );
  }

  return mappingsCache.getOption(key);
}


bool hasMultipleUnlinkedStereopermutations(
  const Shape shape,
  unsigned nIdenticalLigands
) {
  static Temple::MinimalCache<
    Shape,
    std::vector<bool>
  > hasMultipleUnlinkedCache;

  if(nIdenticalLigands == Shapes::size(shape)) {
    return false;
  }

  // Alias a call with 0 to a call with 1 since that is the first calculated value
  if(nIdenticalLigands == 0) {
    ++nIdenticalLigands;
  }

  if(hasMultipleUnlinkedCache.has(shape)) {
    return hasMultipleUnlinkedCache.get(shape).at(nIdenticalLigands - 1);
  }

  // Generate the cache element using dynamic properties
  std::vector<bool> unlinkedStereopermutations;
  for(unsigned i = 0; i < Shapes::size(shape) - 1; ++i) {
    unlinkedStereopermutations.push_back(
      Properties::hasMultipleUnlinkedStereopermutations(
        shape,
        i + 1
      )
    );
  }

  hasMultipleUnlinkedCache.add(
    shape,
    unlinkedStereopermutations
  );

  return unlinkedStereopermutations.at(nIdenticalLigands - 1);
}

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine
