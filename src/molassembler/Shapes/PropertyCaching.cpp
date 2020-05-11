/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Shapes/PropertyCaching.h"

#include "molassembler/Temple/constexpr/ToStl.h"
#include "molassembler/Temple/constexpr/TupleTypePairs.h"
#include "molassembler/Temple/Functional.h"

namespace Scine {
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

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
template<typename ShapeSource, typename ShapeTarget>
struct mappingCalculationFunctor {
  static constexpr auto value() {
    return ConstexprProperties::calculateMapping<ShapeSource, ShapeTarget>();
  }
};

// Calculate the symmetryMapping for all possible combinations of symmetries
constexpr Temple::UpperTriangularMatrix<
  Temple::Optional<ConstexprProperties::MappingsReturnType>,
  nShapes * (nShapes - 1) / 2
> allMappings = Temple::makeUpperTriangularMatrix(
  Temple::Tuples::mapAllPairs<
    Data::allShapeDataTypes,
    mappingCalculationFunctor
  >()
);
#endif

Temple::MinimalCache<
  std::tuple<Shape, Shape, boost::optional<unsigned>>,
  Properties::ShapeTransitionGroup
> mappingsCache;

const boost::optional<const Properties::ShapeTransitionGroup&> getMapping(
  const Shape a,
  const Shape b,
  const boost::optional<Vertex>& removedIndexOption
) {
  if(a == b) {
    return boost::none;
  }

  using Key = typename decltype(mappingsCache)::key_type;

  const Key key {a, b, removedIndexOption};

  if(mappingsCache.has(key)) {
    return mappingsCache.getOption(key);
  }

  int sizeDiff = static_cast<int>(Shapes::size(b)) - static_cast<int>(Shapes::size(a));

  if(sizeDiff == 1 || sizeDiff == 0) {
#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
    /* Is the desired mapping in the generated list of mappings?
     * It can be that allMappings contains only a limited set of mappings!
     *
     * WARNING: this assumes that the enum containing the names of symmetries has
     * the same order as the tuple specifying all (or some) symmetry data types
     * used at compile-time to generate allMappings
     */
    auto& constexprOption = allMappings.at(
      std::min(
        static_cast<unsigned>(a),
        static_cast<unsigned>(b)
      ),
      std::max(
        static_cast<unsigned>(a),
        static_cast<unsigned>(b)
      )
    );

    if(constexprOption.hasValue()) {
      const auto& constexprMappings = constexprOption.value();

      Properties::ShapeTransitionGroup stlResult;
      stlResult.indexMappings = Temple::map(
        Temple::toSTL(constexprMappings.mappings),
        [&](const auto& indexList) -> std::vector<Vertex> {
          std::vector<Vertex> v;
          for(unsigned i : indexList) {
            v.emplace_back(i);
          }
          return v;
        }
      );

      stlResult.angularDistortion = constexprMappings.angularDistortion;
      stlResult.chiralDistortion = constexprMappings.chiralDistortion;

      mappingsCache.add(
        key,
        stlResult
      );
    } else {
      // Calculate dynamically (relevant for targets of size 9 and higher)
      mappingsCache.add(
        key,
        Properties::selectBestTransitionMappings(
          Properties::shapeTransitionMappings(a, b)
        )
      );
    }
#else
    mappingsCache.add(
      key,
      Properties::selectBestTransitionMappings(
        Properties::shapeTransitionMappings(a, b)
      )
    );
#endif
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

#ifdef USE_CONSTEXPR_HAS_MULTIPLE_UNLINKED_STEREOPERMUTATIONS
template<typename Symmetry>
struct makeAllHasUnlinkedStereopermutationsFunctor {
  static constexpr auto value() {
    Temple::DynamicArray<bool, ConstexprProperties::maxShapeSize> nums;

    /* Value for 0 is equal to value for 1, so calculate one less. When all
     * are equal, there is obviously only one stereopermutation, so there is no
     * need to calculate.
     */
    for(unsigned i = 0; i < Shapes::size - 1; ++i) {
      nums.push_back(
        ConstexprProperties::hasMultipleUnlinkedStereopermutations<Symmetry>(i + 1)
      );
    }

    return nums;
  }
};

constexpr Temple::Array<
  Temple::DynamicArray<bool, ConstexprProperties::maxShapeSize>,
  nShapes
> allHasMultipleUnlinkedStereopermutations = Temple::Tuples::map<
  Data::allShapeDataTypes,
  makeAllHasUnlinkedStereopermutationsFunctor
>();
#endif

Temple::MinimalCache<
  Shape,
  std::vector<bool>
> hasMultipleUnlinkedCache;

bool hasMultipleUnlinkedStereopermutations(
  const Shape shape,
  unsigned nIdenticalLigands
) {
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

#ifdef USE_CONSTEXPR_HAS_MULTIPLE_UNLINKED_STEREOPERMUTATIONS
  // Generate the cache element from constexpr non-STL data
  const auto& dynArrRef = allHasMultipleUnlinkedStereopermutations.at(
    static_cast<unsigned>(shape)
  );

  auto stlMapped = Temple::toSTL(dynArrRef);

  hasMultipleUnlinkedCache.add(
    shape,
    stlMapped
  );

  return stlMapped.at(nIdenticalLigands - 1);
#else
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
#endif
}

} // namespace Shapes
} // namespace Scine
