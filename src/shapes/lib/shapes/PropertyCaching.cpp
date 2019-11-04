/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "shapes/PropertyCaching.h"

#include "temple/constexpr/ToSTL.h"
#include "temple/constexpr/TupleTypePairs.h"
#include "temple/Functional.h"

namespace Scine {

namespace Shapes {

constexpr temple::Array<std::pair<double, double>, nShapes> symmetryAngleBounds = temple::TupleType::map<
  data::allShapeDataTypes,
  constexprProperties::AngleBoundsFunctor
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
template<typename SymmetrySource, typename SymmetryTarget>
struct mappingCalculationFunctor {
  static constexpr auto value() {
    return constexprProperties::calculateMapping<SymmetrySource, SymmetryTarget>();
  }
};

// Calculate the symmetryMapping for all possible combinations of symmetries
constexpr temple::UpperTriangularMatrix<
  temple::Optional<constexprProperties::MappingsReturnType>,
  nShapes * (nShapes - 1) / 2
> allMappings = temple::makeUpperTriangularMatrix(
  temple::TupleType::mapAllPairs<
    data::allShapeDataTypes,
    mappingCalculationFunctor
  >()
);
#endif

temple::MinimalCache<
  std::tuple<Shape, Shape, boost::optional<unsigned>>,
  properties::SymmetryTransitionGroup
> mappingsCache;

const boost::optional<const properties::SymmetryTransitionGroup&> getMapping(
  const Shape a,
  const Shape b,
  const boost::optional<unsigned>& removedIndexOption
) {
  if(a == b) {
    return boost::none;
  }

  auto cacheKey = std::make_tuple(a, b, removedIndexOption);

  if(mappingsCache.has(cacheKey)) {
    return mappingsCache.getOption(cacheKey);
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

      properties::SymmetryTransitionGroup STLResult;
      STLResult.indexMappings = temple::map(
        temple::toSTL(constexprMappings.mappings),
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
        cacheKey,
        STLResult
      );
    } else {
      // Calculate dynamically (relevant for targets of size 9 and higher)
      mappingsCache.add(
        cacheKey,
        properties::selectBestTransitionMappings(
          properties::symmetryTransitionMappings(a, b)
        )
      );
    }
#else
    mappingsCache.add(
      cacheKey,
      properties::selectBestTransitionMappings(
        properties::symmetryTransitionMappings(a, b)
      )
    );
#endif
  } else if(sizeDiff == -1 && removedIndexOption) {
    // Deletion case (always dynamic)
    mappingsCache.add(
      cacheKey,
      properties::selectBestTransitionMappings(
        properties::ligandLossTransitionMappings(a, b, removedIndexOption.value())
      )
    );
  }

  return mappingsCache.getOption(cacheKey);
}

#ifdef USE_CONSTEXPR_HAS_MULTIPLE_UNLINKED_STEREOPERMUTATIONS
template<typename Symmetry>
struct makeAllHasUnlinkedStereopermutationsFunctor {
  static constexpr auto value() {
    temple::DynamicArray<bool, constexprProperties::maxShapeSize> nums;

    /* Value for 0 is equal to value for 1, so calculate one less. When all
     * are equal, there is obviously only one stereopermutation, so there is no
     * need to calculate.
     */
    for(unsigned i = 0; i < Shapes::size - 1; ++i) {
      nums.push_back(
        constexprProperties::hasMultipleUnlinkedStereopermutations<Symmetry>(i + 1)
      );
    }

    return nums;
  }
};

constexpr temple::Array<
  temple::DynamicArray<bool, constexprProperties::maxShapeSize>,
  nShapes
> allHasMultipleUnlinkedStereopermutations = temple::TupleType::map<
  data::allShapeDataTypes,
  makeAllHasUnlinkedStereopermutationsFunctor
>();
#endif

temple::MinimalCache<
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

  auto stlMapped = temple::toSTL(dynArrRef);

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
      properties::hasMultipleUnlinkedStereopermutations(
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
