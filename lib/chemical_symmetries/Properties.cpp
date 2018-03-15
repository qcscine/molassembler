#include "Properties.h"
#include "temple/constexpr/ToSTL.h"

namespace Symmetry {

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
template<typename SymmetrySource, typename SymmetryTarget>
struct mappingCalculationFunctor {
  static constexpr auto value() {
    return constexprProperties::calculateMapping<SymmetrySource, SymmetryTarget>();
  }
};

// Calculate the symmetryMapping for all possible combinations of symmetries
constexpr auto allMappings = temple::makeUpperTriangularMatrix(
  temple::TupleType::mapAllPairs<
    data::allSymmetryDataTypes,
    mappingCalculationFunctor
  >()
);
#endif

temple::MinimalCache<
  std::tuple<Symmetry::Name, Symmetry::Name, boost::optional<unsigned>>,
  properties::SymmetryTransitionGroup
> mappingsCache;

const boost::optional<const properties::SymmetryTransitionGroup&> getMapping(
  const Symmetry::Name& a,
  const Symmetry::Name& b,
  const boost::optional<unsigned>& removedIndexOption
) {
  if(a == b) {
    return boost::none;
  }

  auto cacheKey = std::make_tuple(a, b, removedIndexOption);

  if(mappingsCache.has(cacheKey)) {
    return mappingsCache.getOption(cacheKey);
  }

  int sizeDiff = static_cast<int>(Symmetry::size(b)) - static_cast<int>(Symmetry::size(a));

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
      STLResult.indexMappings = temple::mapToVector(
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

#ifdef USE_CONSTEXPR_NUM_UNLINKED_ASSIGNMENTS
template<typename Symmetry>
struct makeAllNumUnlinkedAssignmentsFunctor {
  static constexpr auto value() {
    temple::DynamicArray<unsigned, constexprProperties::maxSymmetrySize> nums;

    /* Value for 0 is equal to value for 1, so calculate one less. When all
     * are equal, there is obviously only one assignment, no need to calculate
     */
    for(unsigned i = 0; i < Symmetry::size - 1; ++i) {
      nums.push_back(
        constexprProperties::numUnlinkedAssignments<Symmetry>(i + 1)
      );
    }

    return nums;
  }
};

constexpr auto allNumUnlinkedAssignments = temple::TupleType::map<
  data::allSymmetryDataTypes,
  makeAllNumUnlinkedAssignmentsFunctor
>();
#endif

temple::MinimalCache<
  Symmetry::Name,
  std::vector<unsigned>
> numUnlinkedCache;

unsigned getNumUnlinked(
  const Symmetry::Name& symmetryName,
  unsigned nIdenticalLigands
) {
  if(nIdenticalLigands == Symmetry::size(symmetryName)) {
    return 1u;
  }

  if(nIdenticalLigands == 0) {
    ++nIdenticalLigands;
  }

  if(numUnlinkedCache.has(symmetryName)) {
    return numUnlinkedCache.get(symmetryName).at(nIdenticalLigands - 1);
  }

#ifdef USE_CONSTEXPR_NUM_UNLINKED_ASSIGNMENTS
  // Generate the cache element from constexpr non-STL data
  const auto& dynArrRef = allNumUnlinkedAssignments.at(
    static_cast<unsigned>(symmetryName)
  );

  auto stlMapped = temple::toSTL(dynArrRef);

  numUnlinkedCache.add(
    symmetryName,
    stlMapped
  );

  return stlMapped.at(nIdenticalLigands - 1);
#else
  // Generate the cache element using dynamic properties
  std::vector<unsigned> unlinkedAssignments;
  for(unsigned i = 0; i < Symmetry::size(symmetryName) - 1; ++i) {
    unlinkedAssignments.push_back(
      dynamicProperties::numUnlinkedAssignments(i + 1)
    );
  }

  numUnlinkedCache.add(
    symmetryName,
    unlinkedAssignments
  );

  return unlinkedAssignments.at(nIdenticalLigands - 1);
#endif
}

} // namespace Symmetry
