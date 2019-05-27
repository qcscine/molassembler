/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "chemical_symmetries/Properties.h"

#include "temple/constexpr/ToSTL.h"
#include "temple/constexpr/TupleTypePairs.h"
#include "temple/Functional.h"

namespace Scine {

namespace Symmetry {

constexpr temple::Array<std::pair<double, double>, nSymmetries> symmetryAngleBounds = temple::TupleType::map<
  data::allSymmetryDataTypes,
  constexprProperties::AngleBoundsFunctor
>();

double minimumAngle(const Symmetry::Name symmetryName) {
  return symmetryAngleBounds.at(
    static_cast<
      std::underlying_type_t<Symmetry::Name>
    >(symmetryName)
  ).first;
}

double maximumAngle(const Symmetry::Name symmetryName) {
  return symmetryAngleBounds.at(
    static_cast<
      std::underlying_type_t<Symmetry::Name>
    >(symmetryName)
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
  nSymmetries * (nSymmetries - 1) / 2
> allMappings = temple::makeUpperTriangularMatrix(
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
  const Symmetry::Name a,
  const Symmetry::Name b,
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
    temple::DynamicArray<bool, constexprProperties::maxSymmetrySize> nums;

    /* Value for 0 is equal to value for 1, so calculate one less. When all
     * are equal, there is obviously only one stereopermutation, so there is no
     * need to calculate.
     */
    for(unsigned i = 0; i < Symmetry::size - 1; ++i) {
      nums.push_back(
        constexprProperties::hasMultipleUnlinkedStereopermutations<Symmetry>(i + 1)
      );
    }

    return nums;
  }
};

constexpr temple::Array<
  temple::DynamicArray<bool, constexprProperties::maxSymmetrySize>,
  nSymmetries
> allHasMultipleUnlinkedStereopermutations = temple::TupleType::map<
  data::allSymmetryDataTypes,
  makeAllHasUnlinkedStereopermutationsFunctor
>();
#endif

temple::MinimalCache<
  Symmetry::Name,
  std::vector<bool>
> hasMultipleUnlinkedCache;

bool hasMultipleUnlinkedStereopermutations(
  const Symmetry::Name symmetryName,
  unsigned nIdenticalLigands
) {
  if(nIdenticalLigands == Symmetry::size(symmetryName)) {
    return false;
  }

  // Alias a call with 0 to a call with 1 since that is the first calculated value
  if(nIdenticalLigands == 0) {
    ++nIdenticalLigands;
  }

  if(hasMultipleUnlinkedCache.has(symmetryName)) {
    return hasMultipleUnlinkedCache.get(symmetryName).at(nIdenticalLigands - 1);
  }

#ifdef USE_CONSTEXPR_HAS_MULTIPLE_UNLINKED_STEREOPERMUTATIONS
  // Generate the cache element from constexpr non-STL data
  const auto& dynArrRef = allHasMultipleUnlinkedStereopermutations.at(
    static_cast<unsigned>(symmetryName)
  );

  auto stlMapped = temple::toSTL(dynArrRef);

  hasMultipleUnlinkedCache.add(
    symmetryName,
    stlMapped
  );

  return stlMapped.at(nIdenticalLigands - 1);
#else
  // Generate the cache element using dynamic properties
  std::vector<bool> unlinkedStereopermutations;
  for(unsigned i = 0; i < Symmetry::size(symmetryName) - 1; ++i) {
    unlinkedStereopermutations.push_back(
      properties::hasMultipleUnlinkedStereopermutations(
        symmetryName,
        i + 1
      )
    );
  }

  hasMultipleUnlinkedCache.add(
    symmetryName,
    unlinkedStereopermutations
  );

  return unlinkedStereopermutations.at(nIdenticalLigands - 1);
#endif
}

} // namespace Symmetry

} // namespace Scine
