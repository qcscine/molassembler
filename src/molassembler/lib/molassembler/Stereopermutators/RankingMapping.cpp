#include "molassembler/Stereopermutators/RankingMapping.h"

#include "temple/Functional.h"

namespace Scine {
namespace molassembler {
namespace {

std::vector<AtomIndex> collectSubstituents(const RankingInformation::SiteListType& sites) {
  std::vector<AtomIndex> substituents;
  for(const auto& site : sites) {
    for(const AtomIndex index : site) {
      substituents.push_back(index);
    }
  }
  temple::inplace::sort(substituents);
  return substituents;
}

boost::optional<SiteMapping::Changes> determineChanges(
  const RankingInformation& a,
  const RankingInformation& b
) {
  auto aSubstituents = collectSubstituents(a.sites);
  auto bSubstituents = collectSubstituents(b.sites);

  const int substituentCountChange = bSubstituents.size() - aSubstituents.size();

  if(std::abs(substituentCountChange) > 1) {
    throw std::logic_error("SiteMappings can handle only a single substituent change at once");
  }

  if(substituentCountChange != 0) {
    std::vector<AtomIndex> changedSubstituents;

    std::set_symmetric_difference(
      std::begin(aSubstituents),
      std::end(aSubstituents),
      std::begin(bSubstituents),
      std::end(bSubstituents),
      std::back_inserter(changedSubstituents)
    );

    if(changedSubstituents.size() != 1) {
      throw std::logic_error("Multiple changed substituents");
    }

    if(substituentCountChange == +1) {
      // b has a substituent more
      return SiteMapping::Changes {
        changedSubstituents.front(),
        b.getSiteIndexOf(changedSubstituents.front())
      };
    }

    if(substituentCountChange == -1) {
      // a has a subtituent more
      return SiteMapping::Changes {
        changedSubstituents.front(),
        a.getSiteIndexOf(changedSubstituents.front())
      };
    }
  }

  return boost::none;
}

unsigned symmetricDifferenceSetSize(
  const std::vector<AtomIndex>& a,
  const std::vector<AtomIndex>& b
) {
  assert(std::is_sorted(std::begin(a), std::end(a)));
  assert(std::is_sorted(std::begin(b), std::end(b)));

  std::vector<AtomIndex> symmetricDifference;

  std::set_symmetric_difference(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::back_inserter(symmetricDifference)
  );

  return symmetricDifference.size();
}

} // namespace

SiteMapping SiteMapping::from(
  const RankingInformation& a,
  const RankingInformation& b
) {
  SiteMapping mapping;

  mapping.changes = determineChanges(a, b);

  const unsigned A = a.sites.size();
  const unsigned B = b.sites.size();

  std::vector<SiteIndex> mappedBs;

  // Map all sites that are identical
  for(unsigned i = 0; i < A; ++i) {
    const auto& atomSet = a.sites[i];
    auto findIter = temple::find(b.sites, atomSet);
    if(findIter != std::end(b.sites)) {
      const auto matchedB = SiteIndex(findIter - std::begin(b.sites));
      mapping.map.emplace(i, matchedB);
      mappedBs.push_back(matchedB);
    }
  }

  // If all sites are mapped by this, you're all done
  if(mapping.map.size() == std::min(A, B)) {
    return mapping;
  }

  /* Map all remaining sites using symmetric difference set size as indicator
   * of similarity
   */
  std::vector<SiteIndex> unmappedAs;
  for(SiteIndex i {0}; i < A; ++i) {
    if(mapping.map.count(i) == 0) {
      unmappedAs.push_back(i);
    }
  }
  const unsigned UA = unmappedAs.size();

  std::vector<SiteIndex> unmappedBs;
  for(SiteIndex i {0}; i < B; ++i) {
    if(temple::find(mappedBs, i) == std::end(mappedBs)) {
      unmappedBs.push_back(i);
    }
  }
  const unsigned UB = unmappedBs.size();

  auto calculateDistance = [&](const std::vector<unsigned>& permutation) -> unsigned {
    unsigned distanceSum = 0;
    for(unsigned i = 0; i < std::min(UA, UB); ++i) {
      distanceSum += symmetricDifferenceSetSize(
        a.sites.at(unmappedAs.at(i)),
        b.sites.at(unmappedBs.at(permutation.at(i)))
      );
    }
    return distanceSum;
  };

  auto permutation = temple::iota<unsigned>(std::max(UA, UB));
  auto bestPermutation = permutation;
  unsigned bestDistance = calculateDistance(permutation);

  while(std::next_permutation(std::begin(permutation), std::end(permutation))) {
    unsigned distance = calculateDistance(permutation);
    if(distance < bestDistance) {
      bestPermutation = permutation;
      bestDistance = distance;
    }
  }

  // Use the best permutation to finish the map
  for(unsigned i = 0; i < std::min(UA, UB); ++i) {
    mapping.map.emplace(
      unmappedAs.at(i),
      unmappedBs.at(bestPermutation.at(i))
    );
  }

  assert(mapping.map.size() == std::min(A, B));

  return mapping;
}

} // namespace molassembler
} // namespace Scine
