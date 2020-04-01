/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/RankingMapping.h"

#include "molassembler/Temple/Functional.h"

namespace Scine {
namespace molassembler {
namespace {

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

boost::optional<SiteIndex> determineChangedSite(
  const unsigned A,
  const unsigned B,
  const SiteMapping::Map& map,
  const std::vector<SiteIndex> mappedBs
) {
  if(A == B) {
    return boost::none;
  }

  if(A < B) {
    for(SiteIndex i {0}; i < B; ++i) {
      if(temple::find(mappedBs, i) == std::end(mappedBs)) {
        return i;
      }
    }
  }

  if(A > B) {
    for(SiteIndex i {0}; i < A; ++i) {
      if(map.count(i) == 0) {
        return i;
      }
    }
  }

  throw std::logic_error("Failed to deduce changed site");
}

boost::optional<SiteIndex> determineChangedSite(
  const unsigned A,
  const unsigned B,
  const SiteMapping::Map& map,
  const std::vector<SiteIndex> possiblyUnmappedAs,
  const std::vector<SiteIndex> mappedBs
) {
  if(A == B) {
    return boost::none;
  }

  if(A < B) {
    for(SiteIndex i {0}; i < B; ++i) {
      if(temple::find(mappedBs, i) == std::end(mappedBs)) {
        return i;
      }
    }
  }

  if(A > B) {
    for(const SiteIndex unmappedA : possiblyUnmappedAs) {
      if(map.count(unmappedA) == 0) {
        return unmappedA;
      }
    }
  }

  throw std::logic_error("Failed to deduce changed site");
}

} // namespace

SiteMapping SiteMapping::from(
  const RankingInformation& a,
  const RankingInformation& b
) {
  SiteMapping mapping;

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

  // If all sites are mapped by this, you're almost done
  if(mapping.map.size() == std::min(A, B)) {
    mapping.changedSite = determineChangedSite(A, B, mapping.map, mappedBs);
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

  auto permutation = temple::iota<unsigned>(UB);
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
    const SiteIndex mappedB = unmappedBs.at(bestPermutation.at(i));
    mapping.map.emplace(unmappedAs.at(i), mappedB);
    mappedBs.push_back(mappedB);
  }

  assert(mapping.map.size() == std::min(A, B));

  mapping.changedSite = determineChangedSite(A, B, mapping.map, unmappedAs, mappedBs);
  return mapping;
}

} // namespace molassembler
} // namespace Scine
