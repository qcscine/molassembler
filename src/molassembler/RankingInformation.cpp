/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "molassembler/RankingInformation.h"

#include "molassembler/Temple/Adaptors/Zip.h"
#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/TinySet.h"
#include "molassembler/Temple/Poset.h"

#include "molassembler/Containers/OrderDiscoveryHelper.h"
#include "molassembler/Cycles.h"

namespace Scine {

namespace Molassembler {

RankingInformation::Link::Link() = default;

RankingInformation::Link::Link(
  std::pair<SiteIndex, SiteIndex> siteIndices,
  std::vector<AtomIndex> sequence,
  const AtomIndex source
) {
  /* Fix degrees of freedom of the underlying information so we can
   * efficiently implement operator <. sites can be an ordered pair:
   */
  sites = std::move(siteIndices);
  if(sites.first == sites.second) {
    throw std::runtime_error("Link site indices are identical!");
  }

  if(sites.first > sites.second) {
    std::swap(sites.first, sites.second);
  }

  // The cycle sequence should be centralized on the source vertex
  cycleSequence = centralizeRingIndexSequence(std::move(sequence), source);

  /* After centralization, the source vertex is first. We need to fix the
   * remaining degree of freedom, which is if that the cycle sequence in
   * between can be reversed. We choose to fix it by making it ascending if
   * there are at least three vertices in the sequence between the second and
   * second-to-last vertices
   */

  if(
    cycleSequence.size() >= 3
    && cycleSequence.at(1) > cycleSequence.back()
  ) {
    // Reverse is [first, last), and cycleSequence is a vector, so:
    std::reverse(
      std::begin(cycleSequence) + 1,
      std::end(cycleSequence)
    );
  }
}

void RankingInformation::Link::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // A link's site indices do not change, but the sequence does
  for(AtomIndex& atomIndex : cycleSequence) {
    atomIndex = permutation.at(atomIndex);
  }

  // Conditional reverse of the sequence now depends on the new vertex indices
  if(
    cycleSequence.size() >= 3
    && cycleSequence.at(1) > cycleSequence.back()
  ) {
    // Reverse is [first, last), and cycleSequence is a vector, so:
    std::reverse(
      std::begin(cycleSequence) + 1,
      std::end(cycleSequence)
    );
  }
}

bool RankingInformation::Link::operator == (const Link& other) const {
  return std::tie(sites, cycleSequence) == std::tie(other.sites, other.cycleSequence);
}

bool RankingInformation::Link::operator != (const Link& other) const {
  return !(*this == other);
}

bool RankingInformation::Link::operator < (const Link& other) const {
  return std::tie(sites, cycleSequence) < std::tie(other.sites, other.cycleSequence);
}


std::vector<unsigned> RankingInformation::siteConstitutingAtomsRankedPositions(
  const std::vector<AtomIndex>& siteAtomList,
  const RankingInformation::RankedSubstituentsType& substituentRanking
) {
  auto positionIndices = Temple::map(
    siteAtomList,
    [&substituentRanking](AtomIndex constitutingIndex) -> unsigned {
      return Temple::find_if(
        substituentRanking,
        [&constitutingIndex](const auto& equalPrioritySet) -> bool {
          return Temple::find(equalPrioritySet, constitutingIndex) != std::end(equalPrioritySet);
        }
      ) - std::begin(substituentRanking);
    }
  );

  std::sort(
    std::begin(positionIndices),
    std::end(positionIndices),
    std::greater<> {}
  );

  return positionIndices;
}

RankingInformation::RankedSitesType RankingInformation::rankSites(
  const SiteListType& sites,
  const RankedSubstituentsType& substituentRanking
) {
  const unsigned sitesSize = sites.size();

  Temple::Poset<SiteIndex> siteIndexPoset {Temple::iota<SiteIndex>(sitesSize)};

  // Order by site size
  siteIndexPoset.orderUnordered(
    [&](const SiteIndex i, const SiteIndex j) -> bool {
      return sites.at(i).size() < sites.at(j).size();
    }
  );

  // Then sub-order by ranked positions of their site constituting atoms
  siteIndexPoset.orderUnordered(
    [&](const SiteIndex i, const SiteIndex j) -> bool {
      return (
        siteConstitutingAtomsRankedPositions(sites.at(i), substituentRanking)
        < siteConstitutingAtomsRankedPositions(sites.at(j), substituentRanking)
      );
    }
  );

  return siteIndexPoset.extract();
}

void RankingInformation::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // .substituentRanking is mapped by applying the vertex permutation
  for(auto& group : substituentRanking) {
    for(AtomIndex& atomIndex : group) {
      atomIndex = permutation.at(atomIndex);
    }
  }
  // .sites too
  for(auto& group : sites) {
    for(AtomIndex& atomIndex : group) {
      atomIndex = permutation.at(atomIndex);
    }
  }
  // .siteRanking is unchanged as it is index based into .sites
  // .links do have to be mapped, though
  for(Link& link : links) {
    link.applyPermutation(permutation);
  }

  // Sort links to re-establish ordering
  std::sort(
    std::begin(links),
    std::end(links)
  );
}

SiteIndex RankingInformation::getSiteIndexOf(const AtomIndex i) const {
  // Find the atom index i within the set of ligand definitions
  auto findIter = Temple::find_if(
    sites,
    [&](const auto& siteAtomList) -> bool {
      return Temple::any_of(
        siteAtomList,
        [&](const AtomIndex siteConstitutingIndex) -> bool {
          return siteConstitutingIndex == i;
        }
      );
    }
  );

  if(findIter == std::end(sites)) {
    throw std::out_of_range("Specified atom index is not part of any ligand");
  }

  const unsigned idx = findIter - std::begin(sites);
  return SiteIndex(idx);
}

unsigned RankingInformation::getRankedIndexOfSite(const SiteIndex i) const {
  auto findIter = Temple::find_if(
    siteRanking,
    [&](const auto& equallyRankedSiteIndices) -> bool {
      return Temple::any_of(
        equallyRankedSiteIndices,
        [&](const unsigned siteIndex) -> bool {
          return siteIndex == i;
        }
      );
    }
  );

  if(findIter == std::end(siteRanking)) {
    throw std::out_of_range("Specified ligand index is not ranked.");
  }

  return findIter - std::begin(siteRanking);
}

bool RankingInformation::hasHapticSites() const {
  return Temple::any_of(
    sites,
    [](const auto& siteAtomList) -> bool {
      return siteAtomList.size() > 1;
    }
  );
}

bool RankingInformation::operator == (const RankingInformation& other) const {
  /* This is a nontrivial operator since there is some degree of freedom in how
   * sites are chosen
   */

  // Check all sizes
  if(
    substituentRanking.size() != other.substituentRanking.size()
    || sites.size() != other.sites.size()
    || siteRanking.size() != other.siteRanking.size()
    || links.size() != other.links.size()
  ) {
    return false;
  }

  // substituentRanking can be compared lexicographically
  if(substituentRanking != other.substituentRanking) {
    return false;
  }

  // Combined comparison of siteRanking with sites
  if(
    !Temple::all_of(
      Temple::Adaptors::zip(
        siteRanking,
        other.siteRanking
      ),
      [&](
        const auto& thisEqualSiteGroup,
        const auto& otherEqualSiteGroup
      ) -> bool {
        Temple::TinyUnorderedSet<AtomIndex> thisSiteGroupVertices;

        for(const auto siteIndex : thisEqualSiteGroup) {
          for(const auto siteConstitutingIndex : sites.at(siteIndex)) {
            thisSiteGroupVertices.insert(siteConstitutingIndex);
          }
        }

        for(const auto siteIndex : otherEqualSiteGroup) {
          for(const auto siteConstitutingIndex : other.sites.at(siteIndex)) {
            if(thisSiteGroupVertices.count(siteConstitutingIndex) == 0) {
              return false;
            }
          }
        }

        return true;
      }
    )
  ) {
    return false;
  }

  // Compare the links (these should be kept sorted throughout execution)
  assert(std::is_sorted(std::begin(links), std::end(links)));
  assert(std::is_sorted(std::begin(other.links), std::end(other.links)));
  return links == other.links;
}

bool RankingInformation::operator != (const RankingInformation& other) const {
  return !(*this == other);
}

} // namespace Molassembler

} // namespace Scine
