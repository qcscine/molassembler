/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/RankingInformation.h"

#include "temple/Adaptors/Zip.h"
#include "temple/Functional.h"
#include "temple/TinySet.h"

#include "molassembler/Containers/OrderDiscoveryHelper.h"
#include "molassembler/Cycles.h"

namespace Scine {

namespace molassembler {

LinkInformation::LinkInformation() = default;

LinkInformation::LinkInformation(
  std::pair<unsigned, unsigned> ligandIndices,
  std::vector<AtomIndex> sequence,
  const AtomIndex source
) {
  /* Fix degrees of freedom of the underlying information so we can
   * efficiently implement operator <. indexPair can be an ordered pair:
   */
  indexPair = std::move(ligandIndices);
  if(indexPair.first > indexPair.second) {
    std::swap(indexPair.first, indexPair.second);
  }

  // The cycle sequence should be centralized on the source vertex
  cycleSequence = centralizeRingIndexSequence(std::move(sequence), source);

  /* After centralization, the source vertex is first and last. We need to
   * fix the remaining degree of freedom, which is if that the cycle sequence
   * in between can be reversed. We choose to fix it by making it ascending
   * if there are at least four vertices in the sequence between the second
   * and second-to-last vertices
   */

  if(
    cycleSequence.size() >= 4
    && cycleSequence.at(1) > cycleSequence.at(cycleSequence.size() - 2)
  ) {
    // Reverse is [first, last), and cycleSequence is a vector, so:
    std::reverse(
      std::begin(cycleSequence) + 1,
      std::end(cycleSequence) - 1
    );
  }
}

void LinkInformation::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // A link's ligand indices do not change, but the sequence does
  for(AtomIndex& atomIndex : cycleSequence) {
    atomIndex = permutation.at(atomIndex);
  }

  // Conditional reverse of the sequence now depends on the new vertex indices
  if(
    cycleSequence.size() >= 4
    && cycleSequence.at(1) > cycleSequence.at(cycleSequence.size() -2)
  ) {
    std::reverse(
      std::begin(cycleSequence) + 1,
      std::end(cycleSequence) - 1
    );
  }
}

bool LinkInformation::operator == (const LinkInformation& other) const {
  return (
    indexPair == other.indexPair
    && cycleSequence == other.cycleSequence
  );
}

bool LinkInformation::operator != (const LinkInformation& other) const {
  return !(*this == other);
}

bool LinkInformation::operator < (const LinkInformation& other) const {
  return std::tie(indexPair, cycleSequence) < std::tie(other.indexPair, other.cycleSequence);
}


std::vector<unsigned> RankingInformation::ligandConstitutingAtomsRankedPositions(
  const std::vector<AtomIndex>& ligand,
  const RankingInformation::RankedType& sortedSubstituents
) {
  auto positionIndices = temple::map(
    ligand,
    [&sortedSubstituents](AtomIndex constitutingIndex) -> unsigned {
      return temple::find_if(
        sortedSubstituents,
        [&constitutingIndex](const auto& equalPrioritySet) -> bool {
          return temple::find(equalPrioritySet, constitutingIndex) != equalPrioritySet.end();
        }
      ) - sortedSubstituents.begin();
    }
  );

  std::sort(
    positionIndices.begin(),
    positionIndices.end(),
    std::greater<> {}
  );

  return positionIndices;
}

RankingInformation::RankedLigandsType RankingInformation::rankLigands(
  const LigandsType& ligands,
  const RankedType& sortedSubstituents
) {
  // Use partial ordering helper
  OrderDiscoveryHelper<unsigned> ligandOrdering;
  ligandOrdering.setUnorderedValues(
    temple::iota<unsigned>(
      ligands.size()
    )
  );

  // Pairwise compare all ligands, adding relational discoveries to the ordering helper
  for(unsigned i = 0; i < ligands.size(); ++i) {
    for(unsigned j = i + 1; j < ligands.size(); ++j) {
      const auto& a = ligands.at(i);
      const auto& b = ligands.at(j);

      if(a.size() < b.size()) {
        ligandOrdering.addLessThanRelationship(i, j);
      } else if(a.size() > b.size()) {
        ligandOrdering.addLessThanRelationship(j, i);
      } else {
        auto aPositions = ligandConstitutingAtomsRankedPositions(a, sortedSubstituents);
        auto bPositions = ligandConstitutingAtomsRankedPositions(b, sortedSubstituents);

        // lexicographical comparison
        if(aPositions < bPositions) {
          ligandOrdering.addLessThanRelationship(i, j);
        } else if(aPositions > bPositions) {
          ligandOrdering.addLessThanRelationship(j, i);
        }

        // No further differentiation.
      }
    }
  }

  // OrderingDiscoveryHelper's getSets() returns DESC order, so we reverse it to ASC
  auto finalSets = ligandOrdering.getSets();

  std::reverse(
    finalSets.begin(),
    finalSets.end()
  );

  return finalSets;
}

void RankingInformation::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // .sortedSubstituents is mapped by applying the vertex permutation
  for(auto& group : sortedSubstituents) {
    for(AtomIndex& atomIndex : group) {
      atomIndex = permutation.at(atomIndex);
    }
  }
  // .ligands too
  for(auto& group : ligands) {
    for(AtomIndex& atomIndex : group) {
      atomIndex = permutation.at(atomIndex);
    }
  }
  // .ligandsRanking is unchanged as it is ligand index based into .ligands
  // links does have to be mapped, though
  for(LinkInformation& link : links) {
    link.applyPermutation(permutation);
  }

  // Re-sort links to establish ordering
  std::sort(
    std::begin(links),
    std::end(links)
  );
}

unsigned RankingInformation::getLigandIndexOf(const AtomIndex i) const {
  // Find the atom index i within the set of ligand definitions
  auto findIter = temple::find_if(
    ligands,
    [&](const auto& ligandIndexList) -> bool {
      return temple::any_of(
        ligandIndexList,
        [&](const AtomIndex ligandConstitutingIndex) -> bool {
          return ligandConstitutingIndex == i;
        }
      );
    }
  );

  if(findIter == std::end(ligands)) {
    throw std::out_of_range("Specified atom index is not part of any ligand");
  }

  return findIter - std::begin(ligands);
}

unsigned RankingInformation::getRankedIndexOfLigand(const unsigned i) const {
  auto findIter = temple::find_if(
    ligandsRanking,
    [&](const auto& equallyRankedLigands) -> bool {
      return temple::any_of(
        equallyRankedLigands,
        [&](const unsigned ligandIndex) -> bool {
          return ligandIndex == i;
        }
      );
    }
  );

  if(findIter == std::end(ligandsRanking)) {
    throw std::out_of_range("Specified ligand index is not ranked.");
  }

  return findIter - std::begin(ligandsRanking);
}

bool RankingInformation::hasHapticLigands() const {
  return temple::any_of(
    ligands,
    [](const auto& ligandIndices) -> bool {
      return ligandIndices.size() > 1;
    }
  );
}

bool RankingInformation::operator == (const RankingInformation& other) const {
  /* This is a nontrivial operator since there is some degree of freedom in how
   * ligands are chosen
   */

  // Check all sizes
  if(
    sortedSubstituents.size() != other.sortedSubstituents.size()
    || ligands.size() != other.ligands.size()
    || ligandsRanking.size() != other.ligandsRanking.size()
    || links.size() != other.links.size()
  ) {
    return false;
  }

  // sortedSubstituents can be compare lexicographically
  if(sortedSubstituents != other.sortedSubstituents) {
    return false;
  }

  // Combined comparison of ligandsRanking with ligands
  if(
    !temple::all_of(
      temple::adaptors::zip(
        ligandsRanking,
        other.ligandsRanking
      ),
      [&](
        const auto& thisEqualLigandsGroup,
        const auto& otherEqualLigandsGroup
      ) -> bool {
        temple::TinyUnorderedSet<AtomIndex> thisLigandGroupVertices;

        for(const auto ligandIndex : thisEqualLigandsGroup) {
          for(const auto ligandConstitutingIndex : ligands.at(ligandIndex)) {
            thisLigandGroupVertices.insert(ligandConstitutingIndex);
          }
        }

        for(const auto ligandIndex : otherEqualLigandsGroup) {
          for(const auto ligandConstitutingIndex : other.ligands.at(ligandIndex)) {
            if(thisLigandGroupVertices.count(ligandConstitutingIndex) == 0) {
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


  // Compare the links (these are sorted since substituentLinks yields them like that)
  return links == other.links;
}

bool RankingInformation::operator != (const RankingInformation& other) const {
  return !(*this == other);
}

} // namespace molassembler

} // namespace Scine
