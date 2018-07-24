#include "RankingInformation.h"

#include "OrderDiscoveryHelper.h"
#include "temple/Containers.h"

namespace molassembler {

std::vector<unsigned> RankingInformation::ligandConstitutingAtomsRankedPositions(
  const std::vector<AtomIndexType>& ligand,
  const RankingInformation::RankedType& sortedSubstituents
) {
  auto positionIndices = temple::map(
    ligand,
    [&sortedSubstituents](AtomIndexType constitutingIndex) -> unsigned {
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

unsigned RankingInformation::getLigandIndexOf(const AtomIndexType i) const {
  // Find the atom index i within the set of ligand definitions
  auto findIter = temple::find_if(
    ligands,
    [&](const auto& ligandIndexList) -> bool {
      return temple::any_of(
        ligandIndexList,
        [&](const AtomIndexType ligandConstitutingIndex) -> bool {
          return ligandConstitutingIndex == i;
        }
      );
    }
  );

  if(findIter == ligands.end()) {
    throw std::out_of_range(
      "atom index is not part of a ligand of this stereocenter"
    );
  }

  return findIter - ligands.begin();
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

  throw std::logic_error("Not implemented!");

  return false;

  // TODO
}

bool RankingInformation::operator != (const RankingInformation& other) const {
  return !(*this == other);
}

} // namespace molassembler
