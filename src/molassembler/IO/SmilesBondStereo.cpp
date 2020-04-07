/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/IO/SmilesBondStereo.h"

#include "molassembler/Molecule.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Temple/Optionals.h"

namespace Scine {
namespace molassembler {
namespace io {

unsigned SmilesBondStereo::findAssignment(
  BondStereopermutator stereopermutator,
  const Molecule& mol,
  const std::vector<PrivateGraph::Vertex>& indexInComponentMap
) const {
  auto first = mol.stereopermutators().option(stereopermutator.placement().first).value();
  auto second = mol.stereopermutators().option(stereopermutator.placement().second).value();

  if(first.placement() == indexInComponentMap.at(right)) {
    std::swap(first, second);
  }

  auto getSiteIndexLeft = [&](const PrivateGraph::Vertex i) -> SiteIndex {
    return first.getRanking().getSiteIndexOf(indexInComponentMap.at(i));
  };
  auto getSiteIndexRight = [&](const PrivateGraph::Vertex i) -> SiteIndex {
    return second.getRanking().getSiteIndexOf(indexInComponentMap.at(i));
  };

  assert(stereopermutator.numAssignments() == 2);
  for(unsigned i = 0; i < 2; ++i) {
    stereopermutator.assign(i);

    auto upOfLeftSiteIndex = temple::optionals::map(upOfLeft, getSiteIndexLeft);
    auto downOfLeftSiteIndex = temple::optionals::map(downOfLeft, getSiteIndexLeft);
    auto upOfRightSiteIndex = temple::optionals::map(upOfRight, getSiteIndexRight);
    auto downOfRightSiteIndex = temple::optionals::map(downOfRight, getSiteIndexRight);

    if(upOfLeftSiteIndex) {
      if(upOfRightSiteIndex) {
        if(std::fabs(stereopermutator.dihedral(first, *upOfLeftSiteIndex, second, *upOfRightSiteIndex)) > 1e-3) {
          continue;
        }
      }
      if(downOfRightSiteIndex) {
        if(std::fabs(stereopermutator.dihedral(first, *upOfLeftSiteIndex, second, *downOfRightSiteIndex) - M_PI) > 1e-3) {
          continue;
        }
      }
    }
    if(downOfLeftSiteIndex) {
      if(upOfRightSiteIndex) {
        if(std::fabs(stereopermutator.dihedral(first, *downOfLeftSiteIndex, second, *upOfRightSiteIndex) - M_PI) > 1e-3) {
          continue;
        }
      }
      if(downOfRightSiteIndex) {
        if(std::fabs(stereopermutator.dihedral(first, *downOfLeftSiteIndex, second, *downOfRightSiteIndex)) > 1e-3) {
          continue;
        }
      }
    }

    return i;
  }

  throw std::logic_error("Failed to find matching stereopermutation for BondStereo state.");
}

} // namespace io
} // namespace molassembler
} // namespace Scine
