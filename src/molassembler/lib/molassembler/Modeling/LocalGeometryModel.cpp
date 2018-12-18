/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Modeling/LocalGeometryModel.h"

#include "boost/range/iterator_range_core.hpp"
#include "boost/optional.hpp"
#include "Utils/ElementInfo.h"

#include "chemical_symmetries/Symmetries.h"
#include "temple/Functional.h"
#include "temple/Optionals.h"

#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/OuterGraph.h"

namespace Scine {

namespace molassembler {

namespace LocalGeometry {

const std::map<BondType, double> bondWeights {
  {BondType::Single, 1.0},
  {BondType::Double, 2.0},
  {BondType::Triple, 3.0},
  {BondType::Quadruple, 4.0},
  {BondType::Quintuple, 5.0},
  {BondType::Sextuple, 6.0},
  {BondType::Eta, 0.0}
};

boost::optional<Symmetry::Name> vsepr(
  const Scine::Utils::ElementType centerAtomType,
  const unsigned nSites,
  const std::vector<BindingSiteInformation>& sites,
  const int formalCharge
) {
  if(nSites <= 1) {
    throw std::logic_error(
      "Don't use a model on terminal atoms! Single bonds don't "
      "have stereochemistry!"
    );
  }

  if(!AtomInfo::isMainGroupElement(centerAtomType)) {
    return boost::none;
  }

  /* Make sure the ligand set doesn't include multiple atoms on a site.
   * VSEPR shouldn't try to handle haptic ligands.
   */
  if(
    std::any_of(
      sites.begin(),
      sites.end(),
      [](const auto& ligand) -> bool {
        return ligand.elements.size() > 1;
      }
    )
  ) {
    return boost::none;
  }

  // get uncharged VE count, returns none if not a main group element
  auto VEOption = molassembler::AtomInfo::mainGroupVE(centerAtomType);

  if(!VEOption) {
    return boost::none;
  }

  /* calculate X, E (VSEPR parameters). X is the number of connected atoms and E
   * is the number of non-bonding electron pairs
   */
  const unsigned& X = nSites;
  const int E = std::ceil(
    (
      static_cast<double>(VEOption.value())
      - formalCharge
      - std::accumulate(
        sites.begin(),
        sites.end(),
        0.0,
        [](const double carry, const auto& ligand) -> double {
          return carry + bondWeights.at(ligand.bondType);
        }
      )
    ) / 2.0
  );

  if(E < 0) {
    // "For some reason, E is < 0 in VSEPR. That shouldn't happen."
    return boost::none;
  }

  using Symmetry::Name;

  const auto XESum = X + E;

  if(XESum == 2) {
    return Name::Linear;
  }

  if(XESum == 3) {
    if(X == 3) {
      return Name::TrigonalPlanar;
    }

    return Name::Bent;
  }

  if(XESum == 4) {
    if(X == 4) {
      return Name::Tetrahedral;
    }

    if(X == 3) {
      /* In VSEPR naming, this is Trigonal pyramidal, but that symmetry exists
       * too, which is trigonal bipyramidal minus an apical ligand.
       */
      return Name::CutTetrahedral;
    }

    return Name::Bent;
  }

  if(XESum == 5) {
    if(X == 5) {
      return Name::TrigonalBiPyramidal;
    }

    if(X == 4) {
      return Name::Seesaw;
    }

    if(X == 3) {
      return Name::TShaped;
    }

    return Name::Linear;
  }

  if(XESum == 6) {
    if(X == 6) {
      return Name::Octahedral;
    }

    if(X == 5) {
      return Name::SquarePyramidal;
    }

    return Name::SquarePlanar;
  }

  if(XESum == 7) {
    if(X == 7) {
      return Name::PentagonalBiPyramidal;
    }

    if(X == 6) {
      return Name::PentagonalPyramidal;
    }

    return Name::PentagonalPlanar;
  }

  if(XESum == 8) {
    return Name::SquareAntiPrismatic;
  }

  /* "Could not find a fitting symmetry for your X + E case: "
    << "X = " << X << ", E = " << E << ". Maybe your molecular graph is "
    << " too weird for VSEPR. Have another look at it.";*/

  // Never runs, for static analysis
  return boost::none;
}

boost::optional<Symmetry::Name> firstOfSize(const unsigned size) {
  // Pick the first Symmetry of fitting size
  auto findIter = std::find_if(
    std::begin(Symmetry::allNames),
    std::end(Symmetry::allNames),
    [&size](const auto symmetryName) -> bool {
      return Symmetry::size(symmetryName) == size;
    }
  );

  if(findIter == std::end(Symmetry::allNames)) {
    return boost::none;
  }

  return *findIter;
}

std::vector<LocalGeometry::BindingSiteInformation> reduceToSiteInformation(
  const OuterGraph& molGraph,
  const AtomIndex index,
  const RankingInformation& ranking
) {
  /*! @todo
   * - No L, X determination. Although, will L, X even be needed for metals?
   *   Maybe only for OZ and NVE determination...
   */
  /* VSEPR formulation is that geometry is a function of
   * - localized charge of central atom
   * - atom type of central atom, neighbors
   * - bond types to neighbors
   */

  // Ensure this is only called on non-terminal atoms
  assert(ranking.ligands.size() > 1);

  // first basic stuff for VSEPR, later L and X for transition metals
  // geometry inference does not care if the substituents are somehow
  // connected (unless in later models the entire structure is considered)
  std::vector<LocalGeometry::BindingSiteInformation> ligands;
  ligands.reserve(ranking.ligands.size());

  for(const auto& ligand : ranking.ligands) {
    ligands.emplace_back(
      LocalGeometry::BindingSiteInformation {
        0,
        0,
        temple::map(ligand, [&](const AtomIndex i) -> Scine::Utils::ElementType {
          return molGraph.elementType(i);
        }),
        molGraph.bondType(
          BondIndex {index, ligand.front()}
        )
      }
    );
  }

  return ligands;
}

int formalCharge(
  const OuterGraph& graph,
  const AtomIndex index
) {
  int formalCharge = 0;

  if(AtomInfo::isMainGroupElement(graph.elementType(index))) {
    int valenceElectrons = AtomInfo::elementData.at(
      static_cast<unsigned>(graph.elementType(index))
    ).valenceElectrons();

    for(const AtomIndex adjacent : boost::make_iterator_range(graph.adjacents(index))) {
      valenceElectrons -= bondWeights.at(
        graph.bondType(
          *graph.bond(index, adjacent)
        )
      );
    }

    if(valenceElectrons > 0) {
      // Make any electrons that we can pair off non-bonding electron pairs
      int freeElectronPairs = valenceElectrons / 2;
      valenceElectrons -= 2 * freeElectronPairs;
    }

    // Assign the result of our calculation to formal charge
    formalCharge = valenceElectrons;
  }

  return formalCharge;
}

Symmetry::Name determineLocalGeometry(
  const OuterGraph& graph,
  const AtomIndex index,
  const RankingInformation& ranking
) {
  auto ligandsVector = reduceToSiteInformation(graph, index, ranking);
  unsigned nSites = ligandsVector.size();

  auto symmetryOptional = LocalGeometry::vsepr(
    graph.elementType(index),
    nSites,
    ligandsVector,
    formalCharge(graph, index)
  ) | temple::callIfNone(LocalGeometry::firstOfSize, nSites);

  if(!symmetryOptional) {
    throw std::logic_error(
      "Could not determine a geometry! Perhaps you have more substituents "
      "than the largest symmety can handle?"
    );
  }

  return symmetryOptional.value();
}

} // namespace LocalGeometry

} // namespace molassembler

} // namespace Scine
