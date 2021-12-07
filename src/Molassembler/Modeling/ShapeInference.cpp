/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Modeling/ShapeInference.h"

#include "boost/optional.hpp"
#include "Utils/Geometry/ElementInfo.h"

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Temple/Functional.h"

#include "Molassembler/Modeling/AtomInfo.h"
#include "Molassembler/Graph.h"

namespace Scine {
namespace Molassembler {
namespace ShapeInference {

double bondWeight(const BondType bond) {
  switch(bond) {
    case BondType::Single: return 1.0;
    case BondType::Double: return 2.0;
    case BondType::Triple: return 3.0;
    case BondType::Quadruple: return 4.0;
    case BondType::Quintuple: return 5.0;
    case BondType::Sextuple: return 6.0;
    default: return 0.0;
  }
}

boost::optional<Shapes::Shape> vsepr(
  const Utils::ElementType centerAtomType,
  const std::vector<BindingSite>& sites,
  const int formalCharge
) {
  const unsigned nSites = sites.size();

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
  auto veOption = Molassembler::AtomInfo::mainGroupVE(centerAtomType);

  if(!veOption) {
    return boost::none;
  }

  /* calculate X, E (VSEPR parameters). X is the number of connected atoms and E
   * is the number of non-bonding electron pairs
   */
  const unsigned X = nSites;
  const int E = std::ceil(
    (
      static_cast<double>(veOption.value())
      - formalCharge
      - std::accumulate(
        sites.begin(),
        sites.end(),
        0.0,
        [](const double carry, const auto& ligand) -> double {
          return carry + bondWeight(ligand.bondType);
        }
      )
    ) / 2.0
  );

  if(E < 0) {
    // "For some reason, E is < 0 in VSEPR. That shouldn't happen."
    return boost::none;
  }

  using Shapes::Shape;

  switch(X + E) {
    case 2:
      return Shape::Line;

    case 3:
      switch(X) {
        case 3: return Shape::EquilateralTriangle;
        default: return Shape::Bent;
      }

    case 4:
      switch(X) {
        case 4: return Shape::Tetrahedron;
        case 3: return Shape::VacantTetrahedron;
        default: return Shape::Bent;
      }

    case 5:
      switch(X) {
        case 5: return Shape::TrigonalBipyramid;
        case 4: return Shape::Seesaw;
        case 3: return Shape::T;
        default: return Shape::Line;
      }

    case 6:
      switch(X) {
        case 6: return Shape::Octahedron;
        case 5: return Shape::SquarePyramid;
        default: return Shape::Square;
      }

    case 7:
      switch(X) {
        case 7: return Shape::PentagonalBipyramid;
        case 6: return Shape::PentagonalPyramid;
        default: return Shape::Pentagon;
      }

    case 8:
      return Shape::SquareAntiprism;

    default: return boost::none;
  }
}

Shapes::Shape firstOfSize(const unsigned size) {
  // Pick the first shape of fitting size
  const auto findIter = std::find_if(
    std::begin(Shapes::allShapes),
    std::end(Shapes::allShapes),
    [&size](const auto shape) -> bool {
      return Shapes::size(shape) == size;
    }
  );

  if(findIter == std::end(Shapes::allShapes)) {
    throw std::runtime_error("No shapes of that size!");
  }

  return *findIter;
}

std::vector<BindingSite> reduceToSiteInformation(
  const Graph& molGraph,
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
  assert(ranking.sites.size() > 1);

  // first basic stuff for VSEPR, later L and X for transition metals
  // geometry inference does not care if the substituents are somehow
  // connected (unless in later models the entire structure is considered)
  std::vector<BindingSite> sites;
  sites.reserve(ranking.sites.size());

  for(const auto& ligand : ranking.sites) {
    sites.emplace_back(
      BindingSite {
        0,
        0,
        Temple::map(ligand, [&](const AtomIndex i) -> Utils::ElementType {
          return molGraph.elementType(i);
        }),
        molGraph.bondType(
          BondIndex {index, ligand.front()}
        )
      }
    );
  }

  return sites;
}

int formalCharge(
  const Graph& graph,
  const AtomIndex index
) {
  int formalCharge = 0;

  if(AtomInfo::isMainGroupElement(graph.elementType(index))) {
    int valenceElectrons = AtomInfo::elementData().at(
      Utils::ElementInfo::Z(graph.elementType(index))
    ).valenceElectrons();

    for(const AtomIndex adjacent : graph.adjacents(index)) {
      valenceElectrons -= static_cast<int>(
        bondWeight(
          graph.bondType(*graph.bond(index, adjacent))
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

boost::optional<Shapes::Shape> inferShape(
  const Graph& graph,
  const AtomIndex index,
  const RankingInformation& ranking
) {
  // So long as only VSEPR is implemented, this function is very simple:
  return vsepr(
    graph.elementType(index),
    reduceToSiteInformation(graph, index, ranking),
    formalCharge(graph, index)
  );
}

} // namespace ShapeInference
} // namespace Molassembler
} // namespace Scine
