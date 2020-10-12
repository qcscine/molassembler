/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Algorithms to determine local shape from graph information
 *
 * Declarations for the general interface with which a number of classes can
 * determine the local geometry that a specific arrangement of atoms should
 * have.
 */

#ifndef INCLUDE_MOLASSEMBLER_LOCAL_GEOMETRY_MODEL_H
#define INCLUDE_MOLASSEMBLER_LOCAL_GEOMETRY_MODEL_H

#include "boost/optional/optional_fwd.hpp"

#include "Utils/Geometry/ElementTypes.h"
#include "Molassembler/Shapes/Shapes.h"

#include "Molassembler/RankingInformation.h"

#include <map>

namespace Scine {
namespace Molassembler {

// Forward-declarations
class Graph;

/**
 * @brief Methods to determine local shapes of atoms based on graph information
 */
namespace ShapeInference {

/**
 * @brief Type used to represent minimal binding site information
 */
struct BindingSite {
  unsigned L, X;

  std::vector<Utils::ElementType> elements;
  /* Only one bond type is needed - If the ligand consists of a single atom,
   * then we only need to store one bond. If the ligand consists of multiple
   * atoms, then the BondType is Eta.
   */
  BondType bondType;

  BindingSite() = default;
  BindingSite(
    const unsigned passL,
    const unsigned passX,
    std::vector<Utils::ElementType> passElements,
    const BondType passBondType
  ) : L(passL), X(passX), elements(std::move(passElements)), bondType(passBondType) {}
};

/*! @brief Mapping of bond type to a floating-point weight
 *
 * @todo Eta bonds have weight 0, is this a good idea?
 */
extern const std::map<BondType, double> bondWeights;

/** @brief Calculates the formal charge on a main group-element atom.
 *
 * @complexity{@math{\Theta(1)}}
 *
 * @param graph The graph to which the atom index belongs
 * @param index The atom index for which to calculate the formal charge
 *
 * @parblock@warning This is an awful concept and should be avoided in any
 * sort of important calculation.
 * @endparblock
 *
 * @parblock@warning This will yield nonsense if the bond orders in your graph
 * are not correct.
 * @endparblock
 *
 * @parblock@warning This function may work sometimes for organic
 * neighborhoods. That's how confident we are in this function.
 * @endparblock
 *
 * @return 0 for non-main group elements, possibly the formal charge otherwise
 */
int formalCharge(
  const Graph& graph,
  AtomIndex index
);

/* Models */
/*! @brief Applies very basic VSEPR theory to derive a shape based on graph
 *   information
 *
 * @complexity{@math{\Theta(1)}}
 */
boost::optional<Shapes::Shape> vsepr(
  Utils::ElementType centerAtomType,
  const std::vector<BindingSite>& sites,
  int formalCharge
);

/*! @brief Yields the first shape of required size.
 *
 * @complexity{@math{\Theta(1)}}
 * @throws std::runtime_error If no shapes of @p size exist
 */
Shapes::Shape firstOfSize(unsigned size);

/*! @brief Reduces a ranking to binding site information
 */
std::vector<BindingSite> reduceToSiteInformation(
  const Graph& molGraph,
  AtomIndex index,
  const RankingInformation& ranking
);

/*! @brief Forwards inference to appropriate model depending on environment
 *
 * Currently just a stub forwarding to VSEPR
 */
boost::optional<Shapes::Shape> inferShape(
  const Graph& graph,
  AtomIndex index,
  const RankingInformation& ranking
);

} // namespace ShapeInference
} // namespace Molassembler
} // namespace Scine

#endif
