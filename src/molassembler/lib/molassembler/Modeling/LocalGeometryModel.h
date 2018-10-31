// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_LOCAL_GEOMETRY_MODEL_H
#define INCLUDE_MOLASSEMBLER_LOCAL_GEOMETRY_MODEL_H

#include "boost/optional/optional_fwd.hpp"

#include "Delib/ElementTypes.h"
#include "chemical_symmetries/Names.h"

#include "molassembler/RankingInformation.h"

#include <map>

/*! @file
 *
 * @brief Algorithms to determine local symmetry from graph information
 *
 * Declarations for the general interface with which a number of classes can
 * determine the local geometry that a specific arrangement of atoms should
 * have.
 */

namespace molassembler {

// Forward-declarations
class OuterGraph;

namespace LocalGeometry {

/* Typedefs */
struct BindingSiteInformation {
  unsigned L, X;

  std::vector<Delib::ElementType> elements;
  /* Only one bond type is needed - If the ligand consists of a single atom,
   * then we only need to store one bond. If the ligand consists of multiple
   * atoms, then the BondType is Eta.
   */
  BondType bondType;

  BindingSiteInformation() = default;
  BindingSiteInformation(
    const unsigned passL,
    const unsigned passX,
    std::vector<Delib::ElementType> passElements,
    const BondType passBondType
  ) : L(passL), X(passX), elements {std::move(passElements)}, bondType {passBondType} {}
};

// Mapping of bond type to a floating-point weight
extern const std::map<BondType, double> bondWeights;

/**
 * @brief Calculates the formal charge on a main group-element atom.
 *
 * @param graph The graph to which the atom index belongs
 * @param index The atom index for which to calculate the formal charge
 *
 * @warning This is an awful function and should be avoided in any sort of
 *   important calculation.
 *
 * @warning This will yield nonsense if the bond orders in your graph are unset.
 *
 * @warning This function may work sometimes for organic surroundings. That's
 *   how confident we are in this function.
 *
 * @return 0 for non-main group elements, the formal charge otherwise
 */
int formalCharge(
  const OuterGraph& graph,
  AtomIndex index
);

/* Models */
boost::optional<Symmetry::Name> vsepr(
  Delib::ElementType centerAtomType,
  unsigned nSites,
  const std::vector<BindingSiteInformation>& sites,
  int formalCharge
);

boost::optional<Symmetry::Name> firstOfSize(unsigned size);


/* Tiered geometry determination function */
std::vector<LocalGeometry::BindingSiteInformation> reduceToSiteInformation(
  const OuterGraph& molGraph,
  AtomIndex index,
  const RankingInformation& ranking
);

Symmetry::Name determineLocalGeometry(
  const OuterGraph& graph,
  AtomIndex index,
  const RankingInformation& ranking
);

} // namespace LocalGeometry

} // namespace molassembler

#endif
