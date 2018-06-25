#ifndef INCLUDE_GEOMETRY_MODELING_MODEL_H
#define INCLUDE_GEOMETRY_MODELING_MODEL_H

#include "chemical_symmetries/Symmetries.h"

#include "RankingInformation.h"

/*! @file
 *
 * Declarations for the general interface with which a number of classes can
 * determine the local geometry that a specific arrangement of atoms should
 * have.
 */

namespace molassembler {

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
    unsigned L,
    unsigned X,
    std::vector<Delib::ElementType> elements,
    BondType bondType
  ) : L {L}, X {X}, elements {elements}, bondType {bondType} {}
};

// Mapping of bond type to a floating-point weight
extern const std::map<BondType, double> bondWeights;

/* Models */
boost::optional<Symmetry::Name> vsepr(
  const Delib::ElementType centerAtomType,
  const unsigned nSites,
  const std::vector<BindingSiteInformation>& sites,
  const int formalCharge
);

boost::optional<Symmetry::Name> firstOfSize(const unsigned& size);


/* Tiered geometry determination function */
std::vector<LocalGeometry::BindingSiteInformation> reduceToSiteInformation(
  const GraphType& molGraph,
  const AtomIndexType index,
  const RankingInformation& ranking
);

Symmetry::Name determineLocalGeometry(
  const GraphType& graph,
  const AtomIndexType index,
  const RankingInformation& ranking
);

} // namespace LocalGeometry

} // namespace molassembler

#endif
