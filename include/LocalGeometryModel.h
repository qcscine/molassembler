#ifndef INCLUDE_GEOMETRY_MODELING_MODEL_H
#define INCLUDE_GEOMETRY_MODELING_MODEL_H

#include "Delib/ElementTypes.h"
#include "symmetry_information/Symmetries.h"

#include "common_typedefs.h"
#include "boost/variant.hpp"

namespace LocalGeometry {

/* Typedefs */

using BondType = MoleculeManip::BondType;

using ElementAndBondPair = std::pair<
  Delib::ElementType,
  BondType
>;

using LigandType = std::tuple<
  unsigned, // L
  unsigned, // X
  std::vector<
    ElementAndBondPair
  >
>;

// Mapping of bond type to a floating-point weight
extern const std::map<BondType, double> bondWeights;

/* Model class */
struct Model {
  // do not instantiate
  Model() = delete;

  static Symmetry::Name determineGeometry(
    const Delib::ElementType& centerAtomType,
    const unsigned& nSites,
    const std::vector<LigandType>& ligands,
    const int& charge
  );
};

} // eo namespace GeometryModeling

#endif
