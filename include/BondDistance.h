#ifndef INCLUDE_BOND_DISTANCE_H
#define INCLUDE_BOND_DISTANCE_H

#include "AtomInfo.h"
#include "common_typedefs.h"

#include <cmath>
#include <array>

/*! @file
 *
 * Bond distance modelling functions.
 */

namespace MoleculeManip {

namespace Bond {

//! Bond order definition for bond types as defined in common_typedefs
static constexpr std::array<double, 8> bondOrderMap {{
  1, 2, 3, 4, 5, 6, 1.5, 0.5
}};

//! UFF bond distance correction constant lambda
constexpr double bondOrderCorrectionLambda = 0.1332;

/*!
 * Calculate the bond distance as modelled in the original UFF paper:
 *
 * TODO FIX CITATION
 * Rappé, Goddard et al. UFF, a full periodic table force field for ...
 */
double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
);

} // namespace Bond

} // namespace MoleculeManip

#endif

