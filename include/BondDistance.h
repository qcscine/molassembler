#ifndef INCLUDE_BOND_DISTANCE_H
#define INCLUDE_BOND_DISTANCE_H

#include <cmath>
#include "ElementTypes.h" // Delib

#include "common_typedefs.h"

namespace MoleculeManip {

namespace Bond {

/* From the original UFF paper
 * Rappé, Goddard et al. UFF, a full periodic table force field for ...
 */
static const std::map<Delib::ElementType, double> bondRadii({
  {Delib::ElementType::none, 0},
  {Delib::ElementType::H, 0.354}, // H_
  {Delib::ElementType::He, 0.849},
  {Delib::ElementType::Li,  1.336},
  {Delib::ElementType::Be, 1.074},
  {Delib::ElementType::B, 0.833}, // average(B_2, B_3)
  {Delib::ElementType::C, 0.828}, // C_3
  {Delib::ElementType::N, 0.7}, // N_3
  {Delib::ElementType::O, 0.658}, // O_3
  {Delib::ElementType::F, 0.668},
  {Delib::ElementType::Ne, 0.92},
  {Delib::ElementType::Na, 1.539},
  {Delib::ElementType::Mg, 1.421},
  {Delib::ElementType::Al, 1.244},
  {Delib::ElementType::Si, 1.117},
  {Delib::ElementType::P, 1.101}, // P_3+3
  {Delib::ElementType::S, 1.064}, // S_3+2
  {Delib::ElementType::Cl, 1.044},
  {Delib::ElementType::Ar, 1.032},
  {Delib::ElementType::K, 1.953},
  {Delib::ElementType::Ca, 1.761},
  {Delib::ElementType::Sc, 1.513},
  {Delib::ElementType::Ti, 1.412},
  {Delib::ElementType::V, 1.402},
  {Delib::ElementType::Cr, 1.345},
  {Delib::ElementType::Mn, 1.382},
  {Delib::ElementType::Fe, 1.270},
  {Delib::ElementType::Co, 1.241},
  {Delib::ElementType::Ni, 1.164},
  {Delib::ElementType::Cu, 1.302},
  {Delib::ElementType::Zn, 1.193},
  {Delib::ElementType::Ga, 1.26},
  {Delib::ElementType::Ge, 1.197},
  {Delib::ElementType::As, 1.211},
  {Delib::ElementType::Se, 1.19},
  {Delib::ElementType::Br, 1.192},
  {Delib::ElementType::Kr, 1.147},
  {Delib::ElementType::Rb, 2.260},
  {Delib::ElementType::Sr, 2.052},
  {Delib::ElementType::Y, 1.698},
  {Delib::ElementType::Zr, 1.564},
  {Delib::ElementType::Nb, 1.473},
  {Delib::ElementType::Mo, 1.467}, // Mo6+6
  {Delib::ElementType::Tc, 1.322},
  {Delib::ElementType::Ru, 1.478},
  {Delib::ElementType::Rh, 1.332},
  {Delib::ElementType::Pd, 1.338},
  {Delib::ElementType::Ag, 1.386},
  {Delib::ElementType::Cd, 1.403},
  {Delib::ElementType::In, 1.459},
  {Delib::ElementType::Sn, 1.398},
  {Delib::ElementType::Sb, 1.407},
  {Delib::ElementType::Te, 1.386},
  {Delib::ElementType::I, 1.382},
  {Delib::ElementType::Xe, 1.267},
  {Delib::ElementType::Cs, 2.570},
  {Delib::ElementType::Ba, 2.277},
  {Delib::ElementType::La, 1.943},
  {Delib::ElementType::Ce, 1.841},
  {Delib::ElementType::Pr, 1.823},
  {Delib::ElementType::Nd, 1.816},
  {Delib::ElementType::Pm, 1.801},
  {Delib::ElementType::Sm, 1.780},
  {Delib::ElementType::Eu, 1.771},
  {Delib::ElementType::Gd, 1.735},
  {Delib::ElementType::Tb, 1.732},
  {Delib::ElementType::Dy, 1.710},
  {Delib::ElementType::Ho, 1.696},
  {Delib::ElementType::Er, 1.673},
  {Delib::ElementType::Tm, 1.660},
  {Delib::ElementType::Yb, 1.637},
  {Delib::ElementType::Lu, 1.671},
  {Delib::ElementType::Hf, 1.611},
  {Delib::ElementType::Ta, 1.511},
  {Delib::ElementType::W, 1.392}, // W_6+6
  {Delib::ElementType::Re, 1.372}, // Re_6+5
  {Delib::ElementType::Os, 1.372},
  {Delib::ElementType::Ir, 1.371},
  {Delib::ElementType::Pt, 1.364},
  {Delib::ElementType::Au, 1.262},
  {Delib::ElementType::Hg, 1.340},
  {Delib::ElementType::Tl, 1.518},
  {Delib::ElementType::Pb, 1.459},
  {Delib::ElementType::Bi, 1.512},
  {Delib::ElementType::Po, 1.5},
  {Delib::ElementType::At, 1.545},
  {Delib::ElementType::Rn, 1.420},
  {Delib::ElementType::Fr, 2.880},
  {Delib::ElementType::Ra, 2.512},
  {Delib::ElementType::Ac, 1.983},
  {Delib::ElementType::Th, 1.721},
  {Delib::ElementType::Pa, 1.711},
  {Delib::ElementType::U, 1.684},
  {Delib::ElementType::Np, 1.666},
  {Delib::ElementType::Pu, 1.657},
  {Delib::ElementType::Am, 1.660},
  {Delib::ElementType::Cm, 1.801},
  {Delib::ElementType::Bk, 1.761},
  {Delib::ElementType::Cf, 1.750},
  {Delib::ElementType::Es, 1.724},
  {Delib::ElementType::Fm, 1.712},
  {Delib::ElementType::Md, 1.689},
  {Delib::ElementType::No, 1.679},
  {Delib::ElementType::Lr, 1.698},
  {Delib::ElementType::Rf, 1.6}, // NO PARAMETERS from here to Mt
  {Delib::ElementType::Db, 1.6},
  {Delib::ElementType::Sg, 1.6},
  {Delib::ElementType::Bh, 1.6},
  {Delib::ElementType::Hs, 1.6},
  {Delib::ElementType::Mt, 1.6}
});

static const std::map<BondType, double> bondOrderMap({
  {BondType::Single, 1},
  {BondType::Double, 2},
  {BondType::Triple, 3},
  {BondType::Quadruple, 4},
  {BondType::Quintuple, 5},
  {BondType::Sextuple, 6},
  {BondType::Aromatic,  1.5}
});

static const double bondOrderCorrectionLambda = 0.1332;

static double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  return (
    bondRadii.at(a)
    + bondRadii.at(b)
    - ( // BO correction
      bondOrderCorrectionLambda
      * (
        bondRadii.at(a)
        + bondRadii.at(b)
      )
      * log( bondOrderMap.at(bondType) )
    )
  );
}

} // eo namespace Bond

} // eo namespace

#endif

