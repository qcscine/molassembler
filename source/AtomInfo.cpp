#include "AtomInfo.h"

namespace MoleculeManip {

namespace AtomInfo {

/* From the original UFF paper
 * Rappé, Goddard et al. UFF, a full periodic table force field for ...
 */
const std::map<Delib::ElementType, double> bondRadii {
  // in Angstrom
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
};

/* ElementData is populated with the following data:
 * - VdW radius, from online CRC Handbook of Chemistry and Physics, Nov. 16,
 *   97. ed.
 * - Electron occupation of orbital types s, p, d, f above largest noble core, 
 *   populated using Madelung's rule and corrected using Meek, Allen: 
 *   Configuration Irregularities, Chem Phys Lett 362, 5-6, August 2002
 *   http://dx.doi.org/10.1016/S0009-2614(02)00919-3
 */
std::map<
  Delib::ElementType,
  ElementInfo
> elementData { // 
  // in Angstrom
  {Delib::ElementType::H,  {1.09, 1, 0, 0, 0}},
  {Delib::ElementType::He, {1.40, 2, 0, 0, 0}},
  {Delib::ElementType::Li, {1.82, 1, 0, 0, 0}},
  {Delib::ElementType::Be, {1.53, 2, 0, 0, 0}},
  {Delib::ElementType::B,  {1.92, 2, 1, 0, 0}},
  {Delib::ElementType::C,  {1.70, 2, 2, 0, 0}},
  {Delib::ElementType::N,  {1.55, 2, 3, 0, 0}},
  {Delib::ElementType::O,  {1.52, 2, 4, 0, 0}},
  {Delib::ElementType::F,  {1.47, 2, 5, 0, 0}},
  {Delib::ElementType::Ne, {1.54, 2, 6, 0, 0}},
  {Delib::ElementType::Na, {2.27, 1, 0, 0, 0}},
  {Delib::ElementType::Mg, {1.73, 2, 0, 0, 0}},
  {Delib::ElementType::Al, {1.84, 2, 1, 0, 0}},
  {Delib::ElementType::Si, {2.10, 2, 2, 0, 0}},
  {Delib::ElementType::P,  {1.80, 2, 3, 0, 0}},
  {Delib::ElementType::S,  {1.80, 2, 4, 0, 0}},
  {Delib::ElementType::Cl, {1.75, 2, 5, 0, 0}},
  {Delib::ElementType::Ar, {1.88, 2, 6, 0, 0}},
  {Delib::ElementType::K,  {2.75, 1, 0, 0, 0}},
  {Delib::ElementType::Ca, {2.31, 2, 0, 0, 0}},
  {Delib::ElementType::Sc, {2.15, 2, 0, 1, 0}},
  {Delib::ElementType::Ti, {2.11, 2, 0, 2, 0}},
  {Delib::ElementType::V,  {2.07, 2, 0, 3, 0}},
  {Delib::ElementType::Cr, {2.06, 1, 0, 5, 0}}, // Madelung: 2, 0, 4, 0
  {Delib::ElementType::Mn, {2.05, 2, 0, 5, 0}},
  {Delib::ElementType::Fe, {2.04, 2, 0, 6, 0}},
  {Delib::ElementType::Co, {2.00, 2, 0, 7, 0}},
  {Delib::ElementType::Ni, {1.97, 2, 0, 8, 0}},
  {Delib::ElementType::Cu, {1.96, 1, 0, 10, 0}}, // Madelung: 2, 0, 9, 0
  {Delib::ElementType::Zn, {2.01, 2, 0, 10, 0}},
  {Delib::ElementType::Ga, {1.87, 2, 1, 10, 0}},
  {Delib::ElementType::Ge, {2.11, 2, 2, 10, 0}},
  {Delib::ElementType::As, {1.85, 2, 3, 10, 0}},
  {Delib::ElementType::Se, {1.90, 2, 4, 10, 0}},
  {Delib::ElementType::Br, {1.85, 2, 5, 10, 0}},
  {Delib::ElementType::Kr, {2.02, 2, 6, 10, 0}},
  {Delib::ElementType::Rb, {3.03, 1, 0, 0, 0}},
  {Delib::ElementType::Sr, {2.49, 2, 0, 0, 0}},
  {Delib::ElementType::Y,  {2.32, 2, 0, 1, 0}},
  {Delib::ElementType::Zr, {2.23, 2, 0, 2, 0}},
  {Delib::ElementType::Nb, {2.18, 1, 0, 4, 0}}, // Madelung: 2, 0, 3, 0
  {Delib::ElementType::Mo, {2.17, 1, 0, 5, 0}}, // Madelung: 2, 0, 4, 0
  {Delib::ElementType::Tc, {2.16, 2, 0, 5, 0}},
  {Delib::ElementType::Ru, {2.13, 1, 0, 7, 0}}, // Madelung: 2, 0, 6, 0
  {Delib::ElementType::Rh, {2.10, 1, 0, 8, 0}}, // Madelung: 2, 0, 7, 0
  {Delib::ElementType::Pd, {2.10, 0, 0, 10, 0}}, // Madelung: 2, 0, 8, 0
  {Delib::ElementType::Ag, {2.11, 1, 0, 10, 0}}, // Madelung: 2, 0, 9, 0
  {Delib::ElementType::Cd, {2.18, 2, 0, 10, 0}},
  {Delib::ElementType::In, {1.93, 2, 1, 10, 0}},
  {Delib::ElementType::Sn, {2.17, 2, 2, 10, 0}},
  {Delib::ElementType::Sb, {2.06, 2, 3, 10, 0}},
  {Delib::ElementType::Te, {2.06, 2, 4, 10, 0}},
  {Delib::ElementType::I,  {1.98, 2, 5, 10, 0}},
  {Delib::ElementType::Xe, {2.16, 2, 6, 10, 0}},
  {Delib::ElementType::Cs, {3.43, 1, 0, 0, 0}},
  {Delib::ElementType::Ba, {2.68, 2, 0, 0, 0}},
  // La: Madelung: [Xe] 6s² 4f¹ -> 2, 0, 0, 1, found [Xe] 6s² 5d¹
  {Delib::ElementType::La, {2.43, 2, 0, 1, 0}},
  // Ce: Madelung: [Xe] 6s² 4f² -> 2, 0, 0, 2, found [Xe] 6s² 4f¹ 5d¹
  {Delib::ElementType::Ce, {2.42, 2, 0, 1, 1}},
  {Delib::ElementType::Pr, {2.40, 2, 0, 0, 3}},
  {Delib::ElementType::Nd, {2.39, 2, 0, 0, 4}},
  {Delib::ElementType::Pm, {2.38, 2, 0, 0, 5}},
  {Delib::ElementType::Sm, {2.36, 2, 0, 0, 6}},
  {Delib::ElementType::Eu, {2.35, 2, 0, 0, 7}},
  // Gd: Madelung: [Xe] 6s² 4f⁸ -> 2, 0, 0, 8, found [Xe] 6s² 4f⁷ 5d¹
  {Delib::ElementType::Gd, {2.34, 2, 0, 1, 7}},
  {Delib::ElementType::Tb, {2.33, 2, 0, 0, 9}},
  {Delib::ElementType::Dy, {2.31, 2, 0, 0, 10}},
  {Delib::ElementType::Ho, {2.30, 2, 0, 0, 11}},
  {Delib::ElementType::Er, {2.29, 2, 0, 0, 12}},
  {Delib::ElementType::Tm, {2.27, 2, 0, 0, 13}},
  {Delib::ElementType::Yb, {2.26, 2, 0, 0, 14}},
  {Delib::ElementType::Lu, {2.24, 2, 0, 1, 14}},
  {Delib::ElementType::Hf, {2.23, 2, 0, 2, 14}},
  {Delib::ElementType::Ta, {2.22, 2, 0, 3, 14}},
  {Delib::ElementType::W,  {2.18, 2, 0, 4, 14}},
  {Delib::ElementType::Re, {2.16, 2, 0, 5, 14}},
  {Delib::ElementType::Os, {2.16, 2, 0, 6, 14}},
  {Delib::ElementType::Ir, {2.13, 2, 0, 7, 14}},
  // Pt: Madelung: [Xe] 6s² (4f¹⁴) 5d⁸ -> 2, 0, 8, 14, found [Xe] 6s¹ (4f¹⁴) 5d⁹
  {Delib::ElementType::Pt, {2.13, 1, 0, 9, 14}},
  // Au: Madelung: [Xe] 6s² (4f¹⁴) 5d⁹ -> 2, 0, 9, 14, found [Xe] 6s¹ (4f¹⁴) 5d¹⁰
  {Delib::ElementType::Au, {2.14, 1, 0, 10, 14}},
  {Delib::ElementType::Hg, {2.23, 2, 0, 10, 14}},
  {Delib::ElementType::Tl, {1.96, 2, 1, 10, 14}},
  {Delib::ElementType::Pb, {2.02, 2, 2, 10, 14}},
  {Delib::ElementType::Bi, {2.07, 2, 3, 10, 14}},
  {Delib::ElementType::Po, {1.97, 2, 4, 10, 14}},
  {Delib::ElementType::At, {2.02, 2, 5, 10, 14}},
  {Delib::ElementType::Rn, {2.20, 2, 6, 10, 14}},
  {Delib::ElementType::Fr, {3.48, 1, 0, 0, 0}},
  {Delib::ElementType::Ra, {2.83, 2, 0, 0, 0}},
  // Ac: Madelung: [Rn] 7s² 5f¹ -> 2, 0, 0, 1, found [Rn] 7s² 6d¹
  {Delib::ElementType::Ac, {2.47, 2, 0, 1, 0}},
  // Th: Madelung: [Rn] 7s² 5f² -> 2, 0, 0, 2, found [Rn] 7s² 6d²
  {Delib::ElementType::Th, {2.45, 2, 0, 2, 0}},
  // Pa: Madelung: [Rn] 7s² 5f³ -> 2, 0, 0, 3, found [Rn] 7s² 5f² 6d¹
  {Delib::ElementType::Pa, {2.43, 2, 0, 1, 2}},
  //  U: Madelung: [Rn] 7s² 5f⁴ -> 2, 0, 0, 4, found [Rn] 7s² 5f³ 6d¹
  {Delib::ElementType::U,  {2.41, 2, 0, 3, 1}},
  // Np: Madelung: [Rn] 7s² 5f⁵ -> 2, 0, 0, 5, found [Rn] 7s² 5f⁴ 6d¹
  {Delib::ElementType::Np, {2.39, 2, 0, 4, 1}},
  {Delib::ElementType::Pu, {2.43, 2, 0, 0, 6}},
  {Delib::ElementType::Am, {2.44, 2, 0, 0, 7}},
  // Cm: Madelung: [Rn] 7s² 5f⁸ -> 2, 0, 0, 8, found [Rn] 7s² 5f⁷ 6d¹
  {Delib::ElementType::Cm, {2.45, 2, 0, 7, 1}},
  {Delib::ElementType::Bk, {2.44, 2, 0, 0, 9}},
  {Delib::ElementType::Cf, {2.45, 2, 0, 0, 10}},
  {Delib::ElementType::Es, {2.45, 2, 0, 0, 11}},
  {Delib::ElementType::Fm, {2.45, 2, 0, 0, 12}},
  {Delib::ElementType::Md, {2.46, 2, 0, 0, 13}},
  {Delib::ElementType::No, {2.46, 2, 0, 0, 14}},
  {Delib::ElementType::Lr, {2.46, 2, 0, 1, 14}},
  {Delib::ElementType::Rf, {0.0, 2, 0, 2, 14}},
  {Delib::ElementType::Db, {0.0, 2, 0, 3, 14}},
  {Delib::ElementType::Sg, {0.0, 2, 0, 4, 14}},
  {Delib::ElementType::Bh, {0.0, 2, 0, 5, 14}},
  {Delib::ElementType::Hs, {0.0, 2, 0, 6, 14}},
  {Delib::ElementType::Mt, {0.0, 2, 0, 7, 14}}
};

bool isMainGroupElement(const Delib::ElementType& elementType) {
  return elementData.at(elementType).shellsFullOrEmpty({'d', 'f'});
}

boost::optional<unsigned> mainGroupVE(const Delib::ElementType& elementType) {
  if(isMainGroupElement(elementType)) {
    return elementData.at(elementType).valenceElectrons({'s', 'p'});
  } else {
    return {};
  }
}

unsigned dElectronCount(const Delib::ElementType& elementType) {
  if(isMainGroupElement(elementType)) {
    return 0;
  } else {
    return elementData.at(elementType).valenceElectrons('d');
  }
}

double vdwRadius(const Delib::ElementType& elementType) {
  return elementData.at(elementType).vdwRadius;
}

} // namespace AtomInfo

} // namespace MoleculeManip
