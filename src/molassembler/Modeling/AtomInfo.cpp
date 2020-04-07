/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Modeling/AtomInfo.h"

#include "Utils/Geometry/ElementInfo.h"

#include "boost/optional.hpp"

namespace Scine {

namespace molassembler {

namespace atom_info {

ElementInfo::ElementInfo(
  const double passVdwRadius,
  const unsigned sValenceElectrons,
  const unsigned pValenceElectrons,
  const unsigned dValenceElectrons,
  const unsigned fValenceElectrons
) :
  valenceElectrons_ {
    sValenceElectrons,
    pValenceElectrons,
    dValenceElectrons,
    fValenceElectrons
  },
  vdwRadius_ (passVdwRadius)
{}

unsigned ElementInfo::maxOccupancy(const char shell) {
  switch(shell) {
    case 's':
      return 2;
    case 'p':
      return 6;
    case 'd':
      return 10;
    case 'f':
      return 14;
    default:
      return 0;
  }
}

unsigned ElementInfo::valenceElectrons(const char shell) const {
  switch(shell) {
    case 's':
      return valenceElectrons_[0];
    case 'p':
      return valenceElectrons_[1];
    case 'd':
      return valenceElectrons_[2];
    case 'f':
      return valenceElectrons_[3];
    default:
      return 0;
  }
}

unsigned ElementInfo::valenceElectrons(const std::vector<char>& shells) const {
  unsigned sum = 0;
  for(const char shell : shells) {
    sum += valenceElectrons(shell);
  }
  return sum;
}

bool ElementInfo::shellFullOrEmpty(const char shell) const {
  return (
    valenceElectrons(shell) == 0
    || valenceElectrons(shell) == maxOccupancy(shell)
  );
}

bool ElementInfo::shellsFullOrEmpty(const std::vector<char>& shells) const {
  for(const char shell : shells) {
    if(!shellFullOrEmpty(shell)) {
      return false;
    }
  }

  return true;
}

//! Returns the total valence electrons
unsigned ElementInfo::valenceElectrons() const {
  unsigned sum = 0;
  for(unsigned i = 0; i < 4; i++) {
    sum += valenceElectrons_[i];
  }
  return sum;
}

double ElementInfo::vdwRadius() const {
  return vdwRadius_;
}

/* From the original UFF paper
 * Rappé, Goddard et al. UFF, a full periodic table force field for ...
 */
const std::array<double, 110> bondRadii {
  0,  // none
  0.354,  // H
  0.849,  // He
  1.336,  // Li
  1.074,  // Be
  0.833,  // average(B_2 B_3)
  0.828,  // C_3
  0.7,  // N_3
  0.658,  // O_3
  0.668,  // F
  0.92,  // Ne
  1.539,  // Na
  1.421,  // Mg
  1.244,  // Al
  1.117,  // Si
  1.101,  // P_3+3
  1.064,  // S_3+2
  1.044,  // Cl
  1.032,  // Ar
  1.953,  // K
  1.761,  // Ca
  1.513,  // Sc
  1.412,  // Ti
  1.402,  // V
  1.345,  // Cr
  1.382,  // Mn
  1.270,  // Fe
  1.241,  // Co
  1.164,  // Ni
  1.302,  // Cu
  1.193,  // Zn
  1.26,  // Ga
  1.197,  // Ge
  1.211,  // As
  1.19,  // Se
  1.192,  // Br
  1.147,  // Kr
  2.260,  // Rb
  2.052,  // Sr
  1.698,  // Y
  1.564,  // Zr
  1.473,  // Nb
  1.467,  // Mo6+6
  1.322,  // Tc
  1.478,  // Ru
  1.332,  // Rh
  1.338,  // Pd
  1.386,  // Ag
  1.403,  // Cd
  1.459,  // In
  1.398,  // Sn
  1.407,  // Sb
  1.386,  // Te
  1.382,  // I
  1.267,  // Xe
  2.570,  // Cs
  2.277,  // Ba
  1.943,  // La
  1.841,  // Ce
  1.823,  // Pr
  1.816,  // Nd
  1.801,  // Pm
  1.780,  // Sm
  1.771,  // Eu
  1.735,  // Gd
  1.732,  // Tb
  1.710,  // Dy
  1.696,  // Ho
  1.673,  // Er
  1.660,  // Tm
  1.637,  // Yb
  1.671,  // Lu
  1.611,  // Hf
  1.511,  // Ta
  1.392,  // W_6+6
  1.372,  // Re_6+5
  1.372,  // Os
  1.371,  // Ir
  1.364,  // Pt
  1.262,  // Au
  1.340,  // Hg
  1.518,  // Tl
  1.459,  // Pb
  1.512,  // Bi
  1.5,  // Po
  1.545,  // At
  1.420,  // Rn
  2.880,  // Fr
  2.512,  // Ra
  1.983,  // Ac
  1.721,  // Th
  1.711,  // Pa
  1.684,  // U
  1.666,  // Np
  1.657,  // Pu
  1.660,  // Am
  1.801,  // Cm
  1.761,  // Bk
  1.750,  // Cf
  1.724,  // Es
  1.712,  // Fm
  1.689,  // Md
  1.679,  // No
  1.698,  // Lr
  1.6,  // Rf, NO PARAMETERS from here to Mt
  1.6,  // Db
  1.6,  // Sg
  1.6,  // Bh
  1.6,  // Hs
  1.6  // Mt
};

double bondRadius(const Utils::ElementType elementType) {
  return bondRadii.at(
    Utils::ElementInfo::Z(elementType)
  );
}

/* ElementData is populated with the following data:
 * - VdW radius, from online CRC Handbook of Chemistry and Physics, Nov. 16,
 *   97. ed.
 * - Electron occupation of orbital types s, p, d, f above largest noble core,
 *   populated using Madelung's rule and corrected using Meek, Allen:
 *   Configuration Irregularities, Chem Phys Lett 362, 5-6, August 2002
 *   http://dx.doi.org/10.1016/S0009-2614(02)00919-3
 */


std::array<ElementInfo, 110> elementData {{
  // in Angstrom
  {0.00, 0u, 0u, 0u, 0u}, // none
  {1.09, 1u, 0u, 0u, 0u}, // H
  {1.40, 2u, 0u, 0u, 0u}, // He
  {1.82, 1u, 0u, 0u, 0u}, // Li
  {1.53, 2u, 0u, 0u, 0u}, // Be
  {1.92, 2u, 1u, 0u, 0u}, // B
  {1.70, 2u, 2u, 0u, 0u}, // C
  {1.55, 2u, 3u, 0u, 0u}, // N
  {1.52, 2u, 4u, 0u, 0u}, // O
  {1.47, 2u, 5u, 0u, 0u}, // F
  {1.54, 2u, 6u, 0u, 0u}, // Ne
  {2.27, 1u, 0u, 0u, 0u}, // Na
  {1.73, 2u, 0u, 0u, 0u}, // Mg
  {1.84, 2u, 1u, 0u, 0u}, // Al
  {2.10, 2u, 2u, 0u, 0u}, // Si
  {1.80, 2u, 3u, 0u, 0u}, // P
  {1.80, 2u, 4u, 0u, 0u}, // S
  {1.75, 2u, 5u, 0u, 0u}, // Cl
  {1.88, 2u, 6u, 0u, 0u}, // Ar
  {2.75, 1u, 0u, 0u, 0u}, // K
  {2.31, 2u, 0u, 0u, 0u}, // Ca
  {2.15, 2u, 0u, 1u, 0u}, // Sc
  {2.11, 2u, 0u, 2u, 0u}, // Ti
  {2.07, 2u, 0u, 3u, 0u}, // V
  {2.06, 1u, 0u, 5u, 0u}, // Cr, Madelung: 2, 0, 4, 0
  {2.05, 2u, 0u, 5u, 0u}, // Mn
  {2.04, 2u, 0u, 6u, 0u}, // Fe
  {2.00, 2u, 0u, 7u, 0u}, // Co
  {1.97, 2u, 0u, 8u, 0u}, // Ni
  {1.96, 1u, 0u, 10u, 0u}, // Cu, Madelung: 2, 0, 9, 0
  {2.01, 2u, 0u, 10u, 0u}, // Zn
  {1.87, 2u, 1u, 10u, 0u}, // Ga
  {2.11, 2u, 2u, 10u, 0u}, // Ge
  {1.85, 2u, 3u, 10u, 0u}, // As
  {1.90, 2u, 4u, 10u, 0u}, // Se
  {1.85, 2u, 5u, 10u, 0u}, // Br
  {2.02, 2u, 6u, 10u, 0u}, // Kr
  {3.03, 1u, 0u, 0u, 0u}, // Rb
  {2.49, 2u, 0u, 0u, 0u}, // Sr
  {2.32, 2u, 0u, 1u, 0u}, // Y
  {2.23, 2u, 0u, 2u, 0u}, // Zr
  {2.18, 1u, 0u, 4u, 0u}, // Nb, Madelung: 2, 0, 3, 0
  {2.17, 1u, 0u, 5u, 0u}, // Mo, Madelung: 2, 0, 4, 0
  {2.16, 2u, 0u, 5u, 0u}, // Tc
  {2.13, 1u, 0u, 7u, 0u}, // Ru, Madelung: 2, 0, 6, 0
  {2.10, 1u, 0u, 8u, 0u}, // Rh, Madelung: 2, 0, 7, 0
  {2.10, 0u, 0u, 10u, 0u}, // Pd, Madelung: 2, 0, 8, 0
  {2.11, 1u, 0u, 10u, 0u}, // Ag, Madelung: 2, 0, 9, 0
  {2.18, 2u, 0u, 10u, 0u}, // Cd
  {1.93, 2u, 1u, 10u, 0u}, // In
  {2.17, 2u, 2u, 10u, 0u}, // Sn
  {2.06, 2u, 3u, 10u, 0u}, // Sb
  {2.06, 2u, 4u, 10u, 0u}, // Te
  {1.98, 2u, 5u, 10u, 0u}, // I
  {2.16, 2u, 6u, 10u, 0u}, // Xe
  {3.43, 1u, 0u, 0u, 0u}, // Cs
  {2.68, 2u, 0u, 0u, 0u}, // Ba
  //  La: Madelung: [Xe] 6s² 4f¹ -> 2, 0, 0, 1, found [Xe] 6s² 5d¹
  {2.43, 2u, 0u, 1u, 0u}, // La
  //  Ce: Madelung: [Xe] 6s² 4f² -> 2, 0, 0, 2, found [Xe] 6s² 4f¹ 5d¹
  {2.42, 2u, 0u, 1u, 1u}, // Ce
  {2.40, 2u, 0u, 0u, 3u}, // Pr
  {2.39, 2u, 0u, 0u, 4u}, // Nd
  {2.38, 2u, 0u, 0u, 5u}, // Pm
  {2.36, 2u, 0u, 0u, 6u}, // Sm
  {2.35, 2u, 0u, 0u, 7u}, // Eu
  //  Gd: Madelung: [Xe] 6s² 4f⁸ -> 2, 0, 0, 8, found [Xe] 6s² 4f⁷ 5d¹
  {2.34, 2u, 0u, 1u, 7u}, // Gd
  {2.33, 2u, 0u, 0u, 9u}, // Tb
  {2.31, 2u, 0u, 0u, 10u}, // Dy
  {2.30, 2u, 0u, 0u, 11u}, // Ho
  {2.29, 2u, 0u, 0u, 12u}, // Er
  {2.27, 2u, 0u, 0u, 13u}, // Tm
  {2.26, 2u, 0u, 0u, 14u}, // Yb
  {2.24, 2u, 0u, 1u, 14u}, // Lu
  {2.23, 2u, 0u, 2u, 14u}, // Hf
  {2.22, 2u, 0u, 3u, 14u}, // Ta
  {2.18, 2u, 0u, 4u, 14u}, // W
  {2.16, 2u, 0u, 5u, 14u}, // Re
  {2.16, 2u, 0u, 6u, 14u}, // Os
  {2.13, 2u, 0u, 7u, 14u}, // Ir
  //  Pt: Madelung: [Xe] 6s² (4f¹⁴) 5d⁸ -> 2, 0, 8, 14, found [Xe] 6s¹ (4f¹⁴) 5d⁹
  {2.13, 1u, 0u, 9u, 14u}, // Pt
  //  Au: Madelung: [Xe] 6s² (4f¹⁴) 5d⁹ -> 2, 0, 9, 14, found [Xe] 6s¹ (4f¹⁴) 5d¹⁰
  {2.14, 1u, 0u, 10u, 14u}, // Au
  {2.23, 2u, 0u, 10u, 14u}, // Hg
  {1.96, 2u, 1u, 10u, 14u}, // Tl
  {2.02, 2u, 2u, 10u, 14u}, // Pb
  {2.07, 2u, 3u, 10u, 14u}, // Bi
  {1.97, 2u, 4u, 10u, 14u}, // Po
  {2.02, 2u, 5u, 10u, 14u}, // At
  {2.20, 2u, 6u, 10u, 14u}, // Rn
  {3.48, 1u, 0u, 0u, 0u}, // Fr
  {2.83, 2u, 0u, 0u, 0u}, // Ra
  //  Ac: Madelung: [Rn] 7s² 5f¹ -> 2, 0, 0, 1, found [Rn] 7s² 6d¹
  {2.47, 2u, 0u, 1u, 0u}, // Ac
  //  Th: Madelung: [Rn] 7s² 5f² -> 2, 0, 0, 2, found [Rn] 7s² 6d²
  {2.45, 2u, 0u, 2u, 0u}, // Th
  //  Pa: Madelung: [Rn] 7s² 5f³ -> 2, 0, 0, 3, found [Rn] 7s² 5f² 6d¹
  {2.43, 2u, 0u, 1u, 2u}, // Pa
  //   U: Madelung: [Rn] 7s² 5f⁴ -> 2, 0, 0, 4, found [Rn] 7s² 5f³ 6d¹
  {2.41, 2u, 0u, 3u, 1u}, // U,
  //  Np: Madelung: [Rn] 7s² 5f⁵ -> 2, 0, 0, 5, found [Rn] 7s² 5f⁴ 6d¹
  {2.39, 2u, 0u, 4u, 1u}, // Np
  {2.43, 2u, 0u, 0u, 6u}, // Pu
  {2.44, 2u, 0u, 0u, 7u}, // Am
  //  Cm: Madelung: [Rn] 7s² 5f⁸ -> 2, 0, 0, 8, found [Rn] 7s² 5f⁷ 6d¹
  {2.45, 2u, 0u, 7u, 1u}, // Cm
  {2.44, 2u, 0u, 0u, 9u}, // Bk
  {2.45, 2u, 0u, 0u, 10u}, // Cf
  {2.45, 2u, 0u, 0u, 11u}, // Es
  {2.45, 2u, 0u, 0u, 12u}, // Fm
  {2.46, 2u, 0u, 0u, 13u}, // Md
  {2.46, 2u, 0u, 0u, 14u}, // No
  {2.46, 2u, 0u, 1u, 14u}, // Lr
  {0.0, 2u, 0u, 2u, 14u}, // Rf
  {0.0, 2u, 0u, 3u, 14u}, // Db
  {0.0, 2u, 0u, 4u, 14u}, // Sg
  {0.0, 2u, 0u, 5u, 14u}, // Bh
  {0.0, 2u, 0u, 6u, 14u}, // Hs
  {0.0, 2u, 0u, 7u, 14u} // Mt
}};

bool isMainGroupElement(const Utils::ElementType elementType) {
  unsigned Z = Utils::ElementInfo::Z(elementType);
  return (
    Z <= 20
    || (31 <= Z && Z <= 38)
    || (49 <= Z && Z <= 56)
    || (81 <= Z && Z <= 88)
    || (113 <= Z && Z <= 118)
  );
}

boost::optional<unsigned> mainGroupVE(const Utils::ElementType elementType) {
  if(isMainGroupElement(elementType)) {
    return elementData.at(
      Utils::ElementInfo::Z(elementType)
    ).valenceElectrons({'s', 'p'});
  }

  return {};
}

unsigned dElectronCount(const Utils::ElementType elementType) {
  if(isMainGroupElement(elementType)) {
    return 0;
  }

  return elementData.at(
    Utils::ElementInfo::Z(elementType)
  ).valenceElectrons('d');
}

double vdwRadius(const Utils::ElementType elementType) {
  return elementData.at(
    Utils::ElementInfo::Z(elementType)
  ).vdwRadius();
}

} // namespace atom_info

} // namespace molassembler

} // namespace Scine
