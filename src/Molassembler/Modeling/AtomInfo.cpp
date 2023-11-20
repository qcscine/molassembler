/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Modeling/AtomInfo.h"

#include "Utils/Geometry/ElementInfo.h"

#include "boost/optional.hpp"

#include <algorithm>

namespace Scine {
namespace Molassembler {
namespace AtomInfo {

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
  return std::all_of(
    std::begin(shells),
    std::end(shells),
    [&](const char shell) {
      return shellFullOrEmpty(shell);
    }
  );
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
const std::array<double, 110>& bondRadii() {
  static const std::array<double, 110> data {
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
  return data;
}

double bondRadius(const Utils::ElementType elementType) {
  return bondRadii().at(
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


const std::array<ElementInfo, 110>& elementData() {
  static const std::array<ElementInfo, 110> data {{
    // in Angstrom
    {0.00, 0, 0, 0, 0}, // none
    {1.09, 1, 0, 0, 0}, // H
    {1.40, 2, 0, 0, 0}, // He
    {1.82, 1, 0, 0, 0}, // Li
    {1.53, 2, 0, 0, 0}, // Be
    {1.92, 2, 1, 0, 0}, // B
    {1.70, 2, 2, 0, 0}, // C
    {1.55, 2, 3, 0, 0}, // N
    {1.52, 2, 4, 0, 0}, // O
    {1.47, 2, 5, 0, 0}, // F
    {1.54, 2, 6, 0, 0}, // Ne
    {2.27, 1, 0, 0, 0}, // Na
    {1.73, 2, 0, 0, 0}, // Mg
    {1.84, 2, 1, 0, 0}, // Al
    {2.10, 2, 2, 0, 0}, // Si
    {1.80, 2, 3, 0, 0}, // P
    {1.80, 2, 4, 0, 0}, // S
    {1.75, 2, 5, 0, 0}, // Cl
    {1.88, 2, 6, 0, 0}, // Ar
    {2.75, 1, 0, 0, 0}, // K
    {2.31, 2, 0, 0, 0}, // Ca
    {2.15, 2, 0, 1, 0}, // Sc
    {2.11, 2, 0, 2, 0}, // Ti
    {2.07, 2, 0, 3, 0}, // V
    {2.06, 1, 0, 5, 0}, // Cr, Madelung: 2, 0, 4, 0
    {2.05, 2, 0, 5, 0}, // Mn
    {2.04, 2, 0, 6, 0}, // Fe
    {2.00, 2, 0, 7, 0}, // Co
    {1.97, 2, 0, 8, 0}, // Ni
    {1.96, 1, 0, 10, 0}, // Cu, Madelung: 2, 0, 9, 0
    {2.01, 2, 0, 10, 0}, // Zn
    {1.87, 2, 1, 10, 0}, // Ga
    {2.11, 2, 2, 10, 0}, // Ge
    {1.85, 2, 3, 10, 0}, // As
    {1.90, 2, 4, 10, 0}, // Se
    {1.85, 2, 5, 10, 0}, // Br
    {2.02, 2, 6, 10, 0}, // Kr
    {3.03, 1, 0, 0, 0}, // Rb
    {2.49, 2, 0, 0, 0}, // Sr
    {2.32, 2, 0, 1, 0}, // Y
    {2.23, 2, 0, 2, 0}, // Zr
    {2.18, 1, 0, 4, 0}, // Nb, Madelung: 2, 0, 3, 0
    {2.17, 1, 0, 5, 0}, // Mo, Madelung: 2, 0, 4, 0
    {2.16, 2, 0, 5, 0}, // Tc
    {2.13, 1, 0, 7, 0}, // Ru, Madelung: 2, 0, 6, 0
    {2.10, 1, 0, 8, 0}, // Rh, Madelung: 2, 0, 7, 0
    {2.10, 0, 0, 10, 0}, // Pd, Madelung: 2, 0, 8, 0
    {2.11, 1, 0, 10, 0}, // Ag, Madelung: 2, 0, 9, 0
    {2.18, 2, 0, 10, 0}, // Cd
    {1.93, 2, 1, 10, 0}, // In
    {2.17, 2, 2, 10, 0}, // Sn
    {2.06, 2, 3, 10, 0}, // Sb
    {2.06, 2, 4, 10, 0}, // Te
    {1.98, 2, 5, 10, 0}, // I
    {2.16, 2, 6, 10, 0}, // Xe
    {3.43, 1, 0, 0, 0}, // Cs
    {2.68, 2, 0, 0, 0}, // Ba
    //  La: Madelung: [Xe] 6s² 4f¹ -> 2, 0, 0, 1, found [Xe] 6s² 5d¹
    {2.43, 2, 0, 1, 0}, // La
    //  Ce: Madelung: [Xe] 6s² 4f² -> 2, 0, 0, 2, found [Xe] 6s² 4f¹ 5d¹
    {2.42, 2, 0, 1, 1}, // Ce
    {2.40, 2, 0, 0, 3}, // Pr
    {2.39, 2, 0, 0, 4}, // Nd
    {2.38, 2, 0, 0, 5}, // Pm
    {2.36, 2, 0, 0, 6}, // Sm
    {2.35, 2, 0, 0, 7}, // Eu
    //  Gd: Madelung: [Xe] 6s² 4f⁸ -> 2, 0, 0, 8, found [Xe] 6s² 4f⁷ 5d¹
    {2.34, 2, 0, 1, 7}, // Gd
    {2.33, 2, 0, 0, 9}, // Tb
    {2.31, 2, 0, 0, 10}, // Dy
    {2.30, 2, 0, 0, 11}, // Ho
    {2.29, 2, 0, 0, 12}, // Er
    {2.27, 2, 0, 0, 13}, // Tm
    {2.26, 2, 0, 0, 14}, // Yb
    {2.24, 2, 0, 1, 14}, // Lu
    {2.23, 2, 0, 2, 14}, // Hf
    {2.22, 2, 0, 3, 14}, // Ta
    {2.18, 2, 0, 4, 14}, // W
    {2.16, 2, 0, 5, 14}, // Re
    {2.16, 2, 0, 6, 14}, // Os
    {2.13, 2, 0, 7, 14}, // Ir
    //  Pt: Madelung: [Xe] 6s² (4f¹⁴) 5d⁸ -> 2, 0, 8, 14, found [Xe] 6s¹ (4f¹⁴) 5d⁹
    {2.13, 1, 0, 9, 14}, // Pt
    //  Au: Madelung: [Xe] 6s² (4f¹⁴) 5d⁹ -> 2, 0, 9, 14, found [Xe] 6s¹ (4f¹⁴) 5d¹⁰
    {2.14, 1, 0, 10, 14}, // Au
    {2.23, 2, 0, 10, 14}, // Hg
    {1.96, 2, 1, 10, 14}, // Tl
    {2.02, 2, 2, 10, 14}, // Pb
    {2.07, 2, 3, 10, 14}, // Bi
    {1.97, 2, 4, 10, 14}, // Po
    {2.02, 2, 5, 10, 14}, // At
    {2.20, 2, 6, 10, 14}, // Rn
    {3.48, 1, 0, 0, 0}, // Fr
    {2.83, 2, 0, 0, 0}, // Ra
    //  Ac: Madelung: [Rn] 7s² 5f¹ -> 2, 0, 0, 1, found [Rn] 7s² 6d¹
    {2.47, 2, 0, 1, 0}, // Ac
    //  Th: Madelung: [Rn] 7s² 5f² -> 2, 0, 0, 2, found [Rn] 7s² 6d²
    {2.45, 2, 0, 2, 0}, // Th
    //  Pa: Madelung: [Rn] 7s² 5f³ -> 2, 0, 0, 3, found [Rn] 7s² 5f² 6d¹
    {2.43, 2, 0, 1, 2}, // Pa
    //   U: Madelung: [Rn] 7s² 5f⁴ -> 2, 0, 0, 4, found [Rn] 7s² 5f³ 6d¹
    {2.41, 2, 0, 3, 1}, // U,
    //  Np: Madelung: [Rn] 7s² 5f⁵ -> 2, 0, 0, 5, found [Rn] 7s² 5f⁴ 6d¹
    {2.39, 2, 0, 4, 1}, // Np
    {2.43, 2, 0, 0, 6}, // Pu
    {2.44, 2, 0, 0, 7}, // Am
    //  Cm: Madelung: [Rn] 7s² 5f⁸ -> 2, 0, 0, 8, found [Rn] 7s² 5f⁷ 6d¹
    {2.45, 2, 0, 7, 1}, // Cm
    {2.44, 2, 0, 0, 9}, // Bk
    {2.45, 2, 0, 0, 10}, // Cf
    {2.45, 2, 0, 0, 11}, // Es
    {2.45, 2, 0, 0, 12}, // Fm
    {2.46, 2, 0, 0, 13}, // Md
    {2.46, 2, 0, 0, 14}, // No
    {2.46, 2, 0, 1, 14}, // Lr
    {0.0, 2, 0, 2, 14}, // Rf
    {0.0, 2, 0, 3, 14}, // Db
    {0.0, 2, 0, 4, 14}, // Sg
    {0.0, 2, 0, 5, 14}, // Bh
    {0.0, 2, 0, 6, 14}, // Hs
    {0.0, 2, 0, 7, 14} // Mt
  }};
  return data;
}

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
    return elementData().at(
      Utils::ElementInfo::Z(elementType)
    ).valenceElectrons({'s', 'p'});
  }

  return {};
}

unsigned dElectronCount(const Utils::ElementType elementType) {
  if(isMainGroupElement(elementType)) {
    return 0;
  }

  return elementData().at(
    Utils::ElementInfo::Z(elementType)
  ).valenceElectrons('d');
}

double vdwRadius(const Utils::ElementType elementType) {
  return elementData().at(
    Utils::ElementInfo::Z(elementType)
  ).vdwRadius();
}

} // namespace AtomInfo
} // namespace Molassembler
} // namespace Scine
