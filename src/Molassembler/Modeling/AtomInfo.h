/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Element type information classes
 *
 * A set of electron-counting helper functions. We keep a dataset of s-p-d-f
 * valence electron counts for all elements. These are required in e.g. VSEPR
 * geometry determinations.
 */

#ifndef INCLUDE_MOLASSEMBLER_ATOM_INFO_H
#define INCLUDE_MOLASSEMBLER_ATOM_INFO_H

#include "boost/optional/optional_fwd.hpp"
#include "Utils/Geometry/ElementTypes.h"

#include <array>
#include <map>
#include <vector>

namespace Scine {
namespace Molassembler {

//! Information on particular element types
namespace AtomInfo {

/*!
 * @brief Stores information about an element.
 */
class ElementInfo {
public:
/* Special member functions */
  ElementInfo(
    double passVdwRadius,
    unsigned sValenceElectrons = 0,
    unsigned pValenceElectrons = 0,
    unsigned dValenceElectrons = 0,
    unsigned fValenceElectrons = 0
  );

/* Static functions */
  static unsigned maxOccupancy(char shell);

  //! Returns the valence electrons for a given shell character (s, p, d, f)
  unsigned valenceElectrons(char shell) const;
  unsigned valenceElectrons(const std::vector<char>& shells) const;

  bool shellFullOrEmpty(char shell) const;
  bool shellsFullOrEmpty(const std::vector<char>& shells) const;

  //! Returns the total valence electrons
  unsigned valenceElectrons() const;

  double vdwRadius() const;

private:
  std::array<unsigned, 4> valenceElectrons_;
  double vdwRadius_;
};

/*!
 * Bond radii for each element from the original UFF paper:
 *
 * Rappé, Anthony K., et al. "UFF, a full periodic table force field for
 * molecular mechanics and molecular dynamics simulations." Journal of the
 * American chemical society 114.25 (1992): 10024-10035.
 */
const std::array<double, 110>& bondRadii();

double bondRadius(Utils::ElementType elementType);

/*!
 * ElementData instances for each element type. This is populated with the
 * following data:
 *
 * - VdW radius, from online CRC Handbook of Chemistry and Physics, Nov. 16,
 *   97. ed.
 * - Electron occupation of orbital types s, p, d, f above largest noble core,
 *   populated using Madelung's rule and corrected using Meek, Allen:
 *   Configuration Irregularities, Chem Phys Lett 362, 5-6, August 2002
 *   http://dx.doi.org/10.1016/S0009-2614(02)00919-3
 */
const std::array<ElementInfo, 110>& elementData();

bool isMainGroupElement(Utils::ElementType elementType);

/*!
 * Returns a count of valence electrons if the specified element type is a main
 * group element. Otherwise, returns boost::none.
 */
boost::optional<unsigned> mainGroupVE(Utils::ElementType elementType);

unsigned dElectronCount(Utils::ElementType elementType);

//! Accessor function to fetch the vdw radius directly from elementData
double vdwRadius(Utils::ElementType elementType);

} // namespace AtomInfo
} // namespace Molassembler
} // namespace Scine

#endif
