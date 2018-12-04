// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_ATOM_INFO_H
#define INCLUDE_MOLASSEMBLER_ATOM_INFO_H

#include "boost/optional/optional_fwd.hpp"
#include "Utils/ElementTypes.h"

#include <array>
#include <map>
#include <vector>

/*! @file
 *
 * @brief Element type information classes
 *
 * A set of electron-counting helper functions. We keep a dataset of s-p-d-f
 * valence electron counts for all elements. These are required in e.g. VSEPR
 * geometry determinations.
 */

namespace molassembler {

//! Information on particular element types
namespace AtomInfo {

/*! @file
 *
 * Stores information about an element.
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
  std::array<unsigned, 4> _valenceElectrons;
  double _vdwRadius;
};

/*!
 * Bond radii for each element from the original UFF paper:
 *
 * Rappé, Anthony K., et al. "UFF, a full periodic table force field for
 * molecular mechanics and molecular dynamics simulations." Journal of the
 * American chemical society 114.25 (1992): 10024-10035.
 */
extern const std::array<double, 110> bondRadii;

double bondRadius(Scine::Utils::ElementType elementType);

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
extern std::array<ElementInfo, 110> elementData;

bool isMainGroupElement(Scine::Utils::ElementType elementType);

/*!
 * Returns a count of valence electrons if the specified element type is a main
 * group element. Otherwise, returns boost::none.
 */
boost::optional<unsigned> mainGroupVE(Scine::Utils::ElementType elementType);

unsigned dElectronCount(Scine::Utils::ElementType elementType);

//! Accessor function to fetch the vdw radius directly from elementData
double vdwRadius(Scine::Utils::ElementType elementType);

} // namespace AtomInfo

} // namespace molassembler

#endif
