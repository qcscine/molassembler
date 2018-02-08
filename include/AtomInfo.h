#ifndef INCLUDE_ATOM_INFO_H
#define INCLUDE_ATOM_INFO_H

#include <map>
#include <vector>
#include <boost/optional.hpp>
#include "Delib/ElementTypes.h" 

/*! @file
 * 
 * A set of electron-counting helper functions. We keep a dataset of s-p-d-f
 * valence electron counts for all elements. These are required in e.g. VSEPR
 * geometry determinations.
 */

namespace MoleculeManip {

//! Information on particular element types
namespace AtomInfo {

/*! @file
 * 
 * Stores information about an element.
 */
class ElementInfo {
private:
  unsigned _valenceElectrons[4];

  unsigned _valenceElectronCount(const char& shell) const {
    switch(shell) {
      case 's':
        return _valenceElectrons[0];
      case 'p':
        return _valenceElectrons[1];
      case 'd':
        return _valenceElectrons[2];
      case 'f':
        return _valenceElectrons[3];
      default:
        return 0;
    }
  }

  unsigned _maxOccupancy(const char& shell) const {
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

public:
  double vdwRadius;

  ElementInfo(
    const double& _vdwRadius,
    const unsigned& sValenceElectrons = 0,
    const unsigned& pValenceElectrons = 0,
    const unsigned& dValenceElectrons = 0,
    const unsigned& fValenceElectrons = 0
  ) : 
    _valenceElectrons {
      sValenceElectrons,
      pValenceElectrons,
      dValenceElectrons,
      fValenceElectrons
    },
    vdwRadius (_vdwRadius)
  {};
  //! Returns the valence electrons for a given shell character (s, p, d, f)
  unsigned valenceElectrons(const char& shell) const {
    return _valenceElectronCount(shell);
  }

  unsigned valenceElectrons(const std::vector<char>& shells) const {
    unsigned sum = 0;
    for(const char& shell : shells) {
      sum += _valenceElectronCount(shell);
    }
    return sum;
  }

  bool shellFullOrEmpty(const char& shell) const {
    return (
      _valenceElectronCount(shell) == 0
      || _valenceElectronCount(shell) == _maxOccupancy(shell)
    );
  }

  bool shellsFullOrEmpty(const std::vector<char>& shells) const {
    for(const char& shell : shells) {
      if(!shellFullOrEmpty(shell)) {
        return false;
      }
    }

    return true;
  }

  //! Returns the total valence electrons
  unsigned valenceElectrons() const {
    unsigned sum = 0;
    for(unsigned i = 0; i < 4; i++) {
      sum += _valenceElectrons[i];
    }
    return sum;
  }
};

/*!
 * Bond radii for each element from the original UFF paper:
 *
 * FIX CITATION
 * Rappé, Goddard et al. UFF, a full periodic table force field for ...
 */
extern const std::array<double, 110> bondRadii;

double bondRadius(const Delib::ElementType& elementType);

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

bool isMainGroupElement(const Delib::ElementType& elementType); 

/*! 
 * Returns a count of valence electrons if the specified element type is a main
 * group element. Otherwise, returns boost::none.
 */
boost::optional<unsigned> mainGroupVE(const Delib::ElementType& elementType);

unsigned dElectronCount(const Delib::ElementType& elementType);

//! Accessor function to fetch the vdw radius directly from elementData
double vdwRadius(const Delib::ElementType& elementType);

} // namespace AtomInfo

} // namespace MoleculeManip

#endif
