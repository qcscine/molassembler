#ifndef INCLUDE_ATOM_INFO_H
#define INCLUDE_ATOM_INTO_H

#include <map>
#include <vector>
#include <boost/optional.hpp>
#include "ElementTypes.h"

namespace MoleculeManip {

namespace AtomInfo {

class ElementInfo {
private:
  unsigned _valenceElectrons[4];
  unsigned _valenceElectronCount(const char& shell) {
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
  unsigned valenceElectrons(const char& shell) {
    return _valenceElectronCount(shell);
  }
  unsigned valenceElectrons(const std::vector<char>& shells) {
    unsigned sum = 0;
    for(const char& shell : shells) {
      sum += _valenceElectronCount(shell);
    }
    return sum;
  }
  //! Returns the total valence electrons
  unsigned valenceElectrons() {
    unsigned sum;
    for(unsigned i = 0; i < 4; i++) {
      sum += _valenceElectrons[i];
    }
    return sum;
  }
};

/* From the original UFF paper
 * Rappé, Goddard et al. UFF, a full periodic table force field for ...
 */
extern const std::map<Delib::ElementType, double> bondRadii;

/* ElementData is populated with the following data:
 * - VdW radius, from online CRC Handbook of Chemistry and Physics, Nov. 16,
 *   97. ed.
 * - Electron occupation of orbital types s, p, d, f above largest noble core, 
 *   populated using Madelung's rule and corrected using Meek, Allen: 
 *   Configuration Irregularities, Chem Phys Lett 362, 5-6, August 2002
 *   http://dx.doi.org/10.1016/S0009-2614(02)00919-3
 */
extern std::map<
  Delib::ElementType,
  ElementInfo
> elementData;

bool isMainGroupElement(const Delib::ElementType& elementType); 

boost::optional<unsigned> mainGroupVE(const Delib::ElementType& elementType);

unsigned dElectronCount(const Delib::ElementType& elementType);

double vdwRadius(const Delib::ElementType& elementType);

} // eo namespace AtomInfo

} // eo namespace MoleculeManip

#endif
