#ifndef INCLUDE_GEOMETRY_MODELING_MODEL_VSEPR_H
#define INCLUDE_GEOMETRY_MODELING_MODEL_VSEPR_H

#include "LocalGeometryModel.h"
#include <sstream>

#include "AtomInfo.h"

/*! @file
 *
 * Contains the definition of the VSEPR local geometry modelling scheme.
 */

namespace LocalGeometry {

struct VSEPR : Model {
  static Symmetry::Name determineGeometry(
    const Delib::ElementType& centerAtomType,
    const unsigned& nSites,
    const std::vector<LigandType>& ligands,
    const int& formalCharge
  );
};

} // namespace LocalGeometry

#endif
