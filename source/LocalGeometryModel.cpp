#include "LocalGeometryModel.h"

namespace LocalGeometry {

const std::map<BondType, double> bondWeights {
  {BondType::Single, 1.0},
  {BondType::Double, 2.0},
  {BondType::Triple, 3.0},
  {BondType::Quadruple, 4.0},
  {BondType::Quintuple, 5.0},
  {BondType::Sextuple, 6.0},
  {BondType::Aromatic, 1.5},
  {BondType::Eta, 0.0} // TODO is this wise? (duplicate in AtomInfo!)
};

} // namespace LocalGeometry
