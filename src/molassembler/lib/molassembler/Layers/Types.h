#ifndef INCLUDE_MOLASSEMBLER_SHARED_TYPES_H
#define INCLUDE_MOLASSEMBLER_SHARED_TYPES_H

namespace molassembler {

/*!
 * Bond type enumeration. Besides the classic organic single, double and triple
 * bonds, bond orders up to sextuple are explicitly included.
 *
 * Although currently unused, Aromatic and Eta bonds are included in
 * anticipation of their necessity.
 */
enum class BondType : unsigned {
  Single,
  Double,
  Triple,
  Quadruple,
  Quintuple,
  Sextuple,
  Aromatic,
  Eta
};

enum class LengthUnit {
  Bohr,
  Angstrom
};

} // namespace molassembler

#endif
