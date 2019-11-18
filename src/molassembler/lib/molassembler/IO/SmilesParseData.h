/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Define data accrued in smiles parser for consumption in the
 *   accompanying molecule builder
 */
#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_DATA_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_PARSER_DATA_H

#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include "Utils/Geometry/ElementInfo.h"
#include "shapes/Shapes.h"
#include "boost/optional.hpp"
#include "molassembler/Types.h"

namespace Scine {
namespace molassembler {
namespace IO {

struct ElementData {
  unsigned Z = 0;
  bool aromatic = false;

  ElementData() = default;
  explicit ElementData(Utils::ElementType e) : Z(Utils::ElementInfo::Z(e)) {}

  static inline ElementData aromaticElement(Utils::ElementType e) {
    ElementData d(e);
    d.aromatic = true;
    return d;
  }
};

struct ChiralData {
  Shapes::Shape shape;
  unsigned chiralIndex;
};

struct AtomData {
  unsigned A = 0;
  ElementData partialElement;
  boost::optional<ChiralData> chiralOptional;
  boost::optional<unsigned> hCount;
  boost::optional<int> chargeOptional;
  bool atomBracket = false;

  inline Utils::ElementType getElement() const {
    if(partialElement.Z == 0) {
      return Utils::ElementType::none;
    }

    if(A == 0) {
      return Utils::ElementInfo::element(partialElement.Z);
    }

    return Utils::ElementInfo::isotope(partialElement.Z, A);
  }
};

struct BondData {
  enum class StereoMarker {Forward, Backward};

  boost::optional<BondType> type;
  boost::optional<StereoMarker> ezStereo;
  boost::optional<unsigned> ringNumber;
};

} // namespace IO
} // namespace molassembler
} // namespace Scine

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::ElementData,
  (unsigned, Z),
  (bool, aromatic)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::ChiralData,
  (Scine::Shapes::Shape, shape),
  (unsigned, chiralIndex)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::AtomData,
  (unsigned, A),
  (Scine::molassembler::IO::ElementData, partialElement),
  (boost::optional<Scine::molassembler::IO::ChiralData>, chiralOptional),
  (boost::optional<unsigned>, hCount),
  (boost::optional<int>, chargeOptional),
  (bool, atomBracket)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::BondData,
  (boost::optional<Scine::molassembler::BondType>, type),
  (boost::optional<Scine::molassembler::IO::BondData::StereoMarker>, ezStereo)
  (boost::optional<unsigned>, ringNumber)
)

#endif
