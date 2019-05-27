/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Calculate symmetry position maps
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SYMMETRY_POSITION_MAPS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SYMMETRY_POSITION_MAPS_H

#include "molassembler/RankingInformation.h"

namespace Scine {

namespace stereopermutation {
struct Stereopermutation;
} // namespace stereopermutation

namespace molassembler {

/*!
 * @brief Generates a flat mapping from site indices to symmetry positions
 *
 * Generates a mapping from site indices to symmetry positions according to
 * the ranking character distribution to symmetry positions of a
 * stereopermutation (its characters member) and any defined links between
 * symmetry positions.
 *
 * @code{.cpp}
 * auto mapping = generateSiteToSymmetryPosition(...);
 * unsigned symmetryPositionOfSiteFour = mapping.at(4u);
 * @endcode
 */
std::vector<unsigned> siteToSymmetryPositionMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
);

/*!
 * @brief Generates a flat mapping from symmetry positions to site indices
 *
 * Generates exactly the inverse map to generateSiteToSymmetryPositionMap
 *
 * @code{cpp}
 * auto mapping = generateSymmetryPositionToSiteMap(...);
 * unsigned siteIndexAtSymmetryPositionFive = mapping.at(5u);
 * @endcode
 */
std::vector<unsigned> symmetryPositionToSiteMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
);

} // namespace molassembler
} // namespace Scine

#endif
