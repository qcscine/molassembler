/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Calculate shape position maps
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SYMMETRY_POSITION_MAPS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SYMMETRY_POSITION_MAPS_H

#include "molassembler/RankingInformation.h"

namespace Scine {

namespace stereopermutation {
class Stereopermutation;
} // namespace stereopermutation

namespace molassembler {

/*! @brief Generates a flat mapping from site indices to shape positions
 *
 * Generates a mapping from site indices to shape positions according to
 * the ranking character distribution to shape positions of a
 * stereopermutation (its characters member) and any defined links between
 * shape positions.
 *
 * @code{.cpp}
 * auto mapping = generateSiteToSymmetryPosition(...);
 * unsigned symmetryPositionOfSiteFour = mapping.at(4u);
 * @endcode
 *
 * @complexity{@math{\Theta(N)}}
 */
std::vector<unsigned> siteToSymmetryPositionMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
);

/*! @brief Generates a flat mapping from shape positions to site indices
 *
 * Generates exactly the inverse map to generateSiteToSymmetryPositionMap
 *
 * @code{cpp}
 * auto mapping = generateSymmetryPositionToSiteMap(...);
 * unsigned siteIndexAtSymmetryPositionFive = mapping.at(5u);
 * @endcode
 *
 * @complexity{@math{\Theta(N)}}
 */
std::vector<unsigned> symmetryPositionToSiteMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
);

} // namespace molassembler
} // namespace Scine

#endif
