/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Calculate shape position maps
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SYMMETRY_POSITION_MAPS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SYMMETRY_POSITION_MAPS_H

#include "molassembler/RankingInformation.h"

#include "temple/StrongIndexMap.h"
#include "shapes/Data.h"

namespace Scine {

namespace stereopermutation {
class Stereopermutation;
} // namespace stereopermutation

namespace molassembler {

using SiteToShapeVertexMap = temple::StrongIndexFlatMap<SiteIndex, shapes::Vertex>;

/*! @brief Generates a flat mapping from site indices to shape vertices
 *
 * Generates a mapping from site indices to shape vertices according to
 * the ranking character distribution to shape vertex of a stereopermutation
 * (its characters member) and any defined links between shape positions.
 *
 * @code{.cpp}
 * auto mapping = siteToShapeVertexMap(...);
 * unsigned symmetryPositionOfSiteFour = mapping.at(4u);
 * @endcode
 *
 * @complexity{@math{\Theta(N)}}
 */
SiteToShapeVertexMap siteToShapeVertexMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
);

/*! @brief Generates a flat mapping from shape vertices to site indices
 *
 * Generates exactly the inverse map to generateSiteToSymmetryPositionMap
 *
 * @code{cpp}
 * auto mapping = shapeVertexToSiteIndexMap(...);
 * unsigned siteIndexAtSymmetryPositionFive = mapping.at(5u);
 * @endcode
 *
 * @complexity{@math{\Theta(N)}}
 */
temple::StrongIndexFlatMap<shapes::Vertex, SiteIndex> shapeVertexToSiteIndexMap(
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites
);

stereopermutation::Stereopermutation stereopermutationFromSiteToShapeVertexMap(
  const SiteToShapeVertexMap& siteToShapeVertexMap,
  const std::vector<LinkInformation>& links,
  const RankingInformation::RankedSitesType& canonicalSites
);

} // namespace molassembler
} // namespace Scine

#endif
