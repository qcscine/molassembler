/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Calculate shape position maps
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SHAPE_VERTEX_MAPS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_SHAPE_VERTEX_MAPS_H

#include "Molassembler/RankingInformation.h"

#include "Molassembler/Temple/StrongIndexPermutation.h"
#include "Molassembler/Shapes/Data.h"

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {
class Stereopermutation;
} // namespace Stereopermutations

using SiteToShapeVertexMap = Temple::StrongIndexPermutation<SiteIndex, Shapes::Vertex>;

/*! @brief Generates a flat mapping from site indices to shape vertices
 *
 * Generates a mapping from site indices to shape vertices according to
 * the ranking character distribution to shape vertex of a stereopermutation
 * (its characters member) and any defined links between shape positions.
 *
 * If the shape position groups for the sites and vertices are provided, this
 * information is used to select the correct site to vertex map if multiple
 * are isomorphic.
 *
 * @code{.cpp}
 * auto mapping = siteToShapeVertexMap(...);
 * Shapes::Vertex shapeVertexOfSiteFour = mapping.at(SiteIndex {4});
 * @endcode
 *
 * @complexity{@math{\Theta(N)}}
 */
SiteToShapeVertexMap siteToShapeVertexMap(
  const Stereopermutations::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<RankingInformation::Link>& siteLinks,
  std::vector<std::vector<SiteIndex>> siteGroups = {},
  std::vector<std::vector<Shapes::Vertex>> vertexGroups = {}
);

/*! @brief Generates a flat mapping from shape vertices to site indices
 *
 * Generates exactly the inverse map to generateSiteToShapeVertexMap
 *
 * @complexity{@math{\Theta(N)}}
 */
Temple::StrongIndexPermutation<Shapes::Vertex, SiteIndex> shapeVertexToSiteIndexMap(
  const Stereopermutations::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const std::vector<RankingInformation::Link>& siteLinks
);

Stereopermutations::Stereopermutation stereopermutationFromSiteToShapeVertexMap(
  const SiteToShapeVertexMap& siteToShapeVertexMap,
  const std::vector<RankingInformation::Link>& links,
  const RankingInformation::RankedSitesType& canonicalSites
);

} // namespace Molassembler
} // namespace Scine

#endif
