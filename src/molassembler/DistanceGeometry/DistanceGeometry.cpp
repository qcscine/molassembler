/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/DistanceGeometry/DistanceGeometry.h"

namespace Scine {
namespace molassembler {
namespace distance_geometry {

ChiralConstraint::ChiralConstraint(
  SiteSequence passSites,
  const double passLower,
  const double passUpper
) : sites(std::move(passSites)),
    lower(passLower),
    upper(passUpper)
{
  // Must be <= because flat targets have lower = upper = 0
  if(lower > upper) {
    throw std::runtime_error("Bad construction of chiral constraint");
  }
}

DihedralConstraint::DihedralConstraint(
  SiteSequence passSites,
  const double passLower,
  const double passUpper
) : sites(std::move(passSites)),
    lower(passLower),
    upper(passUpper)
{
  if(lower > upper) {
    throw std::runtime_error("Bad construction of chiral constraint");
  }
}

}
}
}
