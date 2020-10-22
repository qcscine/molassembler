/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/DistanceGeometry/DistanceGeometry.h"

namespace Scine {
namespace Molassembler {
namespace DistanceGeometry {

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

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine
