/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Data struct to store chiral constraints for DG
 *
 * Contains some central data class declarations and type definitions for the
 * entire Distance Geometry scheme.
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_H

#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Types.h"

#include <tuple>
#include <vector>
#include <cassert>
#include <array>

namespace Scine {

namespace molassembler {

//! Distance geometry-related classes and functions
namespace distance_geometry {

/**
 * @brief Data struct representing a chiral constraint
 *
 * Contains four atom sets and lower and upper bounds on the signed volume
 * spanned by their average spatial positions.
 */
struct ChiralConstraint {
  using AtomListType = std::vector<AtomIndex>;
  using SiteSequence = std::array<AtomListType, 4>;

  //! Site definition sequence (odd permutations invert sign)
  SiteSequence sites;
  //! Lower bound on signed volume
  double lower;
  //! Upper bound on signed volume
  double upper;
  /*! @brief Weight of the chiral constraint
   *
   * This allows tuning of the relative importance of chiral constraints.
   */
  double weight = 1.0;

  ChiralConstraint(
    SiteSequence passSites,
    const double passLower,
    const double passUpper
  ) : sites(std::move(passSites)),
      lower(passLower),
      upper(passUpper)
  {
    // Must be <= because flat targets have lower = upper = 0
    assert(lower <= upper);
  }
};

/**
 * @brief Data struct representing a dihedral constraint
 *
 * Contains four atom sets and lower and upper bounds on the dihedral spanned
 * by their average spatial positions.
 */
struct DihedralConstraint {
  using AtomListType = std::vector<AtomIndex>;
  using SiteSequence = std::array<AtomListType, 4>;

  //! Site definition sequence (odd permutations invert sign)
  SiteSequence sites;
  //! Lower bound on signed volume
  double lower;
  //! Upper bound on signed volume
  double upper;

  DihedralConstraint(
    SiteSequence passSites,
    const double passLower,
    const double passUpper
  ) : sites(std::move(passSites)),
      lower(passLower),
      upper(passUpper)
  {
    assert(lower <= upper);
  }
};

} // namespace distance_geometry

} // namespace molassembler

} // namespace Scine

#endif
