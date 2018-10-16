// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_H

#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Types.h"

#include <tuple>
#include <vector>
#include <cassert>
#include <array>

/*! @file
 *
 * @brief Data struct to store chiral constraints for DG
 *
 * Contains some central data class declarations and type definitions for the
 * entire Distance Geometry scheme.
 */

namespace molassembler {

//! Distance geometry-related classes and functions
namespace DistanceGeometry {

struct ChiralityConstraint {
  using AtomListType = std::vector<AtomIndex>;
  using LigandSequence = std::array<AtomListType, 4>;

  LigandSequence sites;
  double lower, upper;

  ChiralityConstraint(
    LigandSequence passSites,
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
 * @brief Limit triangle inequality bounds smoothing to a subset of all atoms
 *
 * Usually, after choosing a single conformer's atom-pairwise distances from
 * between the distance bounds that are generated from the spatial model, all
 * other distance bounds are re-smoothed using the triangle inequality. However,
 * the reduction in distance bounds slack all over the matrix diminishes with
 * each successive chosen distance, yielding only diminishing returns.
 *
 * With this option, you can choose to stop re-smoothing the entire matrix after
 * a limited number of one-to-all distance choices. This has the following
 * effects:
 *
 * - Reduces computational effort in generating a distance matrix, speeding up
 *   the conformer generation significantly.
 * - Worsens the quality of embedded coordinates prior to refinement.
 */
enum class Partiality {
  /*!
   * @brief Perform smoothing for four one-to-all distance choices
   *
   * In principle, if the distances from four atoms to all others are known,
   * the overall conformation is fully determined. This is unfortunately not
   * realized in practice in the DG procedure.
   *
   * This yields the most speedup in distance matrix generation, but also
   * worsens the initial embedded coordinates the most.
   */
  FourAtom,
  /*!
   * @brief Perform smoothing for ten percent of one-to-all distance choices
   *
   * @note Still performs bounds smoothing for at least four one-to-all
   *   distance choices if the number of atoms is less than 40.
   */
  TenPercent,
  /*!
   * @brief Perform smoothing after all distance choices
   *
   * This yields the slowest distance matrix generation, but also the best
   * initial embedded coordinates.
   */
  All
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
