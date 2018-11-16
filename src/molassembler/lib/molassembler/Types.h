// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_SHARED_TYPES_H
#define INCLUDE_MOLASSEMBLER_SHARED_TYPES_H

#include <cstddef>

/*!@file
 *
 * @brief Defines basic types widely shared across the project.
 */

namespace molassembler {

/*!
 * @brief Discrete bond type numeration
 *
 * Bond type enumeration. Besides the classic organic single, double and triple
 * bonds, bond orders up to sextuple are explicitly included.
 */
enum class BondType : unsigned {
  Single,
  Double,
  Triple,
  Quadruple,
  Quintuple,
  Sextuple,
  Eta
};

enum class LengthUnit {
  Bohr,
  Angstrom
};

//! Unsigned integer atom index type. Used to refer to particular atoms.
using AtomIndex = std::size_t;

//! Type used to refer to particular bonds. Orders first < second.
struct BondIndex {
  AtomIndex first, second;

  BondIndex();
  BondIndex(AtomIndex a, AtomIndex b) noexcept;

  bool operator < (const BondIndex& other) const;
  bool operator == (const BondIndex& other) const;
};

std::size_t hash_value(const BondIndex& bond);

//! Descriptive name for dlib indices
using dlibIndexType = long;

/*!
 * @brief For bitmasks grouping components of immediate atom environments
 *
 * Differing strictnesses of comparisons may be desirable for various
 * purposes, hence a modular comparison function is provided.
 */
enum class AtomEnvironmentComponents : unsigned {
  ElementTypes,
  BondOrders,
  Symmetries,
  Stereopermutations // Symmetries must be set in conjunction with this
};

namespace DistanceGeometry {

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
