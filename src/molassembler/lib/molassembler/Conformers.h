/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Interface for the generation of new conformations of Molecules
 */

#ifndef INCLUDE_MOLASSEMBLER_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_CONFORMER_GENERATION_H

#include "Utils/Typenames.h"
#include "molassembler/Types.h"
#include "boost_outcome/outcome.hpp"
#include <vector>

namespace Scine {

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

// Forward-declarations
class Molecule;

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

/**
 * @brief A configuration object for distance geometry runs with sane defaults
 */
struct Configuration {
  /**
   * @brief Choose for how many atoms to re-smooth the distance bounds after
   *   a distance choice
   */
  Partiality partiality {Partiality::FourAtom};

  /**
   * @brief Limit the maximum number of refinement steps
   *
   * The default value is typically enough for medium-sized systems, but may
   * need to be incremented for large systems.
   */
  unsigned refinementStepLimit {10000};

  /**
   * @brief Sets the gradient at which a refinement is considered complete
   *
   * The default value is fairly tight, and can be loosened if faster results
   * are desired and looser local symmetries are tolerable.
   */
  double refinementGradientTarget {1e-5};

  /**
   * @brief Sets the maximum allowed ratio of failures / (# desired conformers)
   *
   * The default value is loose, and allows many failures. It can be tightened
   * (towards lower values) to progressively loosen attempted spatial modelling
   * faster and admit defeat quicker.
   */
  double failureRatio {2};

  /**
   * @brief Set fixed positions for a subset of atoms
   *
   * By default does not set any fixed positions.
   *
   * @pre Any fixed atom must have zero, one or all ligand sites fully
   *   fixed. No individual ligand sites may be partially fixed (i.e. the atoms
   *   constituting a haptic ligand binding site must be either completely
   *   unfixed or fixed and may not be mixed).
   *
   * @note Remember Scine::Utils::Positions are in bohr length units!
   */
  std::vector<
    std::pair<AtomIndex, Scine::Utils::Position>
  > fixedPositions;
};

} // namespace DistanceGeometry

/*!
 * @brief Generate multiple sets of positional data for a Molecule
 *
 * In the case of a molecule that does not have unassigned stereopermutators,
 * this is akin to generating a conformational ensemble. If there are
 * unassigned stereopermutators, these are assigned at random (consistent with
 * relative statistical occurrences of stereopermutations) for each structure.
 *
 * @param molecule The molecule for which to generate sets of three-dimensional
 *   positions. This molecule may not contain stereopermutators with zero
 *   assignments.
 * @param numStructures The number of desired structures to generate
 * @param configuration The configuration object to control Distance Geometry
 *   in detail. The defaults are usually fine.
 *
 * @pre @p molecule may not contain stereopermutators with zero assignments.
 * @pre @p configuration's preconditions must be met
 *
 * @throws std::runtime_error if any preconditions are unmet
 *
 * @returns A result type which may or may not contain a vector of
 *   PositionCollections (in Bohr length units). The result type is much like
 *   an optional, except that in the error case it carries data about the error
 *   in order to help diagnose possible mistakes made in the molecular graph
 *   specification.
 */
outcome::result<
  std::vector<Scine::Utils::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  unsigned numStructures,
  const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
);

/*! Generate a 3D structure of a Molecule
 *
 * @param molecule The molecule for which to generate three-dimensional
 *   positions. This molecule may not contain stereopermutators with zero
 *   assignments.
 * @param configuration The configuration object to control Distance Geometry
 *   in detail. The defaults are usually fine.
 *
 * @returns A result type which may or may not contain a PositionCollection (in
 *   Bohr length units). The result type is much like an optional, except that
 *   in the error case it carries data about the error in order to help
 *   diagnose possible mistakes made in the molecular graph specification.
 */
outcome::result<Scine::Utils::PositionCollection> generateConformation(
  const Molecule& molecule,
  const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
);

} // namespace molassembler

} // namespace Scine

#endif
