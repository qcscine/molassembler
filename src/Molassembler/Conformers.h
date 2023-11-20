/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interface for the generation of new conformations of Molecules
 */

#ifndef INCLUDE_MOLASSEMBLER_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_CONFORMER_GENERATION_H

#include "Molassembler/Types.h"
#include "Utils/Typenames.h"
#include "Molassembler/Detail/Outcome.h"
#include <vector>

namespace Scine {
namespace Molassembler {

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
 * tradeoff effects:
 *
 * - Reduces computational effort in generating a distance matrix, speeding up
 *   conformer generation
 * - Worsens the quality of initial embedded coordinates prior to refinement,
 *   requiring more refinement steps
 */
enum class MASM_EXPORT Partiality {
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
struct MASM_EXPORT Configuration {
  /**
   * @brief Choose for how many atoms to re-smooth the distance bounds after
   *   a distance choice
   *
   * FourAtom partiality is default because incurring more refinement steps
   * costs less than resmoothing the distance bounds matrix after each
   * successive distance choice.
   */
  Partiality partiality {Partiality::FourAtom};

  /**
   * @brief Limit the maximum number of refinement steps
   *
   * The default value is typically enough for medium-sized systems, but may
   * need to be incremented for large systems.
   */
  unsigned refinementStepLimit {10'000};

  /**
   * @brief Sets the gradient at which a refinement is considered complete
   *
   * The default value is fairly tight, and can be loosened if faster results
   * are desired and looser local shapes are tolerable.
   */
  double refinementGradientTarget {1e-5};

  /**
   * @brief Sets the loosening of the spatial model
   *
   * If a molecule's distortions are not well encompassed in the spatial model,
   * spatial modeling may claim a molecular graph is impossible to represent
   * in three dimensions. Custom loosening of all bounds (distances, angles and
   * dihedrals) may enable molassembler to generate conformers.
   *
   * The default value does not add in any loosening. It is suggested to stay
   * below 2.0.
   */
  double spatialModelLoosening {1.0};

  /**
   * @brief Set fixed positions for a subset of atoms
   *
   * By default does not set any fixed positions.
   *
   * @pre Any fixed atom must have zero, one or all binding sites fully
   *   fixed. No individual ligand sites may be partially fixed (i.e. the atoms
   *   constituting a haptic ligand binding site must be either completely
   *   unfixed or fixed and may not be mixed).
   *
   * @note Remember Utils::Positions are in bohr length units!
   */
  std::vector<
    std::pair<AtomIndex, Utils::Position>
  > fixedPositions;
};

} // namespace DistanceGeometry

/*! @brief Generate multiple sets of positional data for a Molecule
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
 * @pre @p molecule may not contain stereopermutators with zero assignments as
 *   this means that the molecule is not representable in three dimensions.
 * @pre @p configuration's preconditions must be met
 *
 * @complexity{Roughly @math{O(C \cdot N^3)} where @math{C} is the number of
 * conformers and @math{N} is the number of atoms in @p molecule}
 *
 * @parblock @note The Distance Geometry procedure can fail stochastically
 * without fault in the input. See the documentation of DgError for detailed
 * description of error return codes and how to deal with them.
 * @endparblock
 *
 * @parblock @note This function is parallelized. Use the OMP_NUM_THREADS
 * environment variable to control the number of threads used. Results are
 * sequenced and reproducible.
 * @endparblock
 *
 * @parblock @note This function advances the state of the global PRNG.
 * @endparblock
 *
 * @returns A result type which may or may not contain a vector of
 * PositionCollections (in Bohr length units). The result type is much like an
 * optional, except that in the error case it carries data about the error in
 * order to help diagnose possible mistakes made in the molecular graph
 * specification.
 *
 * @code{.cpp}
 * auto ensemble = generateRandomEnsemble(mol, 10);
 * unsigned conformerIndex = 0;
 * for(auto& conformerResult : ensemble) {
 *   if(conformerResult) {
 *     IO::write(
 *       std::to_string(conformerIndex) + ".mol",
 *       mol,
 *       conformerResult.value()
 *     );
 *   } else {
 *     std::cout << "Conformer " << conformerIndex << " failed: "
 *       << conformerResult.error().message() << "\n";
 *   }
 *   ++conformerIndex;
 * }
 * @endcode
 */
MASM_EXPORT std::vector<
  Result<Utils::PositionCollection>
> generateRandomEnsemble(
  const Molecule& molecule,
  unsigned numStructures,
  const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
);

/*! @brief Generate multiple sets of positional data for a Molecule
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
 * @param seed A number to seed the pseudo-random number generator used in
 *   conformer generation with
 * @param configuration The configuration object to control Distance Geometry
 *   in detail. The defaults are usually fine.
 *
 * @pre @p molecule may not contain stereopermutators with zero assignments as
 *   this means that the molecule is not representable in three dimensions.
 * @pre @p configuration's preconditions must be met
 *
 * @complexity{Roughly @math{O(C \cdot N^3)} where @math{C} is the number of
 * conformers and @math{N} is the number of atoms in @p molecule}
 *
 * @parblock @note The Distance Geometry procedure can fail stochastically
 * without fault in the input. See the documentation of DgError for detailed
 * description of error return codes and how to deal with them.
 * @endparblock
 *
 * @parblock @note This function is parallelized. Use the OMP_NUM_THREADS
 * environment variable to control the number of threads used. Results are
 * sequenced and reproducible.
 * @endparblock
 *
 * @returns A result type which may or may not contain a vector of
 * PositionCollections (in Bohr length units). The result type is much like an
 * optional, except that in the error case it carries data about the error in
 * order to help diagnose possible mistakes made in the molecular graph
 * specification.
 */
MASM_EXPORT std::vector<
  Result<Utils::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  unsigned numStructures,
  unsigned seed,
  const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
);

/*! @brief Generate a 3D structure of a Molecule
 *
 * @param molecule The molecule for which to generate three-dimensional
 *   positions. This molecule may not contain stereopermutators with zero
 *   assignments.
 * @param configuration The configuration object to control Distance Geometry
 *   in detail. The defaults are usually fine.
 *
 * @complexity{Roughly @math{O(N^3)}}
 *
 * @parblock @note This function advances the state of the global PRNG.
 * @endparblock
 *
 * @see generateEnsemble
 *
 * @returns A result type which may or may not contain a PositionCollection (in
 *   Bohr length units). The result type is much like an optional, except that
 *   in the error case it carries data about the error in order to help
 *   diagnose possible mistakes made in the molecular graph specification.
 */
MASM_EXPORT Result<Utils::PositionCollection> generateRandomConformation(
  const Molecule& molecule,
  const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
);

/*! @brief Generate a 3D structure of a Molecule
 *
 * @param molecule The molecule for which to generate three-dimensional
 *   positions. This molecule may not contain stereopermutators with zero
 *   assignments.
 * @param seed A number to seed the pseudo-random number generator used in
 *   conformer generation with
 * @param configuration The configuration object to control Distance Geometry
 *   in detail. The defaults are usually fine.
 *
 * @complexity{Roughly @math{O(N^3)}}
 *
 * @see generateEnsemble
 *
 * @returns A result type which may or may not contain a PositionCollection (in
 *   Bohr length units). The result type is much like an optional, except that
 *   in the error case it carries data about the error in order to help
 *   diagnose possible mistakes made in the molecular graph specification.
 */
MASM_EXPORT Result<Utils::PositionCollection> generateConformation(
  const Molecule& molecule,
  unsigned seed,
  const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
);

} // namespace Molassembler
} // namespace Scine

#endif
