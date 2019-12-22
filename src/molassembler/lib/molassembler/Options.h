/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Centralizes the main customization points of the library's behavior.
 */

#ifndef INCLUDE_MOLASSEMBLER_OPTIONS_H
#define INCLUDE_MOLASSEMBLER_OPTIONS_H

#include "shapes/Shapes.h"
#include "molassembler/Prng.h"

#include "Utils/Geometry/ElementTypes.h"
#include "molassembler/AngstromWrapper.h"

namespace Scine {
namespace molassembler {

// Forward-declarations
class OuterGraph;

/*! @brief Randomness source for the entire library
 *
 * This instance supplies the library with randomness. The default seeding
 * behavior is build-type dependent. In debug builds, the seed is fixed. In
 * release builds, the generator is seeded from std::random_device.
 *
 * @note If you wish to get deterministic behavior in release builds, re-seed
 * the generator.
 *
 * @complexity{@math{\Theta(1)}}
 *
 * @warning Do not use this instance in any static object's destructor!
 */
MASM_EXPORT random::Engine& randomnessEngine();

/*! @brief Modeling temperature regime enumerator
 *
 * Specifies for which temperature regime the Molecule is being modeled.
 *
 * This affects
 * - whether nitrogen atoms with a vacant tetrahedron shape and not part of a
 *   small cycle can be stereocenters (nitrogen inversion)
 * - whether trigonal- and pentagonal bipyramids without linked ligands can be
 *   stereocenters (Berry pseudorotation and Bartell mechanism)
 */
enum class MASM_EXPORT TemperatureRegime {
  /*! @brief No thermalization of stereopermutations at atom stereopermutators
   *
   * Meaning no pyramidal inversion, Berry pseudorotation or Bartell mechanisms.
   */
  Low,
  /*! @brief In certain circumstances, stereopermutations are thermalized
   *
   * Only vacant tetrahedron nitrogen atoms in particularly strained cycles (3,
   * 4) can be stereopermutators. Berry pseudorotation and Bartell mechanisms
   * thermalize all stereopermutations of trigonal bipyramid and pentagonal
   * bipyramid shaped centers if none of its substituents are linked.
   */
  High
};

/*! @brief Specifies the effects of graph modifications on chiral centers
 *
 * In case a substituent is added or removed at a stereopermutator, an attempt
 * is made to transfer chiral information into the new geometry. How this
 * attempt is made can be altered.
 *
 * The following graphs are to illustrate which chiral information transfers
 * are possible under which setting.
 *
 * The first is for the situation of ligand gain. There are very few edges and
 * shapes in this graph because in nearly all cases, there are multiple
 * best mappings, and propagation in these cases are enabled only in
 * RandomFromMultipleBest, in which case the graph gets quite dense, so it's
 * not shown. The green edges here are the ligand gain situations propagated
 * under EffortlessAndUnique.
 *
 * \dot
 *  digraph g {
 *    graph [fontname = "Arial", nodesep="0.2", ranksep="0.8"];
 *    node [fontname = "Arial", style = "filled", fillcolor="white"];
 *    edge [fontname = "Arial", penwidth=2, labelfontsize="10"];
 *    rankdir="LR";
 *
 *    cuttetrahedral [label="cut tetrahedral"];
 *    Tshaped [label="T-shaped"];
 *    tetrahedral [label="tetrahedral"];
 *    squareplanar [label="square planar"];
 *    seesaw [label="seesaw"];
 *    squarepyramidal [label="square pyramidal"];
 *    trigonalbipyramidal [label="trigonal bipyramidal"];
 *    octahedral [label="octahedral"];
 *    pentagonalpyramidal [label="pentagonal pyramidal"];
 *    pentagonalbipyramidal [label="pentagonal bipyramidal"];
 *    cuttetrahedral -> tetrahedral [color="forestgreen"];
 *    Tshaped -> squareplanar [color="forestgreen"];
 *    seesaw -> trigonalbipyramidal [color="forestgreen"];
 *    squarepyramidal -> octahedral [color="forestgreen"];
 *    pentagonalpyramidal -> pentagonalbipyramidal [color="forestgreen"];
 *  }
 * \enddot
 *
 * The second is for the situation of ligand loss. The edge label indicates
 * the group of shape ligands from which a ligand is being removed. Green
 * edges represent situtations propagated under EffortlessAndUnique, black
 * edges those under Unique. RandomFromMultipleBest is not shown due to graph
 * density.
 *
 * \dot
 *  digraph g {
 *    graph [fontname = "Arial", nodesep="0.8", ranksep="0.8"];
 *    node [fontname = "Arial", style = "filled", fillcolor="white"];
 *    edge [fontname = "Arial", penwidth=2, labelfontsize="10"];
 *    rankdir="LR";
 *
 *    cuttetrahedral [label="cut tetrahedral"];
 *    Tshaped [label="T-shaped"];
 *    tetrahedral [label="tetrahedral"];
 *    squareplanar [label="square planar"];
 *    seesaw [label="seesaw"];
 *    squarepyramidal [label="square pyramidal"];
 *    trigonalbipyramidal [label="trigonal bipyramidal"];
 *    pentagonalplanar [label="pentagonal planar"];
 *    octahedral [label="octahedral"];
 *    trigonalprismatic [label="trigonal prismatic"];
 *    pentagonalpyramidal [label="pentagonal pyramidal"];
 *    pentagonalbipyramidal [label="pentagonal bipyramidal"];
 *
 *    tetrahedral -> cuttetrahedral [color="forestgreen", label="any"];
 *    squareplanar -> Tshaped [color="forestgreen", label="any"];
 *    seesaw -> cuttetrahedral [color="black", label="axial"];
 *    seesaw -> Tshaped [color="black", label="axial"];
 *    seesaw -> Tshaped [color="forestgreen", label="equat."];
 *    squarepyramidal -> tetrahedral [color="black", label="equat."];
 *    squarepyramidal -> squareplanar [color="black", label="equat."];
 *    squarepyramidal -> squareplanar [color="forestgreen", label="axial"];
 *    squarepyramidal -> seesaw [color="black", label="equat."];
 *    trigonalbipyramidal -> tetrahedral [color="black", label="equat."];
 *    trigonalbipyramidal -> tetrahedral [color="black", label="axial"];
 *    trigonalbipyramidal -> squareplanar [color="black", label="equat."];
 *    trigonalbipyramidal -> seesaw [color="forestgreen", label="equat."];
 *    pentagonalplanar -> squareplanar [color="black", label="any"];
 *    octahedral -> squarepyramidal [color="forestgreen", label="any"];
 *    trigonalprismatic -> pentagonalplanar [color="black", label="any"];
 *    pentagonalpyramidal -> squarepyramidal [color="black", label="equat."];
 *    pentagonalpyramidal -> pentagonalplanar [color="black", label="equat."];
 *    pentagonalpyramidal -> pentagonalplanar [color="forestgreen", label="apical"];
 *    pentagonalbipyramidal -> octahedral [color="black", label="equat."];
 *    pentagonalbipyramidal -> pentagonalpyramidal [color="forestgreen", label="axial"];
 *  }
 *   \enddot
 */
enum class MASM_EXPORT ChiralStatePreservation {
  //! Don't try to preserve chiral state
  None,
  /*!
   * Use only completely unambiguous zero-effort mappings. Those are the green
   * edges in the graphs below. Note that the ligand gain situation from square
   * planar to square pyramidal is not unique, and therefore not shown as
   * green.
   */
  EffortlessAndUnique,
  /*!
   * Propagates if the best shape mapping is unique, i.e. there are no other
   * mappings with the same quality measures. This enables all green and black
   * edges.
   */
  Unique,
  /*!
   * Chooses randomly from the set of best mappings, permitting chiral state
   * propagation in all cases. So propagating chiral state from square planar
   * to square pyramidal is now possible -- there are two ways of placing the
   * new apical ligand -- but you only get one of them.
   */
  RandomFromMultipleBest
};

/**
 * @brief Influences the choice of shape in substituent additions and
 *   removals that lead to increases or decreases of ligand size
 */
enum class MASM_EXPORT ShapeTransition {
  /*!
   * Try to infer a shape from graph information first. Supplant with best
   * shape for chiral state preservation.
   */
  PrioritizeInferenceFromGraph,
  /*!
   * Always choose the shape so that the maximum amount of chiral
   * information can be preserved (as per ChiralStatePreservation setting).
   */
  MaximizeChiralStatePreservation
};

/**
 * @brief Contains all global settings for the library
 */
struct MASM_EXPORT Options {
  /*!
   * @brief Sets the temperature regime to be used for all Molecules
   *
   * Defaults to high temperature approximation.
   */
  static TemperatureRegime temperatureRegime;
  /*!
   * @brief Sets the manner in which chiral state is preserved for all Molecules
   *
   * Defaults to EffortlessAndUnique.
   */
  static ChiralStatePreservation chiralStatePreservation;

  /**
   * @brief Specifies AtomStereopermutator shape choice behavior on ligand
   *   additions or removals
   *
   * Defaults to MaximizeChiralStatePreservation
   */
  static ShapeTransition shapeTransition;
};

} // namespace molassembler
} // namespace Scine

#endif
