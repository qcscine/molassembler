// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_OPTIONS_H
#define INCLUDE_MOLASSEMBLER_OPTIONS_H

#include "chemical_symmetries/Names.h"
#include "molassembler/PRNG.h"

#include "Utils/ElementTypes.h"
#include "molassembler/AngstromWrapper.h"

/*!@file
 *
 * @brief Centralizes the main customization points of the library's behavior.
 */

namespace molassembler {

// Forward-declarations
class OuterGraph;

/*! Randomness source for the entire library
 *
 * This instance supplies the library with randomness. The default seeding
 * behavior is build-type dependent. In debug builds, the seed is fixed. In
 * release builds, the generator is seeded from std::random_device.
 *
 * If you wish to get deterministic behavior in release builds, re-seed the
 * generator.
 *
 * @warning Do not use this instance in any static object's destructor!
 */
random::Engine& randomnessEngine();

/*!
 * Specifies for which temperature regime the Molecule is being modeled.
 * This currently affects only whether nitrogen atoms with a trigonal
 * pyramidal geometry are considered stereopermutators.
 */
enum class TemperatureRegime {
  //! No pyramidal inversion, all nitrogens can be stereopermutators
  Low,
  //! Only nitrogen geometries in particularly strained cycles (3, 4) can be stereopermutators
  High
};

/*!
 * Specifies the effects of graph modifications. In case a substituent is
 * added or removed at a stereopermutator, an attempt is made to transfer chiral
 * information into the new geometry. How this attempt is made can be altered.
 *
 * The following graphs are to illustrate which chiral information transfers
 * are possible under which setting.
 *
 * The first is for the situation of ligand gain. There are very few edges and
 * symmetries in this graph because in nearly all cases, there are multiple
 * best mappings, and propagation in these cases are enabled only in
 * RandomFromMultipleBest, in which case the graph gets quite dense, so it's
 * not shown.
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
 * the group of symmetry ligands from which a ligand is being removed.
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
 *    seesaw -> Tshaped [color="black", label="1 axial"];
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
enum class ChiralStatePreservation {
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
   * Propagates if the best symmetry mapping is unique, i.e. there are no other
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
 * @brief Specifies use of the tau criterion in differentiating between
 *   symmetries of sizes four and five
 *
 * If enabled, enables hard thresholds for tau values to differentiate between
 * some symmetries of sizes four @cite Okuniewski2015 and five @cite Addison1984.
 */
enum class TauCriterion {
  Enable,
  Disable
};

struct Options {
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
   * @brief Specifies whether the tau criterion is used throughout the library
   *   in symmetry fitting.
   *
   * Defaults to Enable
   */
  static TauCriterion tauCriterion;
};

// Forward-declare Cycles and AtomStereopermutator
class Cycles;
class AtomStereopermutator;

/*! Decides whether to keep a stereopermutator or not within a temperature regime
 *
 * Criteria applied are:
 * - Minimum of three adjacent indices
 * - If the high-temperature approximation is invoked, trivalent nitrogen
 *   inverts too rapidly to carry stereoinformation (unless part of a cycle
 *   of size 4 or smaller, where strain hinders inversion)
 */
bool disregardStereopermutator(
  const AtomStereopermutator& stereopermutator,
  Scine::Utils::ElementType centralType,
  const Cycles& cycleData,
  TemperatureRegime temperatureRegimeSetting
);

} // namespace molassembler

#endif
