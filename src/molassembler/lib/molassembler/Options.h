#ifndef INCLUDE_MOLASSEMBLER_OPTIONS_H
#define INCLUDE_MOLASSEMBLER_OPTIONS_H

#include "chemical_symmetries/Symmetries.h"

#include "detail/SharedTypes.h"
#include "detail/AngstromWrapper.h"

namespace molassembler {

/*!
 * Specifies for which temperature regime the Molecule is being modeled.
 * This currently affects only whether nitrogen atoms with a trigonal
 * pyramidal geometry are considered stereocenters. At low temperatures,
 * there is no pyramidal inversion, and these nitrogens are considered
 * valid stereocenters. At high temperatures, only nitrogen geometries
 * that are in particularly strained cycles (3, 4) are considered
 * stereocenters.
 */
enum class TemperatureRegime {Low, High};

/*!
 * Specifies the effects of graph modifications. In case a substituent is
 * added or removed at a stereocenter, an attempt is made to transfer chiral
 * information into the new geometry. How this attempt is made can be altered:
 * - None: Don't try at all (complete loss of information)
 * - Effortless and unique: Use only completely unambiguous zero-effort
 *   mappings. Those are the green edges in the graphs below. Note that the
 *   ligand gain situation from square planar to square pyramidal is not
 *   unique, and therefore not shown as green.
 * - Unique: Propagates if the best symmetry mapping is unique, i.e. there
 *   are no other mappings with the same quality measures. This enables all
 *   green and black edges.
 * - RandomFromMultipleBest: Chooses randomly from the set of best mappings,
 *   permitting chiral state propagation in all cases. So propagating chiral
 *   state from square planar to square pyramidal is now possible -- there
 *   are two ways of placing the new apical ligand -- but you only get one of
 *   them.
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
 *    trigonalpyramidal [label="trigonal pyramidal"];
 *    Tshaped [label="T-shaped"];
 *    tetrahedral [label="tetrahedral"];
 *    squareplanar [label="square planar"];
 *    seesaw [label="seesaw"];
 *    squarepyramidal [label="square pyramidal"];
 *    trigonalbipyramidal [label="trigonal bipyramidal"];
 *    octahedral [label="octahedral"];
 *    pentagonalpyramidal [label="pentagonal pyramidal"];
 *    pentagonalbipyramidal [label="pentagonal bipyramidal"];
 *    trigonalpyramidal -> tetrahedral [color="forestgreen"];
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
 *    trigonalpyramidal [label="trigonal pyramidal"];
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
 *    tetrahedral -> trigonalpyramidal [color="forestgreen", label="any"];
 *    squareplanar -> Tshaped [color="forestgreen", label="any"];
 *    seesaw -> trigonalpyramidal [color="black", label="axial"];
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
  None,
  EffortlessAndUnique,
  Unique,
  RandomFromMultipleBest
};


struct Options {
  /*! Sets the temperature regime to be used for all Molecules
   *
   * Defaults to high temperature approximation.
   */
  static TemperatureRegime temperatureRegime;
  /*! Sets the manner in which chiral state is preserved for all Molecules
   *
   * Defaults to EffortlessAndUnique.
   */
  static ChiralStatePreservation chiralStatePreservation;
};

// Forward-declare Cycles and CNStereocenter
class Cycles;

namespace Stereocenters {
  class CNStereocenter;
} // namespace Stereocenters

/*! Decides whether to keep a stereocenter or not within a temperature regime
 *
 * Criteria applied are:
 * - Minimum of three adjacent indices
 * - If the high-temperature approximation is invoked, trivalent nitrogen
 *   inverts too rapidly to carry stereoinformation (unless part of a cycle
 *   of size 4 or smaller, where strain hinders inversion)
 */
bool disregardStereocenter(
  const Stereocenters::CNStereocenter& stereocenter,
  const Delib::ElementType centralType,
  const Cycles& cycleData,
  const TemperatureRegime temperatureRegimeSetting
);


/*!
 * Fits a stereocenter to a position collection, excluding the seesaw symmetry
 * if a four-coordinate carbon atom is to be fitted to a position collection
 */
void pickyFit(
  Stereocenters::CNStereocenter& stereocenter,
  const GraphType& graph,
  const AngstromWrapper& angstromWrapper,
  const Symmetry::Name expectedSymmetry
);


} // namespace molassembler

#endif
