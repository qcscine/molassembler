/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Editing of molecules past low-level member functions
 */
#ifndef INCLUDE_MOLASSEMBLER_EDITING_H
#define INCLUDE_MOLASSEMBLER_EDITING_H

#include "Molassembler/RankingInformation.h"
#include <tuple>

namespace Scine {
namespace Molassembler {

class Molecule;

namespace Editing {

//! Descriptor for an entire haptic site connected to an atom
using AtomSitePair = std::pair<AtomIndex, SiteIndex>;

/*! @brief Splits a molecule along a bridge edge
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @throws std::logic_error If @p bridge is not a bridge edge, i.e. the graph
 *   is not disconnected by cleaving the indicated bonds. This can be tested
 *   for using @code{.cpp} a.graph().canRemove(bridge) == false) @endcode.
 *
 * @returns A pair of molecules, ordered corresponding to the bond's atoms'
 *   position in the supplied bond index (i.e. .first contains the atom at
 *   bridge.first)
 */
std::pair<Molecule, Molecule> cleave(const Molecule& a, BondIndex bridge);

/*! @brief Splits a molecule at a haptic site
 *
 * @param hapticSite The haptic site to cleave the molecule at. Does not check
 *   if the site is conceptually a bridge of the graph. If it is not, splitting
 *   is not undefined.
 *
 * @note The first molecule always contains the atom to which the haptic site
 *   belonged.
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @throws std::invalid_argument If there is no stereopermutator on the
 *   specified atom.
 *
 */
std::pair<Molecule, Molecule> cleave(const Molecule& a, AtomSitePair hapticSite);

/** @brief Inserts a molecule into a bond of another molecule
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param log The molecule being inserted into
 * @param wedge The molecule being inserted
 * @param logBond The bond in @p log into which @p wedge is being inserted into
 * @param firstWedgeAtom The atom of wedge to bond to the first atom in @p logBond
 * @param secondWedgeAtom The atom of wedge to bond to the second atom in @p logBond
 *
 * @note The bond type of the bond inserted into is reused in the new bonds
 *   to @p firstWedgeAtom and @p secondWedgeAtom.
 *
 * @return A molecule inserted into the bond of another molecule
 */
Molecule insert(
  Molecule log,
  const Molecule& wedge,
  BondIndex logBond,
  AtomIndex firstWedgeAtom,
  AtomIndex secondWedgeAtom
);

/** @brief Fuses two molecules, adding all adjacencies and continuations of one
 * Molecule's atoms to another.
 *
 * Adds all adjacent atoms and continuations of @p bottomAtom in @p bottom to
 * @p topAtom in @p top. @p topAtom's element type is unchanged as it
 * is the 'top' of the superimposition / overlay.
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param top The molecule that is at the top of the superposition.
 * @param bottom The molecule that is set at the bottom of the superposition.
 * @param topAtom An atom index of @p top that is placed 'atop' @p bottom's
 *   @p bottomAtom
 * @param bottomAtom An atom index of @p bottom that is placed 'beneath'
 *   @p top's @p topAtom
 *
 * @return A single molecule fused at the specified atoms
 */
Molecule superpose(
  Molecule top,
  const Molecule& bottom,
  AtomIndex topAtom,
  AtomIndex bottomAtom
);

/**
 * @brief Connect two molecules by substituting away the smaller side of a
 *   pair of bonds of separate molecules
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param left The first molecule
 * @param right The second molecule
 * @param leftBond @p left's bond from which to substitute the lighter part away
 * @param rightBond @p right's bond from which to substitute the lighter part away
 *
 * @note The smaller side is chosen by number of atoms first, then molecular
 *   weight if the number of atoms is equal. Should both sides be equal in
 *   both, which side is picked is undefined.
 */
Molecule substitute(
  const Molecule& left,
  const Molecule& right,
  BondIndex leftBond,
  BondIndex rightBond
);

/**
 * @brief Connect molecules by creating a new bond between two atoms from
 *   separate molecules
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param a The first molecule to connect
 * @param b The second molecule to connect
 * @param aConnectAtom The atom index to connect to in @p a
 * @param bConnectAtom The atom index to connect to in @p b
 * @param bondType The type of the new bond
 *
 * @note This function is symmetric on exchanging the tuples (@p a, @p
 * aConnectAtom), (@p b, @p bConnectAtom), i.e. the order of the molecule and
 * atom index arguments does not matter as long as the atom indices are
 * properly associated to their respective molecule.
 *
 * @return A single molecule connected by the specified bond type at the
 * indicated atoms
 */
Molecule connect(
  Molecule a,
  const Molecule& b,
  AtomIndex aConnectAtom,
  AtomIndex bConnectAtom,
  BondType bondType
);

/**
 * @brief Connects two molecules by connecting multiple atoms from one to a
 *   single atom of the other via single bonds
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param a The molecule the ligand is being connected to
 * @param ligand The ligand molecule being bound to @p complexatingAtom
 * @param complexatingAtom The atom in @p a that @p ligand is being bound to
 * @param ligandBindingAtoms The atoms of @p ligand that should be connected
 *   to @p a
 *
 * @return A joined molecule
 */
Molecule addLigand(
  Molecule a,
  const Molecule& ligand,
  AtomIndex complexatingAtom,
  const std::vector<AtomIndex>& ligandBindingAtoms
);

} // namespace Editing
} // namespace Molassembler
} // namespace Scine

#endif
