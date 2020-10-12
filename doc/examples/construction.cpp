/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/IO.h"

using namespace molassembler;

/* Default-construct a hydrogen molecule. Molecule enforces the concept that a
 * molecule must consist of at least two atoms, and the molecular graph must
 * be a single connected component (i.e. a Molecule cannot actually represent
 * two disconnected graphs).
 */
Molecule a;

/* Read a single molecule from a file. Currently, only three file formats are
 * supported:
 *
 * - XYZ (element types and positional information),
 * - MOLFile V2000 (element types, bonds and positional information)
 * - MASM (custom serialization format which stores the internal representation of a Molecule instance)
 *
 * This will fail if the library determines that the file actually contains
 * multiple disjoint molecules. The filetype is determined from the file suffix.
 */
Molecule b = IO::read("methane.xyz");

/* Read multiple molecules from a file. This anticipates the possibility that a
 * file may contain multiple disjoint molecules. The filetype is determined
 * from the file suffix.
 */
std::vector<Molecule> c = IO::split("methane_and_ethane.mol");
