/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Molecule class interface
 *
 * Contains the Molecule class declaration, which is the central class of the
 * library.
 */

#ifndef INCLUDE_MOLASSEMBLER_MOLECULE_H
#define INCLUDE_MOLASSEMBLER_MOLECULE_H

#include "Molassembler/Graph/GraphInterface.h"
#include "boost/optional.hpp"

#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/Options.h"

namespace Scine {
// External forward declarations
namespace Utils {
using ElementTypeCollection = std::vector<ElementType>;
class AtomCollection;
} // namespace Utils

namespace Molassembler {

// Forward declarations
class Graph;
class StereopermutatorList;
struct RankingInformation;
class AtomStereopermutator;

/*!
 * @brief Models a molecule as a graph (connectivity of atoms) and a list of
 *   stereopermutators.
 *
 * This class models chemical molecules as a combination of a graph
 * and a list of stereopermutators.
 *
 * In the graph, vertices additionally store an element type (and thus
 * represent atoms) and edges store a discretized bond type (and therefore
 * represent bonds).
 *
 * Stereopermutators represent the absolute configuration at an atom or bond in an
 * abstract fashion and do not store coordinate information.
 *
 * @parblock@note
 *   You may be surprised to see that all basic editing of Molecules, even
 *   if it seems to concern the graph only, happens in this class interface. That
 *   is because a graph edit may affect rankings at any stereopermutators in the
 *   molecule due to the algorithm by which substituents are ranked. This means
 *   that for every tiny edit, all stereopermutator substituents are re-ranked and
 *   chiral state, if present, is propagated through a possible ranking change.
 *   For that, the list of stereopermutators is required, which is accessible only
 *   in this class.
 * @endparblock
 *
 * @parblock@note
 *   Some explanation is required to qualify the complexity guarantees. It
 *   is assumed that graphs are sparse, i.e. the average number of substituents
 *   at any atom does not grow with the number of atoms in the graph. Any
 *   function linear in the number of an atom's substituents is therefore of
 *   constant complexity. The variables used in the complexity notations are:
 *   - V = number of atoms/vertices
 *   - E = number of bonds/edges
 *   - A = number of atom stereopermutators
 *   - B = number of bond stereopermutators
 * @note The notations used are the Bachmann-Landau notations:
 *   - @math{O} implies that the function grows asymptotically no faster than
 *   - @math{\Theta} implies that the function grow asymptotically as fast as
 *   - @math{\Omega} implies that the function grows asympotically at least as fast as (Knuth definition)
 * @endparblock
 */
class MASM_EXPORT Molecule final : public GraphInterface {
public:

//!@name Static functions
//!@{
  /**
   * @brief Applies a canonicalization index mapping to an atom collection
   *
   * @param canonicalizationIndexMap Index mapping received from a molecule
   *   canonicalization operation
   * @param atomCollection the atom collection to permute
   *
   * @return A permuted atom collection
   */
  static Utils::AtomCollection applyCanonicalizationMap(
    const std::vector<AtomIndex>& canonicalizationIndexMap,
    const Utils::AtomCollection& atomCollection
  );
//!@}

//!@name Special member functions
//!@{
  /* Rule of five members */
  Molecule(Molecule&& other) noexcept;
  Molecule& operator = (Molecule&& rhs) noexcept;
  Molecule(const Molecule& other);
  Molecule& operator = (const Molecule& rhs);
  ~Molecule() final;

  /*!
   * @brief Default-constructor creates a hydrogen molecule (H2).
   *
   * @complexity{\math{\Theta(1)}}
   */
  Molecule() noexcept;

  /*!
   * @brief Single-element molecule constructor
   *
   * @complexity{@math{\Theta(1)}}
   */
  Molecule(Utils::ElementType element) noexcept;

  /*!
   * @brief Construct a minimal molecule from two element types and a mutual bond type
   *
   * @complexity{@math{\Theta(1)}}
   */
  Molecule(
    Utils::ElementType a,
    Utils::ElementType b,
    BondType bondType = BondType::Single
  ) noexcept;

  /*!
   * @brief Constructs from connectivity alone, inferring the stereopermutators from graph
   *
   * Constructs a molecule from connectivity alone. Local shapes and
   * stereopermutators are inferred from the graph alone.
   *
   * @complexity{\math{O(V + E)} in rankings and stereopermutator
   * instantiations}
   * @throws std::logic_error If the supplied graph has multiple connected
   *   components or there are less than 2 atoms
   */
  explicit Molecule(Graph graph);

  /*! @brief Construct from connectivity and positions
   *
   * Construct an instance from a constituting graph and positional information.
   * Local shapes are deduced from positional information. Stereopermutators
   * are inferred from the graph and assigned using the supplied positional
   * information.
   *
   * @complexity{@math{O(V + B)} where @math{B} is the number of candidate
   * bonds for bond stereopermutators}
   *
   * @param graph The graph from which to construct the molecule
   * @param positions Atom positions to use for the classification of shapes
   *   and interpretation of stereopermutators
   * @param bondStereopermutatorCandidatesOptional If boost::none, all bonds are
   *   candidates for BondStereopermutator. Otherwise, only the specified bonds are
   *   checked for BondStereopermutators.
   *
   * @throws std::logic_error If the supplied graph has multiple connected
   *   components or there are less than 2 atoms
   */
  Molecule(
    Graph graph,
    const AngstromPositions& positions,
    const boost::optional<
      std::vector<BondIndex>
    >& bondStereopermutatorCandidatesOptional = boost::none
  );

  /*! @brief Construct a molecule from underlying data fragments
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @throws std::logic_error If the supplied graph has multiple connected
   *   components or there are less than 2 atoms
   *
   * @warning This function is not intended for library consumers. It is used
   *   internally in implementation details.
   */
  MASM_NO_EXPORT Molecule(
    Graph graph,
    StereopermutatorList stereopermutators,
    boost::optional<AtomEnvironmentComponents> canonicalComponentsOption = boost::none
  );
//!@}

//!@name Modifiers
//!@{
  /*! @brief Adds an atom by attaching it to an existing atom.
   *
   * Adds a new atom, attaching it to an existing atom by a specified bond type.
   *
   * @complexity{Factorial in size of the new shape at @p adjacentTo and
   * linear in re-rankings and propagations of non-terminal positions}
   *
   * @param elementType The element type of the new atom
   * @param adjacentTo The atom to which the new atom is to be attached
   * @param bondType The bond type with which the new atom is to be attached
   *
   * @throws std::out_of_range If adjacentTo is invalid, i.e. >= V()
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  AtomIndex addAtom(
    Utils::ElementType elementType,
    AtomIndex adjacentTo,
    BondType bondType = BondType::Single
  ) final;

  /*!
   * @brief Adds a bond between two existing atoms
   *
   * Adds a bond between two already-existing atoms.
   *
   * @complexity{Factorial in size of the larger new shape at either end of
   * the bond and linear in re-rankings and propagations of non-terminal
   * positions}
   *
   * @param a The first atom index
   * @param b The second atom index
   * @param bondType The bond type with which to connect a and b.
   *
   * @throws std::out_of_range If either atom index is invalid, i.e. >= V()
   * @throws std::logic_error If the atom indices match or the edge already
   *   exists.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  BondIndex addBond(
    AtomIndex a,
    AtomIndex b,
    BondType bondType = BondType::Single
  ) final;

  /*! @brief Add a new BondStereopermutator to the molecule
   *
   * @returns A reference to the added stereopermutator
   * @throws If there is already a bond stereopermutator present
   */
  const BondStereopermutator& addPermutator(
    const BondIndex& bond,
    BondStereopermutator::Alignment alignment = BondStereopermutator::Alignment::Eclipsed
  );

  /** @brief Applies an index permutation to the Molecule state
   *
   * @complexity{@math{\Theta(V + A + B)}}
   *
   * @param permutation A vertex permutation
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*! @brief Sets the stereopermutator assignment at a particular atom
   *
   * This sets the stereopermutator assignment at a specific atom index. For this,
   * a stereopermutator must be instantiated and contained in the StereopermutatorList
   * returned by stereopermutators(). The supplied assignment must be either
   * boost::none or smaller than stereopermutatorPtr->numAssignments().
   *
   * @complexity{@math{O(V)} re-rankings and state propagations}
   *
   * @param a The atom index at which a stereopermutator is to be assigned.
   * @param assignmentOption The new assignment. The special value boost::none
   *   makes the stereopermutator indeterminate. Any indeterminate atom
   *   stereopermutators in a molecule at conformation-generation will be assigned
   *   with a probability according to their sterepermutation's relative
   *   statistical occurence.
   *
   * @throws std::out_of_range If the atom index is invalid, i.e. >= V(),
   *   there is no stereopermutator at this position or the assignment index is
   *   invalid.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators.
   *   Stereopermutators may disappear, change their assignment and number of
   *   assignments, or new stereopermutators can appear as a consequence of the
   *   most minor edit. For procedural safety, consider iterators to
   *   StereopermutatorList members and any stereopermutator state stored
   *   external to a Molecule instance and its members invalidated.
   */
  void assignStereopermutator(
    AtomIndex a,
    const boost::optional<unsigned>& assignmentOption
  );

  /*! @brief Sets the stereopermutator assignment on a bond
   *
   * This sets the stereopermutator assignment at a specific bond index. For
   * this, a stereopermutator must be instantiated and contained in the
   * StereopermutatorList returned by stereopermutators(). The supplied
   * assignment must be either boost::none or smaller than
   * stereopermutatorPtr->numAssignments().
   *
   * @complexity{@math{O(V)} re-rankings and state propagations}
   *
   * @param edge The edge at which a stereopermutator is to be assigned.
   * @param assignmentOption The new assignment. The special value boost::none
   *   makes the stereopermutator indeterminate. Any indeterminate bond
   *   stereopermutators in a molecule at conformation-generation will be
   *   assigned randomly.
   *
   * @throws std::out_of_range If the BondIndex is invalid (i.e. either atom
   *   index >= V()), there is no bond stereopermutator at the supplied edge
   *   or the assignment index is invalid.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators.
   *   Stereopermutators may disappear, change their assignment and number of
   *   assignments, or new stereopermutators can appear as a consequence of the
   *   most minor edit. For procedural safety, consider iterators to
   *   StereopermutatorList members and any stereopermutator state stored
   *   external to a Molecule instance and its members invalidated.
   */
  void assignStereopermutator(
    const BondIndex& edge,
    const boost::optional<unsigned>& assignmentOption
  );

  /*!
   * @brief Assigns a stereopermutator stereopermutation at random
   *
   * This sets the stereopermutator assignment at a specific index, taking relative
   * statistical occurence weights of each stereopermutation into account.
   *
   * @complexity{@math{O(V)} re-rankings and state propagations}
   *
   * @param a The atom index at which the atom stereopermutator is to be assigned
   *   randomly
   * @param engine The PRNG engine to use for random assignment. Defaults to
   *   the library-global random number generator engine
   *
   * @throws std::out_of_range If the atom index is invalid (i.e. is >= V()) or
   *   there is no atom stereopermutator at this bond index.
   *
   * @parblock @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators.
   *   Stereopermutators may disappear, change their assignment and number of
   *   assignments, or new stereopermutators can appear as a consequence of the
   *   most minor edit. For procedural safety, consider iterators to
   *   StereopermutatorList members and any stereopermutator state stored
   *   external to a Molecule instance and its members invalidated.
   * @endparblock
   *
   * @parblock @note This function advances the state of the global PRNG by
   * default if the default argument for @p engine is chosen.
   * @endparblock
   */
  void assignStereopermutatorRandomly(AtomIndex a, Random::Engine& engine = randomnessEngine());

  /*!
   * @brief Assigns a bond stereopermutator to a random assignment
   *
   * @param e The bond index at which to randomly assign a bond stereopermutator
   * @param engine The PRNG engine to use for random assignment. Defaults to
   *   the library-global random number generator engine
   *
   * @complexity{@math{O(V)} re-rankings and state propagations}
   *
   * @throws std::out_of_range If the bond index is invalid (i.e. either atom
   *   index is >= V()) or there is no bond stereopermutator at this bond index.
   *
   * @parblock @note Any molecular edit causes a full re-rank at each
   *   non-terminal atom, and can lead to changes in the list of
   *   stereopermutators. Stereopermutators may disappear, change their
   *   assignment and number of assignments, or new stereopermutators can
   *   appear as a consequence of the most minor edit. For procedural safety,
   *   consider iterators to StereopermutatorList members and any
   *   stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   * @endparblock
   *
   * @parblock @note This function advances the state of the global PRNG by
   * default if the default argument for @p engine is chosen.
   * @endparblock
   */
  void assignStereopermutatorRandomly(const BondIndex& e, Random::Engine& engine = randomnessEngine());

  /** @brief Transform the molecule to a canonical form. Invalidates all atom
   *   and bond indices.
   *
   * @complexity{Theoretically the algorithm falls into the exponential class,
   * but for typical molecules much faster.}
   *
   * @param componentBitmask The components of the molecular graph to include
   *   in the canonicalization procedure.
   *
   * @warning This invalidates all atom indices and bond indices and any
   *   references to constituting members of the molecule.
   *
   * @note Use Molecule::canonicalCompare to compare instances of canonicalized
   *   molecules. If you are using the default value for @p componentBitmask,
   *   Molecule::operator == has a shortcut for fully canonical molecules.
   *
   * @post A call to canonicalComponents() yields @p components supplied here.
   *
   * @return Permutation mapping from old indices to new:
   * @code{.cpp}
   * auto indexMapping = mol.canonicalize();
   * AtomIndex newIndex = indexMapping.at(oldIndex);
   * @endcode
   * You can use this to update invalidated indices.
   */
  std::vector<AtomIndex> canonicalize(
    AtomEnvironmentComponents componentBitmask = AtomEnvironmentComponents::All
  );

  /*! @brief Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * @complexity{@math{O(V + A + B)} stereopermutator updates, re-rankings and
   * propagations}
   *
   * @throws std::out_of_range If the supplied index is invalid, i.e. >= V()
   * @throws std::logic_error If removing the atom disconnects the graph
   *   or removes the final atom.
   *
   * @warning Invalidates **all** atom indices due to renumbering
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void removeAtom(AtomIndex a) final;

  /*! @brief Removes a bond from the graph
   *
   * Removes a bond after checking if removing that bond is safe, i.e. does not
   * disconnect the graph. An example of bonds that can always be removed are
   * ring-closing bonds, since they never disconnect the molecular graph.
   *
   * @complexity{@math{O(V + A + B)} stereopermutator updates, re-rankings and
   * propagations}
   *
   * @throws std::out_of_range If the supplied bond index is invalid, i.e. either
   *   atom index >= V() or the specified bond does not exist.
   * @throws std::logic_error If graph().canRemove() returns false.
   *
   * @note It is not safe to remove a bond just because one of the involved
   *   atoms is terminal, since that atom would then be disconnected from the
   *   rest of the molecule. This function merely removes a bond from the graph.
   *   It is, however, considered safe to remove the terminal vertex, which
   *   involves removing the bond to it.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void removeBond(AtomIndex a, AtomIndex b);

  //!@overload
  void removeBond(const BondIndex& bond) final;

  /*! @brief Removes the BondStereopermutator on a specified edge, if present
   *
   * @returns Whether a permutator was deleted
   */
  bool removePermutator(const BondIndex& bond);

  /*! @brief Changes a bond type. Returns whether the bond already existed
   *
   * Changes the bond type between two atom indices. If the bond does not exist
   * yet, adds the bond.
   *
   * @complexity{@math{\Theta(V)} re-rankings and propagations}
   *
   * @param a The first index of the bond whose type should be changed
   * @param b The second index of the bond whose type should be changed
   * @param bondType The new bond type
   *
   * @return If the bond already existed.
   *
   * @throws out_of_range If a or b are invalid, i.e. >= V()
   * @throws std::logic_error If bondType is specified as BondType::Eta. The
   *   representation of bonding to haptic ligands and its dynamism is handled
   *   internally.
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  bool setBondType(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  ) final;

  /*! @brief Changes an existing atom's element type
   *
   * Changes the element type of an existing atom.
   *
   * @complexity{@math{\Theta(V)} re-rankings and propagations}
   *
   * @param a The atom index of the atom whose element type is to be changed
   * @param elementType The new element type
   *
   * @throws std::out_of_range If a is invalid >= V()
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void setElementType(
    AtomIndex a,
    Utils::ElementType elementType
  ) final;

  /*! @brief Sets the local shape at an atom index
   *
   * This sets the local shape at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a AtomStereopermutator instantiated at this atom index, its
   * underlying shape is altered. If there is no AtomStereopermutator at
   * this index, one is instantiated. In all cases, new or modified
   * stereopermutators are default-assigned if there is only one possible
   * assignment.
   *
   * @complexity{@math{\Theta(V)} re-rankings and propagations}
   *
   * @throws std::out_of_range if the supplied atomic index is invalid
   * @throws std::logic_error if the provided shape is a different size than
   *   that of the existing AtomStereopermutator
   *
   * @note Any molecular edit causes a full re-rank at each non-terminal
   *   atom, and can lead to changes in the list of stereopermutators. Stereopermutators
   *   may disappear, change their assignment and number of assignments, or new
   *   stereopermutators can appear as a consequence of the most minor edit. For
   *   procedural safety, consider iterators to StereopermutatorList members and
   *   any stereopermutator state stored external to a Molecule instance and its
   *   members invalidated.
   */
  void setShapeAtAtom(
    AtomIndex a,
    Shapes::Shape shape
  );

  /*! @brief Change the thermalization of an atom stereopermutator
   *
   * Enables or disables the thermalization of the stereopermutations of an
   * atom stereopermutator. When enabled, a stereopermutator is no longer
   * a stereocenter and has only a single assignment. In conformer generation,
   * a particular stereopermutation is chosen at random.
   *
   * @throws std::out_of_range if the supplied atom index is invalid or
   *   there is no atom stereopermutator at that atom
   */
  void thermalizeStereopermutator(AtomIndex a, bool thermalization = true);
//!@}

//!@name Graph information interface
//!@{
  bool adjacent(AtomIndex a, AtomIndex b) const final;

  boost::optional<BondIndex> bond(AtomIndex a, AtomIndex b) const final;

  BondType bondType(const BondIndex& edge) const final;

  bool canRemove(AtomIndex a) const final;

  bool canRemove(const BondIndex& edge) const final;

  const Cycles& cycles() const final;

  unsigned degree(AtomIndex a) const final;

  Utils::ElementType elementType(AtomIndex a) const final;

  AtomIndex V() const final;

  unsigned E() const final;
//!@}

//!@name Information
//!@{
  /*! @brief Get which components of the graph have been used in canonicalization
   *
   * @complexity{@math{\Theta(1)}}
   */
  boost::optional<AtomEnvironmentComponents> canonicalComponents() const;

  /*! @brief Determines what the local shape at a non-terminal atom ought to
   *   be based on the underlying graph
   *
   * @complexity{@math{\Theta(1)} currently since only basic VSEPR is implemented}
   *
   * Returns the expected shape name at a non-terminal atom by inference
   * from graph information only.
   *
   * Can be extended to apply various levels of theory to determine the
   * shape based on various information contained in the graph, such as
   * ligand field theory.
   *
   * @param index The atom index where the local shape should be determined
   * @param ranking The ranking of all substituents at the supplied index. See
   *   rankPriority()
   *
   * @returns A local shape (of the appropriate size to fit the number of
   *   ligand sites) which may or may not be the shape with the lowest
   *   energy.
   *
   * @throws std::out_of_range If the supplied atomic index is >= V()
   * @throws std::logic_error If the supplied atom is terminal or if the amount
   *   of ligand sites exceeds the size of the largest defined shape.
   *
   * @note Currently applies VSEPR if the element type at the supplied index is
   *   not a transition metal, and returns the first shape of appropriate
   *   size otherwise.
   */
  boost::optional<Shapes::Shape> inferShape(
    AtomIndex index,
    const RankingInformation& ranking
  ) const;

  /*! @brief Returns a graphivz string representation of the molecule
   *
   * @complexity{@math{\Theta(V + E + A + B)}}
   *
   * Creates a graphviz representation of a molecule that can be written into a
   * dotfile and processed with graphviz's `dot` binary to create an image of
   * the molecular graph.
   *
   * Includes tooltip information on stereopermutators when hovering over
   * individual nodes.
   */
  std::string dumpGraphviz() const;

  /*! @brief Provides read-only access to the graph representation
   *
   * @complexity{@math{\Theta(1)}}
   */
  const Graph& graph() const;

  /**
   * @brief Hash function for the molecule
   *
   * Convolutional hash including only the currently canonical components.
   *
   * @parblock
   * @warning Remember that comparing hashes is not a viable alternative for
   * comparing molecules. If molecules have different hashes, they are
   * certainly different, but even if their hashes match, they are not
   * definitively the same due to the reduction of information.
   * @endparblock
   *
   * @parblock
   * @warning Do not compare hashes between molecules with different values of
   *   canonicalComponents().
   * @endparblock
   *
   * @throws std::logic_error If the molecule is acanonical, i.e.
   * canonicalComponents() is None
   *
   * @return A convolutional hash of the currently canonical components.
   */
  std::size_t hash() const;

  /*! @brief Provides read-only access to the list of stereopermutators
   *
   * @complexity{@math{\Theta(1)}}
   */
  const StereopermutatorList& stereopermutators() const;

  /*! @brief Generates stereopermutators from connectivity and positional information
   *
   * Positions are an important source of information for stereopermutators as they
   * will alleviate graph-based shape-determination errors and allow for the
   * determination of stereopermutator assignments through spatial fitting.
   *
   * @complexity{@math{\Theta(S!)} where @math{S} is the largest shape size
   * fitted}
   *
   * @param angstromWrapper Wrapped positions in angstrom length units
   * @param explicitBondStereopermutatorCandidatesOption Permits the specification
   *   of a limited set of bonds on which BondStereopermutator instantiation is
   *   attempted. In Interpret.h, for instance, you can choose not to
   *   instantiate BondStereopermutators below a fractional bond order threshold to
   *   avoid spurious frozen dihedrals. By default, all bonds are candidates.
   *
   * @throws std::out_of_range if a BondIndex in
   *   explicitBondStereopermutatorCandidatesOption does not reference an existing
   *   bond (irrelevant if left default).
   */
  StereopermutatorList inferStereopermutatorsFromPositions(
    const AngstromPositions& angstromWrapper,
    const boost::optional<
      std::vector<BondIndex>
    >& explicitBondStereopermutatorCandidatesOption = boost::none
  ) const;

  //! Returns a command-line interface information string
  std::string str() const;

  /*! @brief Rank substituents of an atom
   *
   * Performs a ranking algorithm that attempts to differentiate branches
   * extending at each substituent atom (haptic ligands are not considered a
   * single unit, rather their component atoms are all individual substituents).
   *
   * Groups substituents into ligands (these are not the same since haptic
   * ligands exist) and ranks those too.
   *
   * @complexity{Theoretically unclear, but for typical cases constant time}.
   *
   * @param a The atom whose substituents are to be ranked
   * @param excludeAdjacent A list of substituent atom indices that should be
   *   excluded from ranking
   * @param positionsOption Positional information can be used to determine
   *   auxiliary stereopermutator assignments that arise in the ranking algorithm
   *   and may have different sub-rankings than at the same position in the
   *   molecule considered on its own. It is preferable to supply this if
   *   positional information is present!
   *
   * @throws std::out_of_range If the supplied atom index is i.e. >= V()
   *
   * @returns a RankingInformation instance that contains all gathered
   * information.
   */
  RankingInformation rankPriority(
    AtomIndex a,
    const std::vector<AtomIndex>& excludeAdjacent = {},
    const boost::optional<AngstromPositions>& positionsOption = boost::none
  ) const;
//!@}

//!@name Comparison
//!@{
  /*! @brief Modular comparison of this Molecule with another, assuming that
   *   both are in a canonical form
   *
   * @complexity{@math{O(V)}}
   *
   * @param other The other canonical molecule to compare against
   * @param componentBitmask The components of an atom's environment to include
   *   in the comparison. You should use the same componentBitmask as when
   *   canonicalizing the molecules you are comparing here. It may be possible
   *   to use a bitmask with fewer components, but certainly not one with more.
   *   May not be None.
   * @returns Whether the molecules are identical (i.e. a special case of
   *   isomorphism in which the vertex mapping is an identity permutation)
   *
   * @throw std::logic_error If @p componentBitmask contains more components
   *   than have been used to canonicalize this molecule instance or @p other.
   */
  bool canonicalCompare(
    const Molecule& other,
    AtomEnvironmentComponents componentBitmask = AtomEnvironmentComponents::All
  ) const;

  /*! @brief Modular comparison of this Molecule with another.
   *
   * This permits detailed specification of which elements of the molecular
   * information you want to use in the comparison.
   *
   * At each atom position, a hash is computed that encompasses all local
   * information that is specified to be used in @p componentsBitmask. This
   * hash is then used during graph isomorphism calculation to avoid finding an
   * isomorphism that does not consider the specified factors.
   *
   * @complexity{@math{O(V_1 \cdot V_2)}}
   *
   * @param other The other molecule to compare against
   * @param componentBitmask Components of an atom's environment to include
   * in isomorphism tests. May not be None.
   *
   * @returns None if the molecules are not isomorphic. Returns an index mapping
   * from @p this to @p other if the molecules are isomorphic.
   *
   * @note The number of stereopermutations that a stereopermutator has is
   * considered part of the shape ComparisonOptions.
   *
   * @note If you choose to discard bond order checking, this merely
   * deactivates bond order hashing and a post-isomorphism-search bond order
   * re-check. Bond order information - if present in the molecule prior to
   * calling this function - is also present in stereopermutator ranking information
   * and hence can influence the number of stereopermutations and the currently
   * set stereopermutation index. This can lead to unexpected but logically
   * consistent comparison behavior.
   *
   * @note This function forwards to canonicalCompare if @p compentBitmask
   * matches the canonical components of both molecules involved.
   */
  boost::optional<std::vector<AtomIndex>> modularIsomorphism(
    const Molecule& other,
    AtomEnvironmentComponents componentBitmask
  ) const;
//!@}

//!@name Operators
//!@{
  /*! @brief Equality operator, performs most strict equality comparison
   *
   * If both molecule instances are fully canonical, calls canonicalCompare().
   * Otherwise calls modularIsomorphism().
   *
   * Implemented as
   * @code{.cpp}
   * if(
   *   canonicalComponents() == AtomEnvironmentComponents::All
   *   && other.canonicalComponents() == AtomEnvironmentComponents::All
   * ) {
   *   return canonicalCompare(other, AtomEnvironmentComponents::All);
   * }
   *
   * return modularIsomorphism(other, AtomEnvironmentComponents::All);
   * @endcode
   */
  bool operator == (const Molecule& other) const;
  //! Inverts Molecule::operator ==
  bool operator != (const Molecule& other) const;
//!@}

private:
  //! Private implementation member
  struct Impl;
  std::unique_ptr<Impl> pImpl_;

  /* Allow access to implementation to editor class that enables more
   * macro-oriented editing as opposed to the low-level editing provided here
   */
  friend struct Editing;
};

} // namespace Molassembler
} // namespace Scine

//! @brief Writes some information about a molecule to a stream
MASM_EXPORT std::ostream& operator << (
  std::ostream& os,
  const Scine::Molassembler::Molecule& molecule
);

#endif
