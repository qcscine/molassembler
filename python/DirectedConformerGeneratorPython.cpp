/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/eigen.h"
#include "pybind11/functional.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/DirectedConformerGenerator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/DistanceGeometry/Error.h"

#include "Utils/Geometry/AtomCollection.h"

using namespace Scine::Molassembler;

using ConformerVariantType = boost::variant<
  Scine::Utils::PositionCollection,
  DgError
>;
extern ConformerVariantType variantCast(outcome::result<Scine::Utils::PositionCollection>);

void init_directed_conformer_generator(pybind11::module& m) {
  pybind11::class_<DirectedConformerGenerator> dirConfGen(
    m,
    "DirectedConformerGenerator",
    R"delim(
      Helper type for directed conformer generation.

      Generates new combinations of BondStereopermutator assignments
      and provides helper functions for the generation of conformers using these
      combinations and the reverse, finding the combinations from conformers.

      It is important that you lower your expectations for the modeling of
      dihedral energy minima, however. Considering that Molassembler neither
      requires you to supply a correct graph, never detects or kekulizes
      aromatic systems nor asks you to supply an overall charge for a molecule,
      it should be understandable that the manner in which Molassembler decides
      where dihedral energy minima are is somewhat underpowered. The manner in
      which shape vertices are aligned in stereopermutation enumeration isn't
      even strictly based on a physical principle. We suggest the following to
      make the most of what the library can do for you:

      - Read the documentation for the various alignments. Consider using not
        just the default
        :class:`~scine_molassembler.BondStereopermutator.Alignment.Staggered`
        alignment, but either
        :class:`~scine_molassembler.BondStereopermutator.Alignemnt.EclipsedAndStaggered`
        or
        :class:`~scine_molassembler.BondStereopermutator.Alignment.BetweenEclipsedAndStaggered`
        to improve your chances of capturing all rotational minima. This will
        likely generate more conformers than strictly required, but should
        capture all minima.
      - Energy minimize all generated conformers with a suitable method and
        then deduplicate.
      - Consider using the
        :class:`~scine_molassembler.DirectedConformerGenerator.Relabeler` to do
        a final deduplication step.

      >>> butane = io.experimental.from_smiles("CCCC")
      >>> generator = DirectedConformerGenerator(butane)
      >>> assert generator.bond_list()
      >>> conformers = []
      >>> while generator.decision_list_set_size() < generator.ideal_ensemble_size():
      ...     conformers.append(
      ...       generator.generate_random_conformation(
      ...         generator.generate_decision_list()
      ...       )
      ...     )
      >>> assert len(conformers) == generator.ideal_ensemble_size()
    )delim"
  );

  pybind11::enum_<DirectedConformerGenerator::IgnoreReason>(
    dirConfGen,
    "IgnoreReason",
    "Reason why a bond is ignored for directed conformer generation"
  ).value("AtomStereopermutatorPreconditionsUnmet", DirectedConformerGenerator::IgnoreReason::AtomStereopermutatorPreconditionsUnmet, "There is not an assigned stereopermutator on both ends of the bond")
    .value("HasAssignedBondStereopermutator", DirectedConformerGenerator::IgnoreReason::HasAssignedBondStereopermutator, "There is already an assigned bond stereopermutator on the bond")
    .value("HasTerminalConstitutingAtom", DirectedConformerGenerator::IgnoreReason::HasTerminalConstitutingAtom, "At least one consituting atom is terminal")
    .value("InCycle", DirectedConformerGenerator::IgnoreReason::InCycle, "The bond is in a cycle (see C++ documentation for details why cycle bonds are excluded)")
    .value("IsEtaBond", DirectedConformerGenerator::IgnoreReason::IsEtaBond, "The bond is an eta bond")
    .value("RotationIsIsotropic", DirectedConformerGenerator::IgnoreReason::RotationIsIsotropic, "Rotation around this bond is isotropic (at least one side's rotating substituents all have the same ranking)");

  dirConfGen.def_static(
    "consider_bond",
    &DirectedConformerGenerator::considerBond,
    pybind11::arg("bond_index"),
    pybind11::arg("molecule"),
    pybind11::arg("alignment") = BondStereopermutator::Alignment::Staggered,
    R"delim(
      Decide whether to consider a bond's dihedral for directed conformer
      generation or not. Returns either an IgnoreReason or an unowned
      stereopermutator instance.

      :param bond_index: Bond index to consider
      :param molecule: The molecule in which bond_index is valid
      :param alignment: Alignment to generate BondStereopermutator instances
        with. Affects stereopermutation counts.
    )delim"
  );

  dirConfGen.def_static(
    "distance",
    &DirectedConformerGenerator::distance,
    pybind11::arg("decision_list_a"),
    pybind11::arg("decision_list_b"),
    pybind11::arg("bounds"),
    R"delim(
      Calculates a distance metric between two decision lists for dihedral
      permutations given bounds on the values at each position of the decision
      lists.

      :param decision_list_a: The first decision list
      :param decision_list_b: The second decision list
      :param bounds: Value bounds on each entry in the decision lists
    )delim"
  );

  dirConfGen.def(
    pybind11::init<Molecule, BondStereopermutator::Alignment, const DirectedConformerGenerator::BondList&>(),
    pybind11::arg("molecule"),
    pybind11::arg("alignment") = BondStereopermutator::Alignment::Staggered,
    pybind11::arg("bonds_to_consider") = std::vector<BondIndex> {},
    R"delim(
      Construct a generator for a particular molecule.

      :param molecule: For which molecule to construct a generator
      :param alignment: Alignment with which to generate BondStereopermutator
        instances on considered bonds
      :param bonds_to_consider: List of bonds that should be considered for
        directed conformer generation. Bonds for which consider_bond yields an
        IgnoreReason will still be ignored.
    )delim"
  );

  dirConfGen.def(
    "generate_decision_list",
    [](DirectedConformerGenerator& cg) {
      return cg.generateNewDecisionList();
    },
    R"delim(
      Generate a new list of discrete dihedral arrangement choices. Guarantees
      that the new list is not yet part of the underlying set. Inserts the
      generated list into the underlying set. Will not generate the same
      decision list twice.

      .. note::
         This function advances ``molassembler``'s global PRNG state.
    )delim"
  );

  dirConfGen.def(
    "insert",
    &DirectedConformerGenerator::insert,
    pybind11::arg("decision_list"),
    R"delim(
      Add a decision list to the underlying set-like data structure.

      :param decision_list: Decision list to insert into the underlying data structure.
    )delim"
  );

  dirConfGen.def(
    "contains",
    &DirectedConformerGenerator::contains,
    pybind11::arg("decision_list"),
    R"delim(
      Checks whether a particular decision list is part of the underlying set

      :param decision_list: Decision list to check for in the underlying data structure
    )delim"
  );

  dirConfGen.def_property_readonly(
    "bond_list",
    &DirectedConformerGenerator::bondList,
    R"delim(
      Get a list of considered bond indices. These are the bonds for which no
      ignore reason was found at construction-time.
    )delim"
  );

  dirConfGen.def(
    "decision_list_set_size",
    &DirectedConformerGenerator::decisionListSetSize,
    R"delim(
      The number of conformer decision lists stored in the underlying set-liked
      data structure
    )delim"
  );

  dirConfGen.def_property_readonly(
    "ideal_ensemble_size",
    &DirectedConformerGenerator::idealEnsembleSize,
    "Returns the number of conformers needed for a full ensemble"
  );


  dirConfGen.def(
    "generate_random_conformation",
    [](
      DirectedConformerGenerator& generator,
      const DirectedConformerGenerator::DecisionList& decisionList,
      const DistanceGeometry::Configuration& configuration
    ) -> ConformerVariantType {
      return variantCast(
        generator.generateRandomConformation(decisionList, configuration)
      );
    },
    pybind11::arg("decision_list"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    R"delim(
      Try to generate a conformer for a particular decision list.

      :param decision_list: Decision list to use in conformer generation
      :param configuration: Distance geometry configurations object. Defaults
        are usually fine.

      .. note::
         This function advances ``molassembler``'s global PRNG state.
    )delim"
  );

  dirConfGen.def(
    "generate_conformation",
    [](
      DirectedConformerGenerator& generator,
      const DirectedConformerGenerator::DecisionList& decisionList,
      const unsigned seed,
      const DistanceGeometry::Configuration& configuration
    ) -> ConformerVariantType {
      return variantCast(
        generator.generateConformation(decisionList, seed, configuration)
      );
    },
    pybind11::arg("decision_list"),
    pybind11::arg("seed"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    R"delim(
      Try to generate a conformer for a particular decision list.

      :param decision_list: Decision list to use in conformer generation
      :param seed: Seed to initialize a PRNG with for use in conformer
        generation.
      :param configuration: Distance geometry configurations object. Defaults
        are usually fine.
    )delim"
  );

  dirConfGen.def(
    "conformation_molecule",
    &DirectedConformerGenerator::conformationMolecule,
    pybind11::arg("decision_list"),
    R"delim(
      Yields a molecule whose bond stereopermutators are set for a particular
      decision list.

      :param decision_list: List of assignments for the considered bonds of the
        generator.
    )delim"
  );

  dirConfGen.def(
    "get_decision_list",
    pybind11::overload_cast<const Scine::Utils::AtomCollection&, BondStereopermutator::FittingMode>(
      &DirectedConformerGenerator::getDecisionList,
      pybind11::const_
    ),
    pybind11::arg("atom_collection"),
    pybind11::arg("fitting_mode") = BondStereopermutator::FittingMode::Nearest,
    R"delim(
      Infer a decision list for the relevant bonds from positions.

      For all bonds considered relevant (i.e. all bonds in bond_list()), fits
      supplied positions to possible stereopermutations and returns the result.
      Entries have a value equal to ``UNKNOWN_DECISION`` if no permutation
      could be recovered. The usual BondStereopermutator fitting tolerances
      apply.

      Assumes several things about your supplied positions:
      - There have only been dihedral changes
      - No atom stereopermutator assignment changes
      - No constitutional rearrangements

      This variant of get_decision_lists checks that the element type sequence
      matches that of the underlying molecule, which holds for conformers
      generated using the underlying molecule.

      :param atom_collection: Positions from which to interpret the decision
        list from.
      :param fitting_mode: Mode altering how decisions are fitted.
    )delim"
  );

  dirConfGen.def(
    "get_decision_list",
    pybind11::overload_cast<const Scine::Utils::PositionCollection&, BondStereopermutator::FittingMode>(
      &DirectedConformerGenerator::getDecisionList,
      pybind11::const_
    ),
    pybind11::arg("positions"),
    pybind11::arg("fitting_mode") = BondStereopermutator::FittingMode::Nearest,
    R"delim(
      Infer a decision list for the relevant bonds from positions.

      For all bonds considered relevant (i.e. all bonds in bond_list()), fits
      supplied positions to possible stereopermutations and returns the result.
      Entries have a value equal to ``UNKNOWN_DECISION`` if no permutation
      could be recovered. The usual BondStereopermutator fitting tolerances
      apply.

      Assumes several things about your supplied positions:
      - There have only been dihedral changes
      - No atom stereopermutator assignment changes
      - No constitutional rearrangements

      :param atom_collection: Positions from which to interpret the decision
        list from.
      :param fitting_mode: Mode altering how decisions are fitted.
    )delim"
  );

  dirConfGen.def_property_readonly_static(
    "UNKNOWN_DECISION",
    [](pybind11::object /* self */) {
      return DirectedConformerGenerator::unknownDecision;
    },
    R"delim(
      Value set in decision list interpretations from positions if no assignment
      could be recovered.
    )delim"
  );

  pybind11::class_<DirectedConformerGenerator::EnumerationSettings> enumerationSettings(
    dirConfGen,
    "EnumerationSettings",
    R"delim(
      Settings for conformer enumeration
    )delim"
  );

  enumerationSettings.def(
    pybind11::init<>(),
    "Default-initialize enumeration settings"
  );

  enumerationSettings.def_readwrite(
    "dihedral_retries",
    &DirectedConformerGenerator::EnumerationSettings::dihedralRetries,
    "Number of attempts to generate the dihedral decision"
  );

  enumerationSettings.def_readwrite(
    "fitting",
    &DirectedConformerGenerator::EnumerationSettings::fitting,
    "Mode for fitting dihedral assignments"
  );

  enumerationSettings.def_readwrite(
    "configuration",
    &DirectedConformerGenerator::EnumerationSettings::configuration,
    "Configuration for conformer generation scheme"
  );

  enumerationSettings.def(
    "__repr__",
    [](pybind11::object settings) -> std::string {
      const std::vector<std::string> members {
        "dihedral_retries",
        "fitting",
        "configuration"
      };

      std::string repr = "(";
      for(const std::string& member : members) {
        const std::string member_repr = pybind11::str(settings.attr(member.c_str()));
        repr += member + "=" + member_repr + ", ";
      }
      repr.pop_back();
      repr.pop_back();
      repr += ")";
      return repr;
    }
  );

  dirConfGen.def(
    "enumerate",
    &DirectedConformerGenerator::enumerate,
    pybind11::arg("callback"),
    pybind11::arg("seed"),
    pybind11::arg("settings") = DirectedConformerGenerator::EnumerationSettings {},
    pybind11::call_guard<pybind11::gil_scoped_release>(),
    R"delim(
      Enumerate all conformers of the captured molecule

      Clears the stored set of decision lists, then enumerates all conformers
      of the molecule in parallel.

      .. note::
         This function is parallelized and will utilize ``OMP_NUM_THREADS``
         threads. Callback invocations are unsequenced but the arguments are
         reproducible.

      :param callback: Function called with decision list and conformer
        positions for each successfully generated pair.
      :param seed: Randomness initiator for decision list and conformer
        generation
      :param settings: Further parameters for enumeration algorithms
    )delim"
  );

  dirConfGen.def(
    "enumerate_random",
    &DirectedConformerGenerator::enumerateRandom,
    pybind11::arg("callback"),
    pybind11::arg("settings") = DirectedConformerGenerator::EnumerationSettings {},
    pybind11::call_guard<pybind11::gil_scoped_release>(),
    R"delim(
      Enumerate all conformers of the captured molecule

      Clears the stored set of decision lists, then enumerates all conformers
      of the molecule in parallel.

      .. note::
         This function is parallelized and will utilize ``OMP_NUM_THREADS``
         threads. Callback invocations are unsequenced but the arguments are
         reproducible given the same global PRNG state.

      .. note::
         This function advances ``molassembler``'s global PRNG state.

      :param callback: Function called with decision list and conformer
        positions for each successfully generated pair.
      :param settings: Further parameters for enumeration algorithms
    )delim"
  );

  dirConfGen.def(
    "relabeler",
    &DirectedConformerGenerator::relabeler,
    "Generate a Relabeler for the underlying molecule and bonds"
  );

  dirConfGen.def(
    "bin_midpoint_integers",
    &DirectedConformerGenerator::binMidpointIntegers,
    "Relabels a decision list into bin midpoint integers"
  );

  pybind11::class_<DirectedConformerGenerator::Relabeler> relabeler(
    dirConfGen,
    "Relabeler",
    R"delim(
      Functionality for relabeling decision lists of minimized structures

      Determines dihedral bins from true dihedral distributions of minimized
      structures and generates bin membership lists for all processed
      structures.
    )delim"
  );

  relabeler.def_static(
    "density_bins",
    &DirectedConformerGenerator::Relabeler::densityBins,
    pybind11::arg("dihedrals"),
    pybind11::arg("delta"),
    pybind11::arg("symmetry_order") = 1,
    R"delim(
      Simplest density-based binning function

      Generates bins for a set of dihedral values by sorting the dihedral
      values and then considering any values within the delta as part of the
      same bin.

      Returns a list of pairs representing bin intervals. It is not guaranteed
      that the start of the interval is smaller than the end of the interval.
      This is because of dihedral angle periodicity. The boundaries of the first
      interval of the bins can have inverted order to indicate wrapping.

      :param dihedrals: List of dihedral values to bin
      :param delta: Maximum dihedral distance between dihedral values to
        include in the same bin in radians

      :raises RuntimeError: If the passed list of dihedrals is empty

      >>> bins = DirectedConformerGenerator.Relabeler.bins
      >>> bins([0.1, 0.2], 0.1)
      [(0.1, 0.2)]
      >>> bins([0.1, 0.2, 0.4], 0.1)
      [(0.1, 0.2), (0.4, 0.4)]
      >>> bins([0.1, 0.2, 3.1, -3.1], 0.1)  # Inverted boundaries with wrap
      [(3.1, -3.1), (0.1, 0.2)]
    )delim"
  );

  relabeler.def(
    "bins",
    &DirectedConformerGenerator::Relabeler::bins,
    pybind11::arg("delta") = M_PI / 6,
    R"delim(
      Generate bins for all observed dihedrals

      :param delta: Maximum dihedral distance between dihedral values to
        include in the same bin in radians
    )delim"
  );

  relabeler.def(
    "add",
    &DirectedConformerGenerator::Relabeler::add,
    pybind11::arg("positions"),
    "Add a particular position to the set to relabel"
  );

  relabeler.def(
    "bin_indices",
    &DirectedConformerGenerator::Relabeler::binIndices,
    pybind11::arg("bins"),
    R"delim(
      Determine relabeling for all added positions

      Returns a list of bin membership indices for each added structure in
      sequence.

      :param bins: Bin intervals for all observed bonds (see bins function)
    )delim"
  );

  relabeler.def(
    "bin_midpoint_integers",
    &DirectedConformerGenerator::Relabeler::binMidpointIntegers,
    pybind11::arg("bin_indices"),
    pybind11::arg("bins"),
    R"delim(
      Relabel bin indices into the rounded dihedral value of their bin midpoint

      :param bin_indices: All structures' bin indices (see bin_indices)
      :param bins: Bin intervals for all observed bonds (see bins function)
    )delim"
  );

  relabeler.def_readonly(
    "sequences",
    &DirectedConformerGenerator::Relabeler::sequences,
    "Dominant index sequences at each considered bond"
  );

  relabeler.def_readonly(
    "dihedrals",
    &DirectedConformerGenerator::Relabeler::observedDihedrals,
    "Observed dihedrals at each bond in added structures"
  );
}
