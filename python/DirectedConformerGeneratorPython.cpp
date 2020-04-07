/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/eigen.h"

#include "molassembler/Molecule.h"
#include "molassembler/DirectedConformerGenerator.h"
#include "molassembler/BondStereopermutator.h"

#include "Utils/Geometry/AtomCollection.h"

void init_directed_conformer_generator(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<DirectedConformerGenerator> dirConfGen(
    m,
    "DirectedConformerGenerator",
    R"delim(
      Helper type for directed conformer generation.

      Generates guaranteed new combinations of BondStereopermutator assignments
      and provides helper functions for the generation of conformers using these
      combinations and the reverse, finding the combinations from conformers.

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
    pybind11::arg("smallest_cycle_dict"),
    R"delim(
      Decide whether to consider a bond's dihedral for directed conformer
      generation or not. Returns either an IgnoreReason or a stereopermutator
      instance.

      :param bond_index: Bond index to consider
      :param molecule: The molecule in which bond_index is valid
      :param smallest_cycle_dict: Dictionary of smallest cycles in molecule.
        Can be generated by the graph cycles instance.
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
    pybind11::init<Molecule, const DirectedConformerGenerator::BondList&>(),
    pybind11::arg("molecule"),
    pybind11::arg("bonds_to_consider") = std::vector<BondIndex> {},
    R"delim(
      Construct a generator for a particular molecule.

      :param molecule: For which molecule to construct a generator
      :param bonds_to_consider: List of bonds that should be considered for
        directed conformer generation. Bonds for which consider_bond yields an
        IgnoreReason will still be ignored.
    )delim"
  );

  dirConfGen.def(
    "generate_decision_list",
    &DirectedConformerGenerator::generateNewDecisionList,
    R"delim(
      Generate a new list of discrete dihedral arrangement choices. Guarantees
      that the new list is not yet part of the underlying set. Inserts the
      generated list into the underlying set. Will not generate the same
      decision list twice.
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

  dirConfGen.def(
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

  dirConfGen.def(
    "ideal_ensemble_size",
    &DirectedConformerGenerator::idealEnsembleSize,
    "Returns the number of conformers needed for a full ensemble"
  );

  using VariantType = boost::variant<
    Scine::Utils::PositionCollection,
    std::string
  >;

  dirConfGen.def(
    "generate_random_conformation",
    [](
      DirectedConformerGenerator& generator,
      const DirectedConformerGenerator::DecisionList& decisionList,
      const distance_geometry::Configuration& configuration
    ) -> VariantType {
      auto result = generator.generateRandomConformation(
        decisionList,
        configuration
      );

      if(result) {
        return result.value();
      }

      return result.error().message();
    },
    pybind11::arg("decision_list"),
    pybind11::arg("configuration") = distance_geometry::Configuration {},
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
      const distance_geometry::Configuration& configuration
    ) -> VariantType {
      auto result = generator.generateConformation(
        decisionList,
        seed,
        configuration
      );

      if(result) {
        return result.value();
      }

      return result.error().message();
    },
    pybind11::arg("decision_list"),
    pybind11::arg("seed"),
    pybind11::arg("configuration") = distance_geometry::Configuration {},
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
      &DirectedConformerGenerator::getDecisionList
    ),
    pybind11::arg("atom_collection"),
    pybind11::arg("fitting_mode") = BondStereopermutator::FittingMode::Thresholded,
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
      &DirectedConformerGenerator::getDecisionList
    ),
    pybind11::arg("positions"),
    pybind11::arg("fitting_mode") = BondStereopermutator::FittingMode::Thresholded,
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
}
