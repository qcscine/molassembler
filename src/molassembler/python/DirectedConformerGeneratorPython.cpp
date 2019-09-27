/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "VariantPython.h"
#include "pybind11/eigen.h"

#include "molassembler/Molecule.h"
#include "molassembler/DirectedConformerGenerator.h"
#include "molassembler/BondStereopermutator.h"

void init_directed_conformer_generator(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<DirectedConformerGenerator> dirConfGen(
    m,
    "DirectedConformerGenerator",
    "Helper type for directed conformer generation.\n\n"
    "Generates guaranteed new combinations of BondStereopermutator assignments "
    "and provides helper functions for the generation of conformers using these"
    "combinations and the reverse, finding the combinations from conformers."
  );

  pybind11::enum_<DirectedConformerGenerator::IgnoreReason>(
    dirConfGen,
    "IgnoreReason",
    "Reason why a bond is ignored for directed conformer generation"
  ).value("AtomStereopermutatorPreconditionsUnmet", DirectedConformerGenerator::IgnoreReason::AtomStereopermutatorPreconditionsUnmet)
    .value("HasAssignedBondStereopermutator", DirectedConformerGenerator::IgnoreReason::HasAssignedBondStereopermutator)
    .value("HasTerminalConstitutingAtom", DirectedConformerGenerator::IgnoreReason::HasTerminalConstitutingAtom)
    .value("InCycle", DirectedConformerGenerator::IgnoreReason::InCycle)
    .value("IsEtaBond", DirectedConformerGenerator::IgnoreReason::IsEtaBond)
    .value("RotationIsIsotropic", DirectedConformerGenerator::IgnoreReason::RotationIsIsotropic);

  dirConfGen.def_static(
    "consider_bond",
    &DirectedConformerGenerator::considerBond,
    pybind11::arg("bond_index"),
    pybind11::arg("molecule"),
    pybind11::arg("smallest_cycle_dict"),
    "Decide whether to consider a bond's dihedral for directed conformer generation or not. Returns either an IgnoreReason or a stereopermutator instance."
  );

  dirConfGen.def_static(
    "distance",
    &DirectedConformerGenerator::distance,
    pybind11::arg("decision_list_a"),
    pybind11::arg("decision_list_b"),
    pybind11::arg("bounds"),
    "Calculates a distance metric between two decision lists for dihedral permutations given bounds on the values at each position of the decision lists."
  );

  dirConfGen.def(
    pybind11::init<Molecule, const DirectedConformerGenerator::BondList&>(),
    pybind11::arg("molecule"),
    pybind11::arg("bonds_to_consider") = std::vector<BondIndex> {},
    "Construct a generator for a particular molecule. You can specify bonds you"
    "want considered for directed conformer generation, though bonds for which"
    "consider_bond yields an IgnoreReason will still be ignored."
  );

  dirConfGen.def(
    "generate_decision_list",
    &DirectedConformerGenerator::generateNewDecisionList,
    "Generate a new list of discrete dihedral arrangement choices. Guarantees"
    "that the new list is not yet part of the underlying set."
  );

  dirConfGen.def(
    "insert",
    &DirectedConformerGenerator::insert,
    "Add a decision list to the underlying set-like data structure"
  );

  dirConfGen.def(
    "contains",
    &DirectedConformerGenerator::contains,
    "Checks whether a particular decision list is part of the underlying set"
  );

  dirConfGen.def(
    "bond_list",
    &DirectedConformerGenerator::bondList,
    "Get a list of considered bond indices"
  );

  dirConfGen.def(
    "decision_list_set_size",
    &DirectedConformerGenerator::decisionListSetSize,
    "The number of conformer decision lists stored in the underlying set-liked data structure"
  );

  dirConfGen.def(
    "ideal_ensemble_size",
    &DirectedConformerGenerator::idealEnsembleSize,
    "The number of conformers needed for a full ensemble"
  );

  using VariantType = boost::variant<
    Scine::Utils::PositionCollection,
    std::string
  >;

  dirConfGen.def(
    "generate_conformation",
    [](
      DirectedConformerGenerator& generator,
      const DirectedConformerGenerator::DecisionList& decisionList,
      const DistanceGeometry::Configuration& configuration
    ) -> VariantType {
      auto result = generator.generateConformation(
        decisionList,
        configuration
      );

      if(result) {
        return result.value();
      }

      return result.error().message();
    },
    pybind11::arg("decision_list"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    "Try to generate a conformer for a particular decision list."
  );

  dirConfGen.def(
    "conformation_molecule",
    &DirectedConformerGenerator::conformationMolecule,
    pybind11::arg("decision_list"),
    "Yields a molecule whose bond stereopermutators are set for a particular decision list"
  );

  dirConfGen.def(
    "get_decision_list",
    &DirectedConformerGenerator::getDecisionList,
    pybind11::arg("positions"),
    "Infer a decision list for relevant bonds from positions"
  );
}
