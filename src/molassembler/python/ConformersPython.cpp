/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "VariantPython.h"
#include "pybind11/eigen.h"

#include "molassembler/Conformers.h"
#include "molassembler/Molecule.h"

using VariantType = boost::variant<
  Scine::Utils::PositionCollection,
  std::string
>;

std::vector<VariantType> generateEnsemble(
  const Scine::molassembler::Molecule& molecule,
  unsigned numStructures,
  const Scine::molassembler::DistanceGeometry::Configuration& config
) {
  auto ensemble = Scine::molassembler::generateEnsemble(
    molecule,
    numStructures,
    config
  );

  std::vector<VariantType> returnList;
  returnList.reserve(ensemble.size());

  for(auto& positionResult : ensemble) {
    if(positionResult) {
      returnList.emplace_back(
        std::move(positionResult.value())
      );
    } else {
      returnList.emplace_back(
        positionResult.error().message()
      );
    }
  }

  return returnList;
}

VariantType generateConformation(
  const Scine::molassembler::Molecule& molecule,
  const Scine::molassembler::DistanceGeometry::Configuration& config
) {
  auto conformerResult = Scine::molassembler::generateConformation(
    molecule,
    config
  );

  if(conformerResult) {
    return conformerResult.value();
  }

  return conformerResult.error().message();
}

void init_conformers(pybind11::module& m) {
  using namespace Scine::molassembler;

  auto dg = m.def_submodule("dg", "Distance geometry");
  dg.doc() = "Distance geometry submodule";

  pybind11::enum_<DistanceGeometry::Partiality>(
    dg,
    "Partiality",
    "Limit triangle inequality bounds smoothing to a subset of all atoms"
  ).value("FourAtom", DistanceGeometry::Partiality::FourAtom, "Resmooth only after each of the first four atom choices")
    .value("TenPercent", DistanceGeometry::Partiality::TenPercent, "Resmooth for the first 10% of all atoms")
    .value("All", DistanceGeometry::Partiality::All, "Resmooth after each distance choice");

  pybind11::class_<DistanceGeometry::Configuration> configuration(
    dg,
    "Configuration",
    "A configuration object for distance geometry runs with sane defaults"
  );

  configuration.def(
    pybind11::init<>(),
    "Default-initialize a Configuration"
  );

  configuration.def_readwrite(
    "partiality",
    &DistanceGeometry::Configuration::partiality,
    "Choose for how many atoms to re-smooth the distance bounds after a "
    "distance choice. Defaults to four-atom partiality."
  );

  configuration.def_readwrite(
    "refinement_step_limit",
    &DistanceGeometry::Configuration::refinementStepLimit,
    "Sets the maximum number of refinement steps. Defaults to 10'000."
  );

  configuration.def_readwrite(
    "refinement_gradient_target",
    &DistanceGeometry::Configuration::refinementGradientTarget,
    "Sets the gradient at which a refinement is considered complete. Defaults to 1e-5."
  );

  configuration.def_readwrite(
    "spatial_model_loosening",
    &DistanceGeometry::Configuration::spatialModelLoosening,
    "Set loosening factor for spatial model (1.0 is no loosening, 2.0 is strong loosening). Defaults to 1.0."
  );

  configuration.def_readwrite(
    "fixed_positions",
    &DistanceGeometry::Configuration::fixedPositions,
    R"delim(
      Set fixed positions for a subset of atoms in bohr units. Defaults to no fixed positions.

      Any fixed atom must have zero, one or all binding sites fully fixed. No
      individual sites may be partially fixed (i.e. the atoms constituting a
      haptic ligand binding site must be either completely free or completely
      fixed and nothing in between).
    )delim"
  );

  dg.def(
    "generate_ensemble",
    &::generateEnsemble,
    pybind11::arg("molecule"),
    pybind11::arg("num_structures"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    R"delim(
      Generate a set of 3D positions for a molecule.

      In the case of a molecule that does not have unassigned
      stereopermutators, this is akin to generating a conformational ensemble.
      If there are unassigned stereopermutators, these are assigned at random
      (consistent with relative statistical occurrences of stereopermutations)
      for each structure. If, for instance, your molecules contains a single
      unassigned asymmetric tetrahedron atom stereopermutator, the ensemble
      will contain confomers of both assignments, akin to a racemic mixture.

      :param molecule: Molecule to generate positions for. May not contain
        stereopermutators with zero assignments (no feasible stereopermutations).
      :param num_structures: Number of desired structures to generate
      :param configuration: Detailed Distance Geometry settings. Defaults are
        usually fine.
      :rtype: Heterogeneous list of either a position result or an error
        string explaining why conformer generation failed.
    )delim"
  );

  dg.def(
    "generate_conformation",
    &::generateConformation,
    pybind11::arg("molecule"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    R"delim(
      Generate 3D positions for a molecule.

      In the case of a molecule that does not have unassigned
      stereopermutators, this is akin to generating a conformer.
      If there are unassigned stereopermutators, these are assigned at random
      (consistent with relative statistical occurrences of stereopermutations)
      for each structure. If, for instance, your molecules contains a single
      unassigned asymmetric tetrahedron atom stereopermutator, the resulting
      conformation will be one of both assignments.

      :param molecule: Molecule to generate positions for. May not contain
        stereopermutators with zero assignments (no feasible stereopermutations).
      :param configuration: Detailed Distance Geometry settings. Defaults are
        usually fine.
      :rtype: Either a position result or an error string explaining why
        conformer generation failed.
    )delim"
  );
}
