/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"

#include "molassembler/Conformers.h"
#include "molassembler/Molecule.h"

std::vector<Scine::Utils::PositionCollection> generateEnsemble(
  const Scine::molassembler::Molecule& molecule,
  unsigned numStructures,
  const Scine::molassembler::DistanceGeometry::Configuration& config
) {
  auto ensembleResult = Scine::molassembler::generateEnsemble(
    molecule,
    numStructures,
    config
  );

  if(ensembleResult) {
    return ensembleResult.value();
  } else {
    std::cout << ensembleResult.error();
    throw std::exception();
  }
}

Scine::Utils::PositionCollection generateConformation(
  const Scine::molassembler::Molecule& molecule,
  const Scine::molassembler::DistanceGeometry::Configuration& config
) {
  auto conformerResult = Scine::molassembler::generateConformation(
    molecule,
    config
  );

  if(conformerResult) {
    return conformerResult.value();
  } else {
    std::cout << conformerResult.error();
    throw std::exception();
  }
}

void init_conformers(pybind11::module& m) {
  using namespace Scine::molassembler;

  auto dg = m.def_submodule("dg", "Distance geometry");

  pybind11::enum_<DistanceGeometry::Partiality>(
    dg,
    "Partiality",
    "Limit triangle inequality bounds smoothing to a subset of all atoms"
  ).value("FourAtom", DistanceGeometry::Partiality::FourAtom)
    .value("TenPercent", DistanceGeometry::Partiality::TenPercent)
    .value("All", DistanceGeometry::Partiality::All);

  pybind11::class_<DistanceGeometry::Configuration> configuration(
    dg,
    "Configuration",
    "A configuration object for distance geometry runs with sane defaults"
  );

  configuration.def_readwrite(
    "partiality",
    &DistanceGeometry::Configuration::partiality,
    "Choose for how man atoms to re-smooth the distance bounds after a "
    "distance choice"
  );

  configuration.def_readwrite(
    "refinement_step_limit",
    &DistanceGeometry::Configuration::refinementStepLimit,
    "Sets the maximum number of refinement steps"
  );

  configuration.def_readwrite(
    "refinement_gradient_target",
    &DistanceGeometry::Configuration::refinementGradientTarget,
    "Sets the gradient at which a refinement is considered complete"
  );

  configuration.def_readwrite(
    "failure_ratio",
    &DistanceGeometry::Configuration::failureRatio,
    "Set maximum allowed ratio of failures / (# desired conformers)"
  );

  configuration.def_readwrite(
    "fixed_positions",
    &DistanceGeometry::Configuration::fixedPositions,
    "Set fixed positions for a subset of atoms"
  );

  dg.def(
    "generate_ensemble",
    &::generateEnsemble,
    pybind11::arg("molecule"),
    pybind11::arg("num_structures"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    "Generate an ensemble of conformations of a molecule of a particular size."
    "The molecule you use to generate conformations may not contain "
    "stereopermutators with zero assignments. Ensemble generation may fail, "
    "which will raise an exception."
  );

  dg.def(
    "generate_conformation",
    &::generateConformation,
    pybind11::arg("molecule"),
    pybind11::arg("configuration") = DistanceGeometry::Configuration {},
    "The molecule you use to generate a conformation may not contain "
    "stereopermutators with zero assignments. Ensemble generation may fail, "
    "which will raise an exception."
  );
}
