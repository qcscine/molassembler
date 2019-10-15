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
    "spatial_model_loosening",
    &DistanceGeometry::Configuration::spatialModelLoosening,
    "Set loosening factor for spatial model (1.0 is no loosening, 2.0 is strong loosening)"
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
