/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/eigen.h"

#include "molassembler/Conformers.h"
#include "molassembler/Molecule.h"

using VariantType = boost::variant<
  Scine::Utils::PositionCollection,
  std::string
>;

std::vector<VariantType> generateRandomEnsemble(
  const Scine::molassembler::Molecule& molecule,
  unsigned numStructures,
  const Scine::molassembler::distance_geometry::Configuration& config
) {
  auto ensemble = Scine::molassembler::generateRandomEnsemble(
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

std::vector<VariantType> generateEnsemble(
  const Scine::molassembler::Molecule& molecule,
  const unsigned numStructures,
  const unsigned seed,
  const Scine::molassembler::distance_geometry::Configuration& config
) {
  auto ensemble = Scine::molassembler::generateEnsemble(
    molecule,
    numStructures,
    seed,
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

VariantType generateRandomConformation(
  const Scine::molassembler::Molecule& molecule,
  const Scine::molassembler::distance_geometry::Configuration& config
) {
  auto conformerResult = Scine::molassembler::generateRandomConformation(
    molecule,
    config
  );

  if(conformerResult) {
    return conformerResult.value();
  }

  return conformerResult.error().message();
}

VariantType generateConformation(
  const Scine::molassembler::Molecule& molecule,
  const unsigned seed,
  const Scine::molassembler::distance_geometry::Configuration& config
) {
  auto conformerResult = Scine::molassembler::generateConformation(
    molecule,
    seed,
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
  dg.doc() = R"delim(
    Conformer generation is based on four-spatial dimension Distance Geometry. This
    library's implementation features the following:

    1. A spatial model generates atom-pairwise bounds on their distance in the final
       conformations and four-atom chiral constraints when distance bounds cannot
       encompass chiral elements of complex shapes. For large shapes, chiral
       information is captured by using multiple chiral constraints.
    2. The distance bounds are smoothed to conform to the triangle inequalities.
       After each successive choice of a fixed distance between the bounds, you can
       choose to re-smooth all bounds (full metrization) or stop re-smoothing after
       a fixed number of chosen distances (partial metrization).
    3. The bounds are embedded in four dimensions and refined in three stages,
       permitting the chiral constraints to invert by expanding into four
       dimensions, and then compressing the fourth dimension back out. Lastly,
       dihedral error terms are minimized.
    4. The refinement error function is modified to enable the placement of haptic
       ligand's bonding atoms' average position at shapes' idealized ligand
       sites.
  )delim";

  pybind11::enum_<distance_geometry::Partiality>(
    dg,
    "Partiality",
    "Limit triangle inequality bounds smoothing to a subset of all atoms"
  ).value("FourAtom", distance_geometry::Partiality::FourAtom, "Resmooth only after each of the first four atom choices")
    .value("TenPercent", distance_geometry::Partiality::TenPercent, "Resmooth for the first 10% of all atoms")
    .value("All", distance_geometry::Partiality::All, "Resmooth after each distance choice");

  pybind11::class_<distance_geometry::Configuration> configuration(
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
    &distance_geometry::Configuration::partiality,
    "Choose for how many atoms to re-smooth the distance bounds after a "
    "distance choice. Defaults to four-atom partiality."
  );

  configuration.def_readwrite(
    "refinement_step_limit",
    &distance_geometry::Configuration::refinementStepLimit,
    "Sets the maximum number of refinement steps. Defaults to 10'000."
  );

  configuration.def_readwrite(
    "refinement_gradient_target",
    &distance_geometry::Configuration::refinementGradientTarget,
    "Sets the gradient at which a refinement is considered complete. Defaults to 1e-5."
  );

  configuration.def_readwrite(
    "spatial_model_loosening",
    &distance_geometry::Configuration::spatialModelLoosening,
    "Set loosening factor for spatial model (1.0 is no loosening, 2.0 is strong loosening). Defaults to 1.0."
  );

  configuration.def_readwrite(
    "fixed_positions",
    &distance_geometry::Configuration::fixedPositions,
    R"delim(
      Set fixed positions for a subset of atoms in bohr units. Defaults to no fixed positions.

      Any fixed atom must have zero, one or all binding sites fully fixed. No
      individual sites may be partially fixed (i.e. the atoms constituting a
      haptic ligand binding site must be either completely free or completely
      fixed and nothing in between).
    )delim"
  );

  dg.def(
    "generate_random_ensemble",
    &::generateRandomEnsemble,
    pybind11::arg("molecule"),
    pybind11::arg("num_structures"),
    pybind11::arg("configuration") = distance_geometry::Configuration {},
    R"delim(
      Generate a set of 3D positions for a molecule.

      In the case of a molecule that does not have unassigned
      stereopermutators, this is akin to generating a conformational ensemble.
      If there are unassigned stereopermutators, these are assigned at random
      (consistent with relative statistical occurrences of stereopermutations)
      for each structure. If, for instance, your molecules contains a single
      unassigned asymmetric tetrahedron atom stereopermutator, the ensemble
      will contain confomers of both assignments, akin to a racemic mixture.

      .. note::
         This function is parallelized and will utilize ``OMP_NUM_THREADS``
         threads. The resulting list is sequenced and reproducible given the
         same global PRNG state.

      .. note::
         This function advances ``molassembler``'s global PRNG state.

      :param molecule: Molecule to generate positions for. May not contain
        stereopermutators with zero assignments (no feasible stereopermutations).
      :param num_structures: Number of desired structures to generate
      :param configuration: Detailed Distance Geometry settings. Defaults are
        usually fine.
      :rtype: Heterogeneous list of either a position result or an error
        string explaining why conformer generation failed.

      >>> # Generate a conformational ensemble
      >>> butane = io.experimental.from_smiles("CCCC")
      >>> results = generate_random_ensemble(butane, 10)
      >>> # Each element in the list can be either a string or a positions matrix
      >>> # So let's see how many failed:
      >>> sum([1 if isinstance(r, str) else 0 for r in results])
      0
    )delim"
  );

  dg.def(
    "generate_ensemble",
    &::generateEnsemble,
    pybind11::arg("molecule"),
    pybind11::arg("num_structures"),
    pybind11::arg("seed"),
    pybind11::arg("configuration") = distance_geometry::Configuration {},
    R"delim(
      Generate a set of 3D positions for a molecule.

      In the case of a molecule that does not have unassigned
      stereopermutators, this is akin to generating a conformational ensemble.
      If there are unassigned stereopermutators, these are assigned at random
      (consistent with relative statistical occurrences of stereopermutations)
      for each structure. If, for instance, your molecules contains a single
      unassigned asymmetric tetrahedron atom stereopermutator, the ensemble
      will contain confomers of both assignments, akin to a racemic mixture.

      .. note::
         This function is parallelized and will utilize ``OMP_NUM_THREADS``
         threads. The resulting list is sequenced and reproducible given the
         same seed.

      :param molecule: Molecule to generate positions for. May not contain
        stereopermutators with zero assignments (no feasible stereopermutations).
      :param num_structures: Number of desired structures to generate
      :param configuration: Detailed Distance Geometry settings. Defaults are
        usually fine.
      :rtype: Heterogeneous list of either a position result or an error
        string explaining why conformer generation failed.

      >>> # Generate a conformational ensemble
      >>> butane = io.experimental.from_smiles("CCCC")
      >>> seed = 1010
      >>> results = generate_ensemble(butane, 10, seed)
      >>> # Each element in the list can be either a string or a positions matrix
      >>> # So let's see how many failed:
      >>> sum([1 if isinstance(r, str) else 0 for r in results])
      0
    )delim"
  );

  dg.def(
    "generate_random_conformation",
    &::generateRandomConformation,
    pybind11::arg("molecule"),
    pybind11::arg("configuration") = distance_geometry::Configuration {},
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

      .. note::
         This function advances ``molassembler``'s global PRNG state.

      >>> # Generate a single conformation
      >>> mol = io.experimental.from_smiles("N[C@](Br)(O)F")
      >>> conformation = generate_random_conformation(mol)
      >>> isinstance(conformation, str) # Did the conformer generation fail?
      False
      >>> type(conformation) # Successful results have matrix type:
      <class 'numpy.ndarray'>
    )delim"
  );

  dg.def(
    "generate_conformation",
    &::generateConformation,
    pybind11::arg("molecule"),
    pybind11::arg("seed"),
    pybind11::arg("configuration") = distance_geometry::Configuration {},
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
      :param seed: Seed with which to initialize a PRNG with for the conformer
        generation procedure.
      :param configuration: Detailed Distance Geometry settings. Defaults are
        usually fine.
      :rtype: Either a position result or an error string explaining why
        conformer generation failed.

      >>> # Generate a single conformation
      >>> mol = io.experimental.from_smiles("N[C@](Br)(O)F")
      >>> conformation = generate_conformation(mol, 110)
      >>> isinstance(conformation, str) # Did the conformer generation fail?
      False
      >>> type(conformation) # Successful results have matrix type:
      <class 'numpy.ndarray'>
    )delim"
  );
}
