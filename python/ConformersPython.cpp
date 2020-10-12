/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/eigen.h"

#include "Molassembler/Conformers.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/DistanceGeometry/Error.h"

using namespace Scine;
using namespace Molassembler;

using ConformerVariantType = boost::variant<Utils::PositionCollection, DgError>;

ConformerVariantType variantCast(outcome::result<Utils::PositionCollection> result) {
  if(result) {
    return std::move(result.value());
  }

  if(result.error().category().name() != Detail::DGError_category().name()) {
    throw std::invalid_argument("Error is not of expected category!");
  }

  return DgError(result.error().value());
}

namespace {

void init_partiality(pybind11::module& dg) {
  pybind11::enum_<DistanceGeometry::Partiality>(
    dg,
    "Partiality",
    "Limit triangle inequality bounds smoothing to a subset of all atoms"
  ).value("FourAtom", DistanceGeometry::Partiality::FourAtom, "Resmooth only after each of the first four atom choices")
    .value("TenPercent", DistanceGeometry::Partiality::TenPercent, "Resmooth for the first 10% of all atoms")
    .value("All", DistanceGeometry::Partiality::All, "Resmooth after each distance choice");
}

void init_configuration(pybind11::module& dg) {
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

  configuration.def(
    "__repr__",
    [](pybind11::object settings) -> std::string {
      const std::vector<std::string> members {
        "partiality",
        "refinement_step_limit",
        "refinement_gradient_target",
        "spatial_model_loosening",
        "fixed_positions"
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
}

void init_error(pybind11::module& dg) {
  pybind11::enum_<DgError> error(
    dg,
    "Error",
    "Things that can go wrong in Distance Geometry"
  );

  error.value(
    "ZeroAssignmentStereopermutators",
    DgError::ZeroAssignmentStereopermutators,
    R"delim(
      The molecule you are trying to generate conformers for has
      zero-assignment stereopermutators, meaning that it is not representable
      in three dimensional space.

      AtomStereopermutators remove those stereopermutations from the
      user-accessible set that it deems obviously impossible. This includes
      overlapping haptic ligand binding cones and multidentate ligand bridges
      that are too short to span the angle needed in a stereopermutation. In
      most cases, this will simply eliminate trans-arranged multidentate
      ligands with too short bridges. It is conservative, however, and may not
      be strict enough to eliminate all stereopermutators with bridges that are
      too short to comply with the spatial modeling. In that case, you may get
      GraphImpossible.

      If you get this error, reconsider whether your input graph is reasonable
      and representable in three dimensions. If you believe the atom
      stereopermutator incorrectly has zero assignments, please contact us and
      open an issue.
    )delim"
  );

  error.value(
    "GraphImpossible",
    DgError::GraphImpossible,
    R"delim(
      The molecule you are trying to generate conformers for is either
      incompatible with the applied spatial model or plain not representable in
      three dimensions.

      The applied spatial model is not very smart and mostly applies simple
      geometric considerations. One one hand, it may be that centers whose
      shapes are heavily distorted due to e.g. multiple small cycles are not
      recognized correctly or modeled loosely enough in order for a conformer
      to be possible. On the other hand, it is also possible to create graphs
      that are not representable in three dimensions. In both circumstances,
      you will get this error to indicate that the spatial model cannot deal
      with your input graph.

      If you get this error, reconsider whether your input graph is reasonable
      and representable in three dimensions. If you are sure it is, please
      contact us and open an issue.
    )delim"
  );

  error.value(
    "RefinementException",
    DgError::RefinementException,
    R"delim(
      An exception occurred during refinement

      The form of the potential during Distance Geometry refinement can be
      exceptional due to e.g. divisions by zero. These exceptions can occur
      completely randomly and have no bearing on validity of input.

      If you get this error, generate some more conformers. If all of your inputs
      yield refinement exceptions, there might be a modeling problem, so please
      contact us and open an issue.
    )delim"
  );

  error.value(
    "RefinementMaxIterationsReached",
    DgError::RefinementMaxIterationsReached,
    R"delim(
      Refinement could not find a minimum in your specified maximum
      number of iterations

      Typically, this may mean that the Molecule you are trying to generate
      conformers is either way too big for the number of iterations in the
      potential minimization or that refinement got stuck in some sort of bad
      back and forwards.

      Try adjusting your number of iterations. If the problem persists, there
      may be a problem with the form of the refinement potential, so please
      contact us and open an issue.
    )delim"
  );

  error.value(
    "RefinedStructureInacceptable",
    DgError::RefinedStructureInacceptable,
    R"delim(
      The result of a refinement did not meet criteria for acceptance

      In Distance Geometry, we generate a list of atom-pairwise distance bounds
      that indicate the minimum and maximum distance that atom pair should have
      in a final conformation. Additionally, chiral constraints are generated
      (similar to improper dihedrals) that indicate whether a chiral element is
      arranged correctly. Refinement, which tries to minimize these errors, may
      end up in a local minimum that still violates some of these bounds.

      This is purely a stochastic problem and should not reflect on your inputs.
      If you get this error, generate some more conformers. If all of your
      inputs yield this error, there may be a problem with the refinement
      potential, so please contact us and open an issue.
    )delim"
  );

  error.value(
    "RefinedChiralsWrong",
    DgError::RefinedChiralsWrong,
    R"delim(
      Chiral constraints on the refined structure are still incorrect

      After refinement, chiral constraints are completely wrong. This is a rare
      stochastic error and should not reflect on your inputs.

      Generate more conformers. If you get this error a lot on your inputs,
      there may be a problem with the refinement potential, so please contact us
      and open an issue.
    )delim"
  );

  error.value(
    "DecisionListMismatch",
    DgError::DecisionListMismatch,
    R"delim(
      In directed conformer generation, failed to generate decision list
    )delim"
  );

  error.value(
    "UnknownException",
    DgError::UnknownException,
    "Unknown exception occurred. Please report this as an issue to the developers!"
  );
}

} // namespace

void init_conformers(pybind11::module& m) {
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

  init_partiality(dg);
  init_configuration(dg);
  init_error(dg);

  dg.def(
    "generate_random_ensemble",
    [](
      const Molecule& molecule,
      unsigned numStructures,
      const DistanceGeometry::Configuration& config
    ) -> std::vector<ConformerVariantType> {
      return Temple::map(
        generateRandomEnsemble(molecule, numStructures, config),
        variantCast
      );
    },
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
      >>> sum([1 if isinstance(r, Error) else 0 for r in results])
      0
    )delim"
  );

  dg.def(
    "generate_ensemble",
    [](
      const Molecule& molecule,
      const unsigned numStructures,
      const unsigned seed,
      const DistanceGeometry::Configuration& config
    ) -> std::vector<ConformerVariantType> {
      return Temple::map(
        generateEnsemble(molecule, numStructures, seed, config),
        variantCast
      );
    },
    pybind11::arg("molecule"),
    pybind11::arg("num_structures"),
    pybind11::arg("seed"),
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
      >>> sum([1 if isinstance(r, Error) else 0 for r in results])
      0
    )delim"
  );

  dg.def(
    "generate_random_conformation",
    [](
      const Molecule& molecule,
      const DistanceGeometry::Configuration& config
    ) -> ConformerVariantType {
      return variantCast(generateRandomConformation(molecule, config));
    },
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

      .. note::
         This function advances ``molassembler``'s global PRNG state.

      >>> # Generate a single conformation
      >>> mol = io.experimental.from_smiles("N[C@](Br)(O)F")
      >>> conformation = generate_random_conformation(mol)
      >>> isinstance(conformation, Error) # Did the conformer generation fail?
      False
      >>> type(conformation) # Successful results have matrix type:
      <class 'numpy.ndarray'>
    )delim"
  );

  dg.def(
    "generate_conformation",
    [](
      const Molecule& molecule,
      const unsigned seed,
      const DistanceGeometry::Configuration& config
    ) -> ConformerVariantType {
      return variantCast(
        generateConformation(molecule, seed, config)
      );
    },
    pybind11::arg("molecule"),
    pybind11::arg("seed"),
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
      :param seed: Seed with which to initialize a PRNG with for the conformer
        generation procedure.
      :param configuration: Detailed Distance Geometry settings. Defaults are
        usually fine.
      :rtype: Either a position result or an error string explaining why
        conformer generation failed.

      >>> # Generate a single conformation
      >>> mol = io.experimental.from_smiles("N[C@](Br)(O)F")
      >>> conformation = generate_conformation(mol, 110)
      >>> isinstance(conformation, Error) # Did the conformer generation fail?
      False
      >>> type(conformation) # Successful results have matrix type:
      <class 'numpy.ndarray'>
    )delim"
  );
}
