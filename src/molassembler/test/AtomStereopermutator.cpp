/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "boost/test/unit_test.hpp"

#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Conformers.h"
#include "molassembler/IO.h"
#include "chemical_symmetries/Symmetries.h"
#include "chemical_symmetries/Properties.h"
#include "Utils/QuaternionFit.h"

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(atomStereopermutatorUpDown) {
  using namespace Symmetry;

  // Copy out preservation mode and set for test
  auto prior = Options::chiralStatePreservation;
  Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;

  std::vector<
    std::pair<Name, Name>
  > upPairs {
    {Name::SquarePyramidal, Name::Octahedral},
    {Name::PentagonalPyramidal, Name::PentagonalBiPyramidal},
    {Name::Seesaw, Name::TrigonalBiPyramidal},
    {Name::TShaped, Name::SquarePlanar},
    {Name::CutTetrahedral, Name::Tetrahedral},
    {Name::SquarePyramidal, Name::Octahedral}
  };

  for(const auto& pair: upPairs) {
    BOOST_CHECK_MESSAGE(
      AtomStereopermutator::up(pair.first) == pair.second,
      "Expected up(" << name(pair.first) << ") == " << name(pair.second)
        << ", got " << name(AtomStereopermutator::up(pair.first)) << " instead."
    );
  }

  std::vector<
    std::tuple<Name, unsigned, Name>
  > downTuples {
    {Name::SquarePyramidal, 4u, Name::SquarePlanar},
    {Name::SquarePyramidal, 0u, Name::Tetrahedral},
  };

  for(const auto& tuple: downTuples) {
    BOOST_CHECK_MESSAGE(
      AtomStereopermutator::down(std::get<0>(tuple), std::get<1>(tuple)) == std::get<2>(tuple),
      "Expected down(" << name(std::get<0>(tuple)) << ", " << std::get<1>(tuple) << ") == " << name(std::get<2>(tuple))
        << ", got " << name(AtomStereopermutator::down(std::get<0>(tuple), std::get<1>(tuple))) << " instead."
    );
  }

  Options::chiralStatePreservation = ChiralStatePreservation::Unique;

  BOOST_CHECK_MESSAGE(
    AtomStereopermutator::down(Name::TrigonalPrismatic, 0) == Name::PentagonalPlanar,
    "Expected down(trig prism, 0) == pentagonal planar, got " << name(AtomStereopermutator::down(Name::TrigonalPrismatic, 0)) << " instead."
  );

  // Restore prior state
  Options::chiralStatePreservation = prior;
}

BOOST_AUTO_TEST_CASE(ligandAdditionPropagatedStateSuperposable) {
  /* For the RSMD argument to work, the mappings must be effortless (i.e. no
   * angular distortion). So we need to make sure the setting is correct.
   */
  BOOST_REQUIRE(Options::chiralStatePreservation == ChiralStatePreservation::EffortlessAndUnique);

  auto trySymmetryPropagation = [](const Symmetry::Name source) -> bool {
    assert(Symmetry::size(source) != Symmetry::constexprProperties::maxSymmetrySize);
    assert(Symmetry::hasMultipleUnlinkedAssignments(source, 0));

    Molecule priorMol {Utils::ElementType::Ru, Utils::ElementType::H, BondType::Single};

    while(priorMol.graph().degree(0) != Symmetry::size(source)) {
      priorMol.addAtom(
        static_cast<Utils::ElementType>(priorMol.graph().N()),
        0u,
        BondType::Single
      );
    }

    // Set the source symmetry
    priorMol.setGeometryAtAtom(0u, source);

    // Check some preconditions
    auto stereopermutatorOption = priorMol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->getSymmetry() == source);
    BOOST_REQUIRE_MESSAGE(
      stereopermutatorOption->numStereopermutations() > 1
      && stereopermutatorOption->numAssignments() > 1,
      "There are not more than one stereopermutations and assignments for a "
      "symmetry that was expected to have them per hasMultipleUnlinkedAssignments"
    );

    // Assign the stereopermutator at random
    priorMol.assignStereopermutatorRandomly(0u);

    // More checks
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    BOOST_REQUIRE(stereopermutatorOption->getSymmetry() == source);

    auto postMol = priorMol;
    // Transition to a larger symmetry
    postMol.addAtom(
      static_cast<Utils::ElementType>(postMol.graph().N()),
      0u,
      BondType::Single
    );

    // Ensure the transition to some larger symmetry has been made.
    stereopermutatorOption = postMol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(
      Symmetry::size(stereopermutatorOption->getSymmetry())
      == Symmetry::size(source) + 1
    );

    /* If the stereopermutator is unassigned, then no state was propagated.
     * This can be algorithmic failure or that there is no suitable mapping
     * for this particular symmetry and the chosen target symmetry. We can't
     * really differentiate, so we let those cases go.
     */
    if(!stereopermutatorOption->assigned()) {
      return false;
    }

    /* So now we have some propagated state. Let's generate conformers and
     * quaternion fit it over the first to see whether propagation was right.
     */

    auto priorConformerResult = generateConformation(priorMol);
    BOOST_REQUIRE(priorConformerResult);
    auto postConformerResult = generateConformation(postMol);
    BOOST_REQUIRE(postConformerResult);

    /* Quaternion fit conformations of source and target, matching
     * element types. Check RMSD fit is < 0.5 (bohr). I visually inspected a
     * fit at 0.48 bohr RMSD that was correct.
     */

    Utils::PositionCollection positions = postConformerResult.value();
    // Drop the last row for the new atom so we can fit 1:1
    positions.conservativeResize(positions.rows() - 1, 3);
    Utils::QuaternionFit fit {priorConformerResult.value(), positions};

    bool passRMSD = fit.getRMSD() < 0.5;

    BOOST_CHECK_MESSAGE(
      passRMSD,
      "RMSD fit from " << Symmetry::name(source)
      << " to " << Symmetry::name(stereopermutatorOption->getSymmetry())
      << " not smaller than 0.5 bohr, is:" << fit.getRMSD()
    );

    if(!passRMSD) {
      IO::write(
        "sym-" + std::to_string(Symmetry::nameIndex(source)) + "-prior.mol",
        priorMol,
        priorConformerResult.value()
      );

      IO::write(
        "sym-" + std::to_string(Symmetry::nameIndex(source)) + "-post.mol",
        postMol,
        postConformerResult.value()
      );
    }

    return true;
  };

  unsigned skippedCount = 0;
  unsigned testedCount = 0;
  for(const Symmetry::Name name : Symmetry::allNames) {
    if(
      Symmetry::size(name) != Symmetry::constexprProperties::maxSymmetrySize
      && Symmetry::hasMultipleUnlinkedAssignments(name, 0u)
    ) {
      if(trySymmetryPropagation(name)) {
        ++testedCount;
      } else {
        ++skippedCount;
      }
    }
  }

  BOOST_CHECK_MESSAGE(
    testedCount > 0,
    "Ligand addition propagation test was unable to test any symmetry "
    "propagations since no propagated symmetries were assigned"
  );
}

BOOST_AUTO_TEST_CASE(atomStereopermutatorContinuity) {
  /* The following symmetry pairs should be fully reversible starting from any
   * stereopermutation in the first symmetry, adding a ligand and removing it
   * again.
   *
   * Part of the test is that the target symmetry is selected from the source
   * symmetry and the source symmetry is selected from the source symmetry.
   */
  std::vector<
    std::pair<Symmetry::Name, Symmetry::Name>
  > reversibleSymmetryPairs {
    {Symmetry::Name::Seesaw, Symmetry::Name::TrigonalBiPyramidal},
    {Symmetry::Name::SquarePyramidal, Symmetry::Name::Octahedral},
    {Symmetry::Name::PentagonalPyramidal, Symmetry::Name::PentagonalBiPyramidal},
    {Symmetry::Name::TShaped, Symmetry::Name::SquarePlanar},
    {Symmetry::Name::CutTetrahedral, Symmetry::Name::Tetrahedral},
    {Symmetry::Name::TrigonalPyramidal, Symmetry::Name::TrigonalBiPyramidal}
  };

  const std::array<Utils::ElementType, 7> substituentElements {
    Utils::ElementType::F,
    Utils::ElementType::Cl,
    Utils::ElementType::Br,
    Utils::ElementType::I,
    Utils::ElementType::N,
    Utils::ElementType::C,
    Utils::ElementType::B
  };

  auto testSymmetryPair = [&substituentElements](
    const Symmetry::Name source,
    const Symmetry::Name target
  ) {
    if(Symmetry::size(source) + 1 != Symmetry::size(target)) {
      throw std::logic_error("Test is not set up right, symmetries are not increment size apart");
    }

    // Initialize a molecule from a metal and the first substituent
    Molecule mol {Utils::ElementType::Ru, Utils::ElementType::H, BondType::Single};

    // Add different substituents until we reach the source size
    while(mol.graph().degree(0) != Symmetry::size(source)) {
      mol.addAtom(substituentElements.at(mol.graph().N() - 2), 0u, BondType::Single);
    }

    // Set the starting symmetry
    mol.setGeometryAtAtom(0u, source);

    // Ensure there is a stereopermutator of the correct symmetry there now
    auto stereopermutatorOption = mol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->getSymmetry() == source);
    BOOST_REQUIRE(stereopermutatorOption->numStereopermutations() > 1);
    BOOST_REQUIRE(stereopermutatorOption->numAssignments() > 1);

    // Assign the stereopermutator.
    mol.assignStereopermutatorRandomly(0u);
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    const unsigned sourceAssignment = stereopermutatorOption->assigned().value();
    BOOST_REQUIRE(stereopermutatorOption->getSymmetry() == source);

    // Transition to the target symmetry
    const AtomIndex lastAddedIndex = mol.addAtom(
      substituentElements.at(mol.graph().N() - 2),
      0u,
      BondType::Single
    );

    /* Ensure the transition to the target symmetry has been made and chiral
     * state was preserved in some form
     */
    stereopermutatorOption = mol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    BOOST_REQUIRE(stereopermutatorOption->getSymmetry() == target);

    // Revert the addition
    mol.removeAtom(lastAddedIndex);

    // We expect the same assignment from before
    stereopermutatorOption = mol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->getSymmetry() == source);
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    BOOST_CHECK(stereopermutatorOption->assigned().value() == sourceAssignment);
  };

  for(const auto& symmetryPair : reversibleSymmetryPairs) {
    testSymmetryPair(symmetryPair.first, symmetryPair.second);
  }
}
