/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "boost/test/unit_test.hpp"

#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Conformers.h"
#include "molassembler/IO.h"
#include "shapes/Data.h"
#include "shapes/PropertyCaching.h"
#include "Utils/Math/QuaternionFit.h"
#include "Utils/Geometry/ElementInfo.h"

#include "Fixtures.h"

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(AtomStereopermutatorUpDown) {
  using namespace shapes;

  // Copy out preservation mode and set for test
  auto prior = Options::chiralStatePreservation;
  Options::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;

  std::vector<
    std::pair<Shape, Shape>
  > upPairs {
    {Shape::SquarePyramid, Shape::Octahedron},
    {Shape::PentagonalPyramid, Shape::PentagonalBipyramid},
    {Shape::Seesaw, Shape::TrigonalBipyramid},
    {Shape::T, Shape::Square},
    {Shape::VacantTetrahedron, Shape::Tetrahedron},
    {Shape::SquarePyramid, Shape::Octahedron}
  };

  for(const auto& pair: upPairs) {
    BOOST_CHECK_MESSAGE(
      AtomStereopermutator::up(pair.first) == pair.second,
      "Expected up(" << name(pair.first) << ") == " << name(pair.second)
        << ", got " << name(AtomStereopermutator::up(pair.first)) << " instead."
    );
  }

  std::vector<
    std::tuple<Shape, unsigned, Shape>
  > downTuples {
    {Shape::SquarePyramid, 4u, Shape::Square},
    {Shape::SquarePyramid, 0u, Shape::Tetrahedron},
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
    AtomStereopermutator::down(Shape::TrigonalPrism, 0) == Shape::Pentagon,
    "Expected down(trig prism, 0) == pentagonal planar, got " << name(AtomStereopermutator::down(Shape::TrigonalPrism, 0)) << " instead."
  );

  // Restore prior state
  Options::chiralStatePreservation = prior;
}

BOOST_FIXTURE_TEST_CASE(LigandAdditionPropagatedStateSuperposable, LowTemperatureFixture) {
  /* For the RSMD argument to work, the mappings must be effortless (i.e. no
   * angular distortion). So we need to make sure the setting is correct.
   */
  BOOST_REQUIRE(Options::chiralStatePreservation == ChiralStatePreservation::EffortlessAndUnique);

  auto trySymmetryPropagation = [](const shapes::Shape source) -> bool {
    assert(shapes::size(source) != shapes::constexpr_properties::maxShapeSize);
    assert(shapes::hasMultipleUnlinkedStereopermutations(source, 0));

    Molecule priorMol {Utils::ElementType::Ru, Utils::ElementType::H, BondType::Single};

    while(priorMol.graph().degree(0) != shapes::size(source)) {
      priorMol.addAtom(
        Utils::ElementInfo::element(priorMol.graph().N()),
        0u,
        BondType::Single
      );
    }

    // Set the source symmetry
    priorMol.setShapeAtAtom(0u, source);

    // Check some preconditions
    auto stereopermutatorOption = priorMol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->getShape() == source);
    BOOST_REQUIRE_MESSAGE(
      stereopermutatorOption->numStereopermutations() > 1
      && stereopermutatorOption->numAssignments() > 1,
      "There are not more than one stereopermutations and assignments for a "
      "shape that was expected to have them per hasMultipleUnlinkedStereopermutations"
    );

    // Assign the stereopermutator at random
    priorMol.assignStereopermutatorRandomly(0u);

    // More checks
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    BOOST_REQUIRE(stereopermutatorOption->getShape() == source);

    auto postMol = priorMol;
    // Transition to a larger shape
    postMol.addAtom(
      Utils::ElementInfo::element(postMol.graph().N()),
      0u,
      BondType::Single
    );

    // Ensure the transition to some larger shape has been made.
    stereopermutatorOption = postMol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(
      shapes::size(stereopermutatorOption->getShape())
      == shapes::size(source) + 1
    );

    /* If the stereopermutator is unassigned, then no state was propagated.
     * This can be algorithmic failure or that there is no suitable mapping
     * for this particular shape and the chosen target shape. We can't
     * really differentiate, so we let those cases go.
     */
    if(!stereopermutatorOption->assigned()) {
      return false;
    }

    /* So now we have some propagated state. Let's generate conformers and
     * quaternion fit it over the first to see whether propagation was right.
     */

    auto priorConformerResult = generateRandomConformation(priorMol);
    BOOST_REQUIRE(priorConformerResult);
    auto postConformerResult = generateRandomConformation(postMol);
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
      "RMSD fit from " << shapes::name(source)
      << " to " << shapes::name(stereopermutatorOption->getShape())
      << " not smaller than 0.5 bohr, is:" << fit.getRMSD()
    );

    if(!passRMSD) {
      io::write(
        "sym-" + std::to_string(shapes::nameIndex(source)) + "-prior.mol",
        priorMol,
        priorConformerResult.value()
      );

      io::write(
        "sym-" + std::to_string(shapes::nameIndex(source)) + "-post.mol",
        postMol,
        postConformerResult.value()
      );
    }

    return true;
  };

  unsigned skippedCount = 0;
  unsigned testedCount = 0;
  for(const shapes::Shape shape : shapes::allShapes) {
    if(shapes::size(shape) > 6) {
      continue;
    }

    if(
      shapes::size(shape) != shapes::constexpr_properties::maxShapeSize
      && shapes::hasMultipleUnlinkedStereopermutations(shape, 0u)
    ) {
      if(trySymmetryPropagation(shape)) {
        ++testedCount;
      } else {
        ++skippedCount;
      }
    }
  }

  BOOST_CHECK_MESSAGE(
    testedCount > 0,
    "Ligand addition propagation test was unable to test any shape "
    "propagations since no propagated shapes were assigned"
  );
}

/* This test needs to happen in the low temperature because otherwise the
 * trigonal bipyramid is thermalized (since none of its ligands are bidentate)
 */
BOOST_FIXTURE_TEST_CASE(AtomStereopermutatorContinuity, LowTemperatureFixture) {
  /* The following shape pairs should be fully reversible starting from any
   * stereopermutation in the first shape, adding a ligand and removing it
   * again.
   *
   * Part of the test is that the target shape is selected from the source
   * shape and the source shape is selected from the source shape.
   */
  std::vector<
    std::pair<shapes::Shape, shapes::Shape>
  > reversibleSymmetryShapes {
    {shapes::Shape::Seesaw, shapes::Shape::TrigonalBipyramid},
    {shapes::Shape::SquarePyramid, shapes::Shape::Octahedron},
    {shapes::Shape::PentagonalPyramid, shapes::Shape::PentagonalBipyramid},
    {shapes::Shape::T, shapes::Shape::Square},
    {shapes::Shape::VacantTetrahedron, shapes::Shape::Tetrahedron},
    {shapes::Shape::TrigonalPyramid, shapes::Shape::TrigonalBipyramid}
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
    const shapes::Shape source,
    const shapes::Shape target
  ) {
    if(shapes::size(source) + 1 != shapes::size(target)) {
      throw std::logic_error("Test is not set up right, symmetries are not increment size apart");
    }

    // Initialize a molecule from a metal and the first substituent
    Molecule mol {Utils::ElementType::Ru, Utils::ElementType::H, BondType::Single};

    // Add different substituents until we reach the source size
    while(mol.graph().degree(0) != shapes::size(source)) {
      mol.addAtom(substituentElements.at(mol.graph().N() - 2), 0u, BondType::Single);
    }

    // Set the starting shape
    mol.setShapeAtAtom(0u, source);

    // Ensure there is a stereopermutator of the correct shape there now
    auto stereopermutatorOption = mol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->getShape() == source);
    BOOST_REQUIRE(stereopermutatorOption->numStereopermutations() > 1);
    BOOST_REQUIRE(stereopermutatorOption->numAssignments() > 1);

    // Assign the stereopermutator.
    mol.assignStereopermutatorRandomly(0u);
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    const unsigned sourceAssignment = stereopermutatorOption->assigned().value();
    BOOST_REQUIRE(stereopermutatorOption->getShape() == source);

    // Transition to the target shape
    const AtomIndex lastAddedIndex = mol.addAtom(
      substituentElements.at(mol.graph().N() - 2),
      0u,
      BondType::Single
    );

    /* Ensure the transition to the target shape has been made and chiral
     * state was preserved in some form
     */
    stereopermutatorOption = mol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    BOOST_REQUIRE(stereopermutatorOption->getShape() == target);

    // Revert the addition
    mol.removeAtom(lastAddedIndex);

    // We expect the same assignment from before
    stereopermutatorOption = mol.stereopermutators().option(0u);
    BOOST_REQUIRE(stereopermutatorOption);
    BOOST_REQUIRE(stereopermutatorOption->getShape() == source);
    BOOST_REQUIRE(stereopermutatorOption->assigned());
    BOOST_CHECK(stereopermutatorOption->assigned().value() == sourceAssignment);
  };

  for(const auto& shapePair : reversibleSymmetryShapes) {
    testSymmetryPair(shapePair.first, shapePair.second);
  }
}
