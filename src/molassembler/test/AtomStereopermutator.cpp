/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "boost/test/unit_test.hpp"

#include "molassembler/Molecule.h"
#include "molassembler/AtomStereopermutator.h"
#include "chemical_symmetries/Symmetries.h"

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(AtomStereopermutatorUpDown) {
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
