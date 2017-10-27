#define BOOST_TEST_MODULE StereocenterListTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "StereocenterList.h"
#include "BoundsFromSymmetry.h"

using namespace MoleculeManip;

BOOST_AUTO_TEST_CASE( basicTests ) {
  auto mol = DGDBM::symmetricMolecule(Symmetry::Name::Octahedral);

  // self-consistency
  BOOST_CHECK(mol.getStereocenterList() == mol.getStereocenterList());
}
