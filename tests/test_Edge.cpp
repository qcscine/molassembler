#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include "Edge.h"

using namespace MoleculeManip;

BOOST_AUTO_TEST_CASE( edge_all ) {
  unsigned a = 0, b = 2, c = 1, d = 4;
  // instantiate
  Edge instance(std::make_tuple(0, 1, BondType::Single));
  Edge biggerInstance(std::make_tuple(a, b, BondType::Double));
  Edge evenBiggerInstance(std::make_tuple(c, d, BondType::Triple));

  // member lexicographical comparisons
  BOOST_CHECK(instance.smallerThan(a, b));
  BOOST_CHECK(instance.smallerThan(c, d));
  BOOST_CHECK(!instance.greaterThan(a, b));
  BOOST_CHECK(!instance.greaterThan(c, d));

  BOOST_CHECK(biggerInstance.greaterThan(0, 1));
  BOOST_CHECK(biggerInstance.smallerThan(c, d));

  // operators
  BOOST_CHECK(instance < biggerInstance);
  BOOST_CHECK(biggerInstance < evenBiggerInstance);

  BOOST_CHECK(!(
    instance < instance
  ));
  BOOST_CHECK(!(
    biggerInstance < biggerInstance
  ));
  BOOST_CHECK(!(
    evenBiggerInstance < evenBiggerInstance
  ));
}
