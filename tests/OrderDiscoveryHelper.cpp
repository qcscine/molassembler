#define BOOST_TEST_MODULE OrderDiscoveryHelperTestModule
#define BOOST_TEST_DYN_LINK

#include "boost/test/unit_test.hpp"

#include "OrderDiscoveryHelper.h"

using namespace MoleculeManip;
using namespace std::string_literals;

template<typename T>
std::string showDiscoveryState(
  const OrderDiscoveryHelper<T>& orderingHelper
) {
  return (
    "Sets: "s 
    + TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& set) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(set) + "}"s;
        }
      )
    ) 
    + ", Undecided"s
    + TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getUndecidedSets(),
        [](const auto& set) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(set) + "}"s;
        }
      )
    ) 
  );
}

BOOST_AUTO_TEST_CASE(sampleTest) {
  auto helper = OrderDiscoveryHelper<unsigned>({1, 2, 3, 6});

  BOOST_CHECK(helper.getSets().size() == 1 && helper.getUndecidedSets().size() == 1);

  helper.addLessThanRelationship(6, 1);
  helper.addLessThanRelationship(6, 2);
  helper.addLessThanRelationship(6, 3);

  BOOST_CHECK_MESSAGE(
    helper.getSets().size() == 2 
    && helper.getUndecidedSets().size() == 1,
    "State: " << showDiscoveryState(helper)
  );

  BOOST_CHECK_MESSAGE(
    !helper.isTotallyOrdered(),
    "Helper believes total order achieved prior to full discovery!"
  );

  helper.addLessThanRelationship(2, 1);
  helper.addLessThanRelationship(3, 1);

  BOOST_CHECK_MESSAGE(
    helper.getSets().size() == 3 
    && helper.getUndecidedSets().size() == 1,
    "State: " << showDiscoveryState(helper)
  );

  BOOST_CHECK_MESSAGE(
    !helper.isTotallyOrdered(),
    "Helper believes total order achieved prior to full discovery!"
  );

  helper.addLessThanRelationship(3, 2);

  BOOST_CHECK_MESSAGE(
    helper.getSets().size() == 4 
    && helper.getUndecidedSets().size() == 0,
    "State: " << showDiscoveryState(helper)
  );

  BOOST_CHECK_MESSAGE(
    helper.isTotallyOrdered(),
    "Does not recognize total order condition after full ordering discovered"
  );

  int x = 4, y = 9, z = 1;

  std::set<const int*> pointerSet {&x, &y, &z};

  auto pointerHelper = OrderDiscoveryHelper<const int*>(pointerSet);

  pointerHelper.addLessThanRelationship(&x, &z);
  pointerHelper.addLessThanRelationship(&x, &y);
  pointerHelper.addLessThanRelationship(&z, &y);

  BOOST_CHECK_MESSAGE(
    pointerHelper.getSets().size() == 3 
    && pointerHelper.getUndecidedSets().size() == 0,
    "State: " << showDiscoveryState(helper)
  );
}
