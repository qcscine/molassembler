#define BOOST_TEST_MODULE OrderDiscoveryHelperTestModule

#include "boost/test/unit_test.hpp"

#include "molassembler/OrderDiscoveryHelper.h"

using namespace molassembler;
using namespace std::string_literals;

template<typename T>
std::string showDiscoveryState(
  const OrderDiscoveryHelper<T>& orderingHelper
) {
  return (
    "Sets: "s
    + temple::condenseIterable(
      temple::map(
        orderingHelper.getSets(),
        [](const auto& set) -> std::string {
          return "{"s + temple::condenseIterable(set) + "}"s;
        }
      )
    )
    + ", Undecided"s
    + temple::condenseIterable(
      temple::map(
        orderingHelper.getUndecidedSets(),
        [](const auto& set) -> std::string {
          return "{"s + temple::condenseIterable(set) + "}"s;
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

BOOST_AUTO_TEST_CASE(stateTransferTests) {
  OrderDiscoveryHelper<unsigned> knowledge {
    {4u, 9u, 13u, 20u}
  };

  knowledge.addLessThanRelationship(4, 9);
  knowledge.addLessThanRelationship(4, 13);
  knowledge.addLessThanRelationship(4, 20);
  knowledge.addLessThanRelationship(9, 13);
  knowledge.addLessThanRelationship(9, 20);
  knowledge.addLessThanRelationship(13, 20);

  OrderDiscoveryHelper<unsigned> partialMatch {
    {4u, 13u, 20u}
  };

  partialMatch.addRelationshipsFromOther(knowledge);

  BOOST_CHECK_MESSAGE(
    partialMatch.isSmaller(4, 13)
    && partialMatch.isSmaller(4, 20)
    && partialMatch.isSmaller(13, 20),
    "Information not transferred by addRelationshipsFromOther"
  );

  OrderDiscoveryHelper<unsigned> toMerge {
    {10u, 13u}
  };


  toMerge.addLessThanRelationship(10, 13);
  toMerge.addAllFromOther(knowledge);

  BOOST_CHECK_MESSAGE(
    toMerge.isSmaller(13, 20),
    "Information not transferred by addAllFromOther"
  );

  BOOST_CHECK_MESSAGE(
    toMerge.isSmaller(10, 20),
    "Transferability edge from 10 -> 20 not found!"
  );
}
