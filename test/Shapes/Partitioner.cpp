/* @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details. for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Shapes/Partitioner.h"

using namespace Scine;
using namespace shapes;

BOOST_AUTO_TEST_CASE(Partitions) {
  for(unsigned i = 1; i < 4; ++i) {
    for(unsigned j = 1; j < 4; ++j) {
      Partitioner partitioner {i, j};
      do {
        BOOST_CHECK(Partitioner::isOrderedMapping(partitioner.map()));
      } while(partitioner.next_partition());
    }
  }

  // Test a few specific counts
  auto countPartitions = [](const unsigned S, const unsigned E) -> unsigned {
    Partitioner partitioner {S, E};
    unsigned partitions = 0;
    do {
      ++partitions;
    } while(partitioner.next_partition());
    return partitions;
  };

  BOOST_CHECK_EQUAL(countPartitions(1, 1), 1);
  BOOST_CHECK_EQUAL(countPartitions(1, 2), 1);
  BOOST_CHECK_EQUAL(countPartitions(2, 1), 1);
  BOOST_CHECK_EQUAL(countPartitions(2, 2), 3);
  BOOST_CHECK_EQUAL(countPartitions(2, 3), 10);
}
