/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Graph/InnerGraph.h"

#include "temple/Stringify.h"
#include <iostream>

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(splitGraph) {
  InnerGraph methane(4);
  methane.elementType(0) = Utils::ElementType::C;
  methane.elementType(1) = Utils::ElementType::H;
  methane.elementType(2) = Utils::ElementType::H;
  methane.elementType(3) = Utils::ElementType::H;
  auto bridge = methane.addEdge(0, 1, BondType::Single);
  methane.addEdge(0, 2, BondType::Single);
  methane.addEdge(0, 3, BondType::Single);

  auto sides = methane.splitAlongBridge(bridge);

  BOOST_CHECK(
    std::find(
      std::begin(sides.first),
      std::end(sides.first),
      boost::source(bridge, methane.bgl())
    ) != std::end(sides.first)
  );

  BOOST_CHECK(
    std::find(
      std::begin(sides.second),
      std::end(sides.second),
      boost::target(bridge, methane.bgl())
    ) != std::end(sides.second)
  );
}
