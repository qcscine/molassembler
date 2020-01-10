/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Graph/InnerGraph.h"

#include "temple/Functional.h"

#include "molassembler/IO/SmilesParser.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(SplitGraph) {
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

BOOST_AUTO_TEST_CASE(BridgeEdges) {
  InnerGraph ethane(8);
  ethane.elementType(0) = Utils::ElementType::C;
  ethane.elementType(1) = Utils::ElementType::C;
  ethane.addEdge(0, 1, BondType::Single);
  for(unsigned i = 2; i < 8; ++i) {
    ethane.elementType(i) = Utils::ElementType::H;
  }
  for(unsigned i = 2; i < 5; ++i) {
    ethane.addEdge(0, i, BondType::Single);
  }
  for(unsigned i = 5; i < 8; ++i) {
    ethane.addEdge(1, i, BondType::Single);
  }

  // All edges must be bridge edges
  BOOST_CHECK(
    temple::all_of(
      boost::make_iterator_range(ethane.edges()),
      [&](const auto& e) -> bool {
        return !ethane.canRemove(e);
      }
    )
  );

  // All carbons must be articulation vertices
  BOOST_CHECK(!ethane.canRemove(0) && !ethane.canRemove(1));

  // And again from SMILES
  auto e = io::experimental::parseSmilesSingleMolecule("CC");
  BOOST_CHECK(e.graph().adjacent(0, 1));
  BOOST_CHECK(!e.graph().canRemove(BondIndex {0, 1}));
}
