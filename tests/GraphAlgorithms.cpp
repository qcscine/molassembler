#define BOOST_TEST_MODULE GraphAlgorithmsTests
#define BOOST_TEST_DYN_LINK

#include "boost/test/unit_test.hpp"
#include "template_magic/Containers.h"

#include "GraphAlgorithms.h"
#include "Molecule.h"
#include "RepeatedElementCollection.h"
#include "StdlibTypeAlgorithms.h"

#include <random>

/* TODO
 */

using namespace MoleculeManip;
using namespace GraphAlgorithms;

using NodeType = BFSVisitors::TreeGenerator::NodeType;

BOOST_AUTO_TEST_CASE( graphAlgorithmTrivialTests ) {

  auto testInstance = Molecule(
    makeRepeatedElementCollection(Delib::ElementType::H, 8),
    Edges({
      {{0, 1}, BondType::Single},
      {{1, 2}, BondType::Single},
      {{1, 4}, BondType::Single},
      {{2, 3}, BondType::Single},
      {{3, 4}, BondType::Single},
      {{4, 5}, BondType::Single},
      {{5, 6}, BondType::Single},
      {{5, 7}, BondType::Single}
    })
  );

  BOOST_CHECK(
    numConnectedComponents(testInstance.getGraph()) == 1
  );

  auto testExpansionCorrectness = [](
    Molecule&& molecule, 
    const std::shared_ptr<NodeType>& comparisonTreePtr
  ) {
    auto treePtr = makeTree(molecule.getGraph());
    bool test = (*treePtr == *comparisonTreePtr);
    BOOST_CHECK(test);
    if(!test) {
      std::cout << "Expected tree: " << std::endl;
      std::cout << comparisonTreePtr << std::endl;
      std::cout << "Got: " << std::endl;
      std::cout << treePtr << std::endl << std::endl;
    }
  };

  // triangle
  testExpansionCorrectness(
    Molecule(
      makeRepeatedElementCollection(Delib::ElementType::H, 3),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 2}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {2ul}),
      Tree::nodePtr(2ul, {1ul})
    })
  );

  // square
  testExpansionCorrectness(
    Molecule(
      makeRepeatedElementCollection(Delib::ElementType::H, 4),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 3}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {3ul}),
      Tree::nodePtr(2ul, {3ul})
    })
  );

  // pentangle
  testExpansionCorrectness(
    Molecule(
      makeRepeatedElementCollection(Delib::ElementType::H, 5),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 4}, BondType::Single},
        {{3, 4}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {
        Tree::nodePtr(3ul, {4ul})
      }),
      Tree::nodePtr(2ul, {
        Tree::nodePtr(4ul, {3ul})
      })
    })
  );

  // spiro for good measure
  testExpansionCorrectness(
    Molecule(
      makeRepeatedElementCollection(Delib::ElementType::H, 7),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 3}, BondType::Single},
        {{3, 4}, BondType::Single},
        {{3, 5}, BondType::Single},
        {{4, 6}, BondType::Single},
        {{5, 6}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {
        Tree::nodePtr(3ul, {
          Tree::nodePtr(4ul, {6ul}),
          Tree::nodePtr(5ul, {6ul})
        })
      }),
      Tree::nodePtr(2ul, {
        Tree::nodePtr(3ul, {
          Tree::nodePtr(4ul, {6ul}),
          Tree::nodePtr(5ul, {6ul})
        })
      })
    })
  );

  // fused triangles
  testExpansionCorrectness(
    Molecule(
      makeRepeatedElementCollection(Delib::ElementType::H, 4),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 3}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {
        Tree::nodePtr(3ul),
        Tree::nodePtr(2ul, {3ul})
      }),
      Tree::nodePtr(2ul, {
        Tree::nodePtr(3ul),
        Tree::nodePtr(1ul, {3ul})
      }),
    })
  );

}

// test order independence of tree generation on sequence in Molecule
BOOST_AUTO_TEST_CASE(treeGenerationOrderIndependence) {
  
  Edges edges({
    {{0, 1}, BondType::Single},
    {{1, 2}, BondType::Single},
    {{1, 5}, BondType::Single},
    {{2, 3}, BondType::Single},
    {{3, 4}, BondType::Single},
    {{4, 5}, BondType::Single},
    {{4, 6}, BondType::Single},
    {{6, 7}, BondType::Single},
    {{6, 8}, BondType::Single}
  });

  auto prototypeTree = makeTree(
    Molecule(
      makeRepeatedElementCollection(Delib::ElementType::H, 9),
      edges
    ).getGraph(), 
    0
  );

  std::vector<unsigned> sequence (edges.size());
  std::iota(
    sequence.begin(),
    sequence.end(),
    0
  );

  std::vector<unsigned> _seeds;
  std::mt19937 _randomEngine;

#ifdef NDEBUG
  std::random_device randomDevice;
  for(unsigned n = 0; n < 5; ++n) _seeds.emplace_back(randomDevice());
#else 
  _seeds.emplace_back(2721813754);
#endif

  std::seed_seq _seedSequence(_seeds.begin(), _seeds.end());
  _randomEngine.seed(_seedSequence);

  unsigned nTests = 100;
  bool pass = true;
  for(unsigned n = 0; n < nTests && pass; n++) {
    std::shuffle(
      sequence.begin(),
      sequence.end(),
      _randomEngine
    );

    // create new Molecule and fill it with edges in the specified manner
    Molecule molecule {
      Delib::ElementType::H,
      Delib::ElementType::H,
      BondType::Single
    };

    for(unsigned i = 2; i < 9; ++i) {
      molecule.addAtom(
        Delib::ElementType::H,
        0,
        BondType::Single
      );
    }

    for(const auto& index : sequence) {
      auto edgesConstIterator = edges.begin();
      std::advance(edgesConstIterator, index);
      molecule.addBond(
        edgesConstIterator -> first.first,
        edgesConstIterator -> first.second,
        BondType::Single
      );
    }

    auto currentTree = makeTree(molecule.getGraph(), 0);

    if(*prototypeTree != *currentTree) {
      pass = false;
      break;
    }
  }

  BOOST_CHECK(pass);
}

using LinkReturnType = std::set<
  std::pair<AtomIndexType, AtomIndexType>
>;

bool testSubstituentConnectivity(
  const Molecule& molecule,
  const AtomIndexType& centralAtom,
  const LinkReturnType& expectedResult
) {
  auto condensePairSet = [](const LinkReturnType& linkSet) -> std::string {
    std::vector<std::string> pairStrings;

    for(const auto& linkedPair : linkSet) {
      pairStrings.emplace_back(
        "{"s + std::to_string(linkedPair.first) + ", "s 
        + std::to_string(linkedPair.second) + "}"s
      );
    }

    return TemplateMagic::condenseIterable(pairStrings);
  };

  auto adjacencies = molecule.getAdjacencies(centralAtom);
  std::set<AtomIndexType> adjacenciesSet {
    adjacencies.begin(),
    adjacencies.end()
  };

  auto result = findSubstituentLinks(
    molecule.getGraph(),
    centralAtom,
    adjacenciesSet
  );

  bool pass = (result == expectedResult);
  
  if(!pass) {
    std::cout << "Expected " << condensePairSet(expectedResult) << ", got: "
      << condensePairSet(result) << std::endl;
  } else {
    std::cout << std::endl;
  }

  return pass;
}



BOOST_AUTO_TEST_CASE(substituentConnectivity) {

  const auto spiroGraph = Molecule(
    makeRepeatedElementCollection(Delib::ElementType::C, 7),
    Edges({
      {{0, 1}, BondType::Single},
      {{0, 2}, BondType::Single},
      {{1, 3}, BondType::Single},
      {{2, 3}, BondType::Single},
      {{3, 4}, BondType::Single},
      {{3, 5}, BondType::Single},
      {{4, 6}, BondType::Single},
      {{5, 6}, BondType::Single}
    })
  );

  testSubstituentConnectivity(
    spiroGraph,
    0,
    LinkReturnType({
      {1, 2}
    })
  );

  testSubstituentConnectivity(
    spiroGraph,
    1,
    LinkReturnType({
      {0, 3}
    })
  );

  testSubstituentConnectivity(
    spiroGraph,
    3,
    LinkReturnType({
      {1, 2},
      {4, 5}
    })
  );

  auto triConnectedLigand  = Molecule(
    makeRepeatedElementCollection(Delib::ElementType::C, 7),
    Edges({
      {{0, 1}, BondType::Single},
      {{0, 2}, BondType::Single},
      {{0, 3}, BondType::Single},
      {{1, 4}, BondType::Single},
      {{1, 5}, BondType::Single},
      {{2, 5}, BondType::Single},
      {{2, 6}, BondType::Single},
      {{3, 4}, BondType::Single},
      {{3, 6}, BondType::Single}
    })
  );

  testSubstituentConnectivity(
    triConnectedLigand,
    0,
    LinkReturnType({
      {1, 2},
      {1, 3},
      {2, 3}
    })
  );
}
