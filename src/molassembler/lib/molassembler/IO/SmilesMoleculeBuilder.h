#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_MOLECULE_BUILDER_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_MOLECULE_BUILDER_H

#include "molassembler/IO/SmilesParseData.h"
#include "molassembler/Types.h"
#include "molassembler/Graph/InnerGraph.h"
#include "boost/variant.hpp"
#include <vector>
#include <stack>

namespace Scine {

namespace Utils {
enum class ElementType : unsigned;
} // namespace Utils


namespace molassembler {

class Molecule;

namespace IO {

struct MoleculeBuilder {
  static bool isValenceFillElement(Utils::ElementType e);

  static unsigned valenceFillElementImplicitHydrogenCount(
    const int valence,
    Utils::ElementType e
  );

  static BondType mutualBondType(
    const boost::optional<BondType>& a,
    const boost::optional<BondType>& b
  );

  // On atom addition
  void addAtom(const AtomData& atom);

  void addRingClosure(const BondData& bond);

  // Trigger on branch open
  inline void branchOpen() {
    vertexStack.push(vertexStack.top());
  }

  inline void branchClose() {
    assert(!vertexStack.empty());
    vertexStack.pop();
  }

  // Trigger on dot bond parse
  inline void setNextAtomUnbonded() {
    lastBondData = SimpleLastBondData::Unbonded;
  }

  // Triggered on non-default bond information after atom addition
  inline void setNextAtomBondInformation(const BondData& bond) {
    lastBondData = bond;
  }

  void setAtomStereo(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<InnerGraph::Vertex>& indexInComponentMap
  );

  void setBondStereo(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<InnerGraph::Vertex>& indexInComponentMap
  );

  // Interpret the graph as possibly distinct molecules
  std::vector<Molecule> interpret();

  enum class SimpleLastBondData {
    Unbonded,
    Unspecified
  };

  //! State for last stored bond data
  boost::variant<SimpleLastBondData, BondData> lastBondData = SimpleLastBondData::Unbonded;

  //! Possibly disconnected tracking graph
  InnerGraph graph;

  //! State to track the vertex a new vertex is bound to
  std::stack<InnerGraph::Vertex> vertexStack;

  //! Storage for bonds marked with stereo indicators ("/" and "\")
  using StereoMarkedBondTuple = std::tuple<InnerGraph::Vertex, InnerGraph::Vertex, BondData::StereoMarker>;
  std::vector<StereoMarkedBondTuple> stereoMarkedBonds;

  //! Storage for ring closure bond indicators
  std::unordered_map<
    unsigned,
    std::pair<InnerGraph::Vertex, boost::optional<BondType>>
  > ringClosures;

  //! AtomData for each created vertex
  std::vector<AtomData> vertexData;
};

} // namespace IO
} // namespace molassembler
} // namespace Scine

#endif
