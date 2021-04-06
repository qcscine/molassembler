/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Semantic interpreter of the smiles grammar, constructs molecules
 *   from the data accrued by the parser
 */
#ifndef INCLUDE_MOLASSEMBLER_IO_SMILES_MOLECULE_BUILDER_H
#define INCLUDE_MOLASSEMBLER_IO_SMILES_MOLECULE_BUILDER_H

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/IO/SmilesParseData.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "boost/variant.hpp"
#include "boost/bimap.hpp"
#include "boost/graph/subgraph.hpp"
#include <stack>

namespace Scine {

namespace Utils {
enum class ElementType : unsigned;
} // namespace Utils


namespace Molassembler {

class Molecule;

namespace IO {

struct PiSubgraph {
  using BaseGraph = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    boost::no_property,
    boost::property<boost::edge_index_t, int>
  >;
  using Graph = boost::subgraph<BaseGraph>;
  using Vertex = typename BaseGraph::vertex_descriptor;
  using Edge = typename BaseGraph::edge_descriptor;
  using IndexMap = boost::bimap<Vertex, Vertex>;
  using VertexSet = std::unordered_set<Vertex>;

  static bool hasUnpairedElectrons(Vertex i, int charge, const PrivateGraph& g);

  static bool permittedElementType(Utils::ElementType e) {
    switch(Utils::ElementInfo::base(e)) {
      case Utils::ElementType::C:
      case Utils::ElementType::N:
      case Utils::ElementType::O:
      case Utils::ElementType::S:
      case Utils::ElementType::P:
      case Utils::ElementType::As:
      case Utils::ElementType::Sb:
      case Utils::ElementType::Se:
      case Utils::ElementType::Te: {
        return true;
      }
      default:
        return false;
    }
  }

  struct ViableOmissible {
    bool viable;
    bool omissible;
  };

  static ViableOmissible viableOmissible(
    Vertex i,
    const PrivateGraph& component,
    const AtomData& atomData
  );

  PiSubgraph() = default;

  PiSubgraph(
    const std::vector<PrivateGraph::Edge>& edges,
    const std::vector<AtomData>& atoms,
    const PrivateGraph& g
  );

  Vertex findOrAdd(Vertex i) {
    const auto iter = index.left.find(i);
    if(iter == index.left.end()) {
      const Vertex a = boost::add_vertex(graph);
      index.insert(IndexMap::relation(i, a));
      return a;
    }
    return iter->second;
  }

  /* NOTE: non-const because creating subgraphs mutates the parent subgraph
   * instances
   */
  boost::optional<VertexSet> match();

  Graph graph;
  IndexMap index;
  VertexSet omissible;
};

/**
 * @brief Semantic interpreter of the smiles grammar
 *
 * Constructs possibly disconnected graph during parsing, then constructs
 * molecules from accrued data and inference where necessary.
 */
class MoleculeBuilder {
public:
//!@name Parsing triggers
//!@{
  /* @brief Parsing trigger on encountering an atom
   *
   * Adds an atom with the last set bond information.
   */
  void addAtom(AtomData atom);

  //! @brief Parsing trigger on encountering a ring closure marker
  void addRingClosure(const BondData& bond);

  //! @brief Parsing trigger on branch open
  inline void branchOpen() {
    vertexStack.push(vertexStack.top());
  }

  //! @brief Parsing trigger on branch close
  inline void branchClose() {
    assert(!vertexStack.empty());
    vertexStack.pop();
  }

  //! @brief Parsing trigger on finding a dot (molecule separator) in place of a bond
  inline void setNextAtomUnbonded() {
    lastBondData = boost::none;
  }

  //! @brief Parsing trigger on encountering non-default bond information
  inline void setNextAtomBondInformation(const BondData& bond) {
    lastBondData = bond;
  }
//!@}

  //! @brief Interpret the collected graph as (possibly multiple molecules)
  std::vector<Molecule> interpret(const std::string& smiles);

private:
//!@name Static private members
//!@{
  /*! Determines the mutual bond type of two bond type optionals
   *
   * For simplification of ring closure double-checking, where the ring-closure
   * bond type may or may not be specified at either or both markers.
   *
   * @throws std::runtime_error If the bond types are mismatched
   */
  static BondType mutualBondType(
    const boost::optional<BondType>& a,
    const boost::optional<BondType>& b
  );

  /*! Fetches a map to help with the atom chiral markers
   *
   * Maps the order in which substituents were specified in the SMILES onto
   * shape vertices. Note that there is rotational freedom to these.
   */
  static std::vector<Shapes::Vertex> shapeMap(const ChiralData& chiralData);
//!@}

//!@name Private member functions
//!@{
  //! Determine valence-incremented atoms by aromatics in each component
  std::unordered_set<PrivateGraph::Vertex> matchAromatics(
    std::vector<PrivateGraph>& precursors,
    const std::vector<unsigned>& componentMap,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap
  ) const;

  //! Set shapes according to specified charges and stereo markers
  void setShapes(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap
  );

  //! Set atom stereo post-parse and conversion to molecules
  void setAtomStereo(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap,
    const std::string& smiles
  );

  //! Set bond stereo post-parse and conversion to molecules
  void setBondStereo(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap,
    const std::string& smiles
  );

  //! Add bond stereopermutators in aromatic cycles
  void addAromaticBondStereo(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap
  );
//!@}

//!@name Private member data
//!@{
  //! State for last stored bond data
  boost::optional<BondData> lastBondData;

  //! Possibly disconnected tracking graph
  PrivateGraph graph;

  //! State to track the vertex a new vertex is bound to
  std::stack<PrivateGraph::Vertex> vertexStack;

  //! Storage for bonds marked with stereo indicators ("/" and "\")
  using StereoMarkedBondTuple = std::tuple<PrivateGraph::Vertex, PrivateGraph::Vertex, SmilesBondType>;
  std::vector<StereoMarkedBondTuple> stereoMarkedBonds;

  //! Storage for pi-subgraph edges
  std::vector<PrivateGraph::Edge> piSubgraphEdges;

  //! Storage for ring closure bond indicators
  std::unordered_map<
    unsigned,
    std::pair<PrivateGraph::Vertex, boost::optional<BondType>>
  > ringClosures;

  //! AtomData for each created vertex
  std::vector<AtomData> vertexData;
//!@}
};

} // namespace IO
} // namespace Molassembler
} // namespace Scine

#endif
