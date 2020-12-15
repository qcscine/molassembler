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
#include <stack>

namespace Scine {

namespace Utils {
enum class ElementType : unsigned;
} // namespace Utils


namespace Molassembler {

class Molecule;

namespace IO {

/**
 * @brief Semantic interpreter of the smiles grammar, constructs molecules
 *   from the data accrued by the parser
 */
class MoleculeBuilder {
public:
//!@name Parsing triggers
//!@{
  /* @brief Parsing trigger on encountering an atom
   *
   * Adds an atom with the last set bond information.
   */
  void addAtom(const AtomData& atom);

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
    lastBondData = SimpleLastBondData::Unbonded;
  }

  //! @brief Parsing trigger on encountering non-default bond information
  inline void setNextAtomBondInformation(const BondData& bond) {
    lastBondData = bond;
  }
//!@}

  //! @brief Interpret the collected graph as (possibly multiple molecules)
  std::vector<Molecule> interpret();

private:
//!@name Static private members
//!@{
  //! Checks whether an element type is valence filled
  static bool isValenceFillElement(Utils::ElementType e);

  //! Determines the implicit hydrogen count of a valence fill element
  static unsigned valenceFillElementImplicitHydrogenCount(
    int valence,
    Utils::ElementType e,
    bool aromatic
  );

  /*! Determines the mutual bond type of two bond type optionals
   *
   * @throws std::runtime_error If the bond types are mismatched
   */
  static BondType mutualBondType(
    const boost::optional<BondType>& a,
    const boost::optional<BondType>& b
  );

  //! Fetches a map to help with the atom chiral markers
  static std::vector<Shapes::Vertex> shapeMap(const ChiralData& chiralData);
//!@}

//!@name Private member functions
//!@{
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
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap
  );

  //! Set bond stereo post-parse and conversion to molecules
  void setBondStereo(
    std::vector<Molecule>& molecules,
    const std::vector<unsigned>& componentMap,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap
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
  enum class SimpleLastBondData {
    Unbonded,
    Unspecified
  };

  //! State for last stored bond data
  boost::variant<SimpleLastBondData, BondData> lastBondData = SimpleLastBondData::Unbonded;

  //! Possibly disconnected tracking graph
  PrivateGraph graph;

  //! State to track the vertex a new vertex is bound to
  std::stack<PrivateGraph::Vertex> vertexStack;

  //! Storage for bonds marked with stereo indicators ("/" and "\")
  using StereoMarkedBondTuple = std::tuple<PrivateGraph::Vertex, PrivateGraph::Vertex, BondData::StereoMarker>;
  std::vector<StereoMarkedBondTuple> stereoMarkedBonds;

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
