/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/PrivateGraph.h"
#include "boost/optional.hpp"

namespace Scine {
namespace Molassembler {

class Molecule;
class BondStereopermutator;

namespace IO {

/**
 * @brief Helper class to represent the stereo configuration of a double bond
 *   as indicated in a SMILES string
 */
struct SmilesBondStereo {
  boost::optional<PrivateGraph::Vertex> left;
  PrivateGraph::Vertex right;
  boost::optional<PrivateGraph::Vertex> upOfLeft;
  boost::optional<PrivateGraph::Vertex> downOfLeft;
  boost::optional<PrivateGraph::Vertex> upOfRight;
  boost::optional<PrivateGraph::Vertex> downOfRight;

  unsigned findAssignment(
    BondStereopermutator stereopermutator,
    const Molecule& mol,
    const std::vector<PrivateGraph::Vertex>& indexInComponentMap
  ) const;
};

} // namespace IO
} // namespace Molassembler
} // namespace Scine
