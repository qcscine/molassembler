#ifndef INCLUDE_ADJACENCY_MATRIX_H
#define INCLUDE_ADJACENCY_MATRIX_H

#include "Molecule.h"

/*! @file
 *
 * Contains a matrix class that represents molecular bonding in a boolean
 * manner.
 */

namespace molassembler {

/*!
 * A class constructed from a Molecule that captures bonding in a boolean
 * manner. No bond type information is present, merely if atoms are bonded or
 * not.
 */
class AdjacencyMatrix {
private:
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> _matrix;

public:
  const unsigned N;

  AdjacencyMatrix() = delete;
  explicit AdjacencyMatrix(const Molecule& molecule);

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& getMatrixRef() {
    return _matrix;
  }

  bool& operator () (
    const AtomIndexType& i,
    const AtomIndexType& j
  );

  bool operator () (
    const AtomIndexType& i,
    const AtomIndexType& j
  ) const;
};

}

#endif
