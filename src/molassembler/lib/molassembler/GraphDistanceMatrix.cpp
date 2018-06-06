#include "GraphDistanceMatrix.h"

#include "detail/StdlibTypeAlgorithms.h"

namespace molassembler {

GraphDistanceMatrix::GraphDistanceMatrix(const AdjacencyMatrix& adjacencyMatrix)
: N(adjacencyMatrix.N) {
  _matrix.resize(N, N);

  // copy in from adjacencyMatrix
  for(AtomIndexType i = 0; i < N; i++) {
    for(AtomIndexType j = i + 1; j < N; j++) {
      this->operator()(i, j) = static_cast<unsigned>(
        adjacencyMatrix(i, j) // is boolean
      );
    }
  }

  _transformToDistances();
}

void GraphDistanceMatrix::_copyInRow(
  const AtomIndexType& sourceRow,
  const AtomIndexType& targetRow
) {
  for(unsigned col = targetRow + 1; col < N; col++) {
    if(
      col != sourceRow
      && this->operator()(sourceRow, col) != 0
    ) {
      // C++17 if-init opportunity!
      auto replacementDistance = (
        this->operator()(sourceRow, col)
        + this->operator()(targetRow, sourceRow)
      );
      if(
        this->operator()(targetRow, col) == 0
        || this->operator()(targetRow, col) > replacementDistance
      ) {
        this->operator()(targetRow, col) = replacementDistance;
      }
    }
  }
}

void GraphDistanceMatrix::_transformToDistances() {
  for(unsigned i = 0; i < N; i++) { // for every row
    // store whether the rows below have been added in
    // addedIn[j] says whether the row i+j+1 has been added in
    std::vector<bool> addedIn (N, false);
    addedIn[i] = true;

    for(
      unsigned counter = 1;
      !temple::all_of(addedIn);
      counter++
    ) {

      // Find rows that have value equal to current counter
      std::vector<unsigned> rowsToCopy;
      for(unsigned j = 0; j < N; j++) {
        if(
          i != j
          && !addedIn[j]
          && this->operator()(i, j) == counter
        ) {
          rowsToCopy.push_back(j);
        }
      }

      if(rowsToCopy.empty()) {
        break;
      }

      // copy them in
      for(const auto& copyIndex : rowsToCopy) {
        _copyInRow(i, copyIndex);
        addedIn[copyIndex] = true;
      }
    }
  }
}

/*!
 * Extract shortest chains of indices from a to b, e.g. in
 *
 *     1
 *    / \
 *   0   3 – 4
 *    \ /
 *     2
 *
 * if a = 0 and b = 4, return vector{vector{0, 1, 3, 4}, vector{0, 2, 3, 4}}.
 *
 * Main Idea behind algorithm: Keep a vector of chains. If, from the last
 * position in the chain, there are multiple adjacent atoms that have lower
 * distance to the target atom, create a chain for each in every iteration.
 */
std::vector<
  std::vector<AtomIndexType>
> GraphDistanceMatrix::extractChains(
  const AtomIndexType& a,
  const AtomIndexType& b
) {
  // return vectors
  std::vector<
    std::vector<AtomIndexType>
  > chains { {a} };

  // a shorthand
  auto distanceToTarget = [this, &b](const AtomIndexType& index) {
    return this->operator()(index, b);
  };

  /* as long as the chains are not at the target
   * (since all shortest chains must have same length, checking one chain is
   * sufficient, and there will always be one chain since the graph is connected)
   */
  while(distanceToTarget(chains[0].back()) > 0) {

    std::vector<
      std::vector<AtomIndexType>
    > iterationChains;

    for(const auto& chain: chains) {
      // find all neighbors of the last atom that have lower distance to target
      for(AtomIndexType i = 0; i < N; i++) {
        if(
          i != chain.back() // avoids (i, i) access
          && this->operator()(chain.back(), i) == 1 // i is a direct neighbor
          && distanceToTarget(i) < distanceToTarget(chain.back())
        ) {
          iterationChains.emplace_back(
            StdlibTypeAlgorithms::copyMerge(chain, {i})
          );
        }
      }
    }

    // avoid copy and clear, instead swap and clear
    std::swap(chains, iterationChains);
    iterationChains.clear();
  }

  return chains;
}

unsigned& GraphDistanceMatrix::operator () (
  const AtomIndexType& i,
  const AtomIndexType& j
) {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}

unsigned GraphDistanceMatrix::operator () (
  const AtomIndexType& i,
  const AtomIndexType& j
) const {
  return _matrix(
    std::min(i, j),
    std::max(i, j)
  );
}


} // namespace molassembler
