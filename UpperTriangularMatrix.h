#ifndef INCLUDE_CONSTEXPR_MAGIC_UPPER_TRIANGULAR_MATRIX_H
#define INCLUDE_CONSTEXPR_MAGIC_UPPER_TRIANGULAR_MATRIX_H

#include "Vector.h"

/*! @file
 * 
 * Provides a \c constexpr class that stores the data of an upper-triangular 
 * matrix via a std::array and provides two-index access.
 */

/* TODO
 * - current constructors of UpperTriangular force you to specify the size of
 *   the resulting matrix, but it could be inferred from the size of the 
 *   array used for construction (and checked for correctness!)
 */

namespace ConstexprMagic {

namespace UpperTriangularMatrixImpl {

namespace index_conversion {
  //!  Converts from (i, j) matrix indices to the linear k index for the array
  template<unsigned long matrixSize, typename UnsignedType>
  constexpr unsigned toSingleIndex(const UnsignedType& i, const UnsignedType& j) {
    return (
      matrixSize * (matrixSize - 1) / 2
      - (matrixSize - i) * ((matrixSize - i) - 1) / 2 
      + j - i - 1
    );
  }

  //! Converts from the linear array to (i, j) matrix indices
  template<unsigned long matrixSize, typename UnsignedType>
  constexpr std::pair<UnsignedType, UnsignedType> toDoubleIndex(const UnsignedType& k) {
    UnsignedType i {
      matrixSize - 2 
      - static_cast<unsigned>(
        Math::floor(
          Math::sqrt(
            0.0
            + 4 * matrixSize * (matrixSize - 1) 
            - 8 * k 
            - 7
          ) / 2.0 
          - 0.5
        )
      )
    };

    return {
      i,
      (
        k + i + 1 
        - matrixSize * (matrixSize - 1) / 2 
        + (matrixSize - i) * ((matrixSize - i) - 1) / 2
      )
    };
  }
} // namespace index_conversion

} // namespace UpperTriangularMatrixImpl

/*!
 * Stores the data of an upper-triangular matrix in a linear array in an 
 * all-constexpr fashion.
 */
template<typename ValueType, unsigned long matrixSize>
class UpperTriangularMatrix {
private:
/* State */
  std::array<ValueType, matrixSize * (matrixSize - 1) / 2> _data;

public:
/* Public information */
  const unsigned N = matrixSize;

/* Constructors */
  constexpr explicit UpperTriangularMatrix(
    const std::array<ValueType, matrixSize * (matrixSize - 1) / 2>& data
  ) : _data(data)
  {}

/* Modification */
  double& at(
    const unsigned& i,
    const unsigned& j
  ) {
    if(i >= j || i >= N || j >= N) {
      throw "Index out of bounds!";
    }

    return _data[
      UpperTriangularMatrixImpl::index_conversion::toSingleIndex<matrixSize>(i, j)
    ];
  }

/* Information */
  constexpr double at(
    const unsigned& i,
    const unsigned& j
  ) const {
    if(i >= j || i >= N || j >= N) {
      throw "Index out of bounds!";
    }

    return _data[
      UpperTriangularMatrixImpl::index_conversion::toSingleIndex<matrixSize>(i, j)
    ];
  }
};

//! Helper constructing function that deduces the required type signature
template<unsigned long matrixSize, typename ValueType>
constexpr UpperTriangularMatrix<ValueType, matrixSize> makeUpperTriangularMatrix(
  const std::array<ValueType, matrixSize * (matrixSize - 1) / 2>& data
) {
  return UpperTriangularMatrix<ValueType, matrixSize>(data);
}

} // namespace ConstexprMagic

#endif
