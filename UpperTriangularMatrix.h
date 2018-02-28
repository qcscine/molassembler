#ifndef INCLUDE_CONSTEXPR_MAGIC_UPPER_TRIANGULAR_MATRIX_H
#define INCLUDE_CONSTEXPR_MAGIC_UPPER_TRIANGULAR_MATRIX_H

#include "Vector.h"
#include "Array.h"
#include "FloatingPointComparison.h"

/*! @file
 * 
 * Provides a \c constexpr class that stores the data of an upper-triangular 
 * matrix via a std::array and provides two-index access.
 */

namespace constable {

namespace UpperTriangularMatrixImpl {

// Can be changed to std::array and dependency on Array removed with C++17 
template<typename T, size_t size>
using ArrayType = Array<T, size>;

namespace index_conversion {
  //!  Converts from (i, j) matrix indices to the linear k index for the array
  template<size_t N, typename UnsignedType>
  constexpr unsigned toSingleIndex(const UnsignedType& i, const UnsignedType& j) {
    return (
      N * (N - 1) / 2
      - (N - i) * ((N - i) - 1) / 2 
      + j - i - 1
    );
  }

  //! Converts from the linear array to (i, j) matrix indices
  template<size_t N, typename UnsignedType>
  constexpr std::pair<UnsignedType, UnsignedType> toDoubleIndex(const UnsignedType& k) {
    UnsignedType i {
      N - 2 
      - static_cast<unsigned>(
        Math::floor(
          Math::sqrt(
            0.0
            + 4 * N * (N - 1) 
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
        - N * (N - 1) / 2 
        + (N - i) * ((N - i) - 1) / 2
      )
    };
  }
} // namespace index_conversion

constexpr bool isValidMatrixSize(const size_t& dataSize) {
  auto root = Math::sqrt(1.0 + 8 * dataSize);

  return(
    floating::isCloseRelative(
      root,
      static_cast<double>(Math::floor(root)),
      1e-6
    )
  );
}

constexpr size_t getMatrixSize(const size_t& dataSize) {
  return (
    1 + Math::sqrt(1.0 + 8 * dataSize)
  ) / 2;
}

} // namespace UpperTriangularMatrixImpl

/*!
 * Stores the data of an upper-triangular matrix in a linear array in an 
 * all-constexpr fashion.
 */
template<typename ValueType, size_t dataSize>
class UpperTriangularMatrix {
private:
/* State */
  UpperTriangularMatrixImpl::ArrayType<
    ValueType,
    dataSize
  > _data;

public:
/* Public information */
  static constexpr unsigned N = UpperTriangularMatrixImpl::getMatrixSize(dataSize);

/* Constructors */
  constexpr UpperTriangularMatrix() : _data {} {}

  template<
    template<typename, size_t> class ArrayType
  > constexpr explicit UpperTriangularMatrix(
    const ArrayType<ValueType, dataSize>& data
  ) : _data(data)
  {
    static_assert(
      UpperTriangularMatrixImpl::isValidMatrixSize(dataSize) && dataSize > 0,
      "Passed data size is not consistent with an upper triangular matrix!"
    );
  }

/* Modification */
  constexpr ValueType& at(
    const unsigned& i,
    const unsigned& j
  ) {
    if(i >= j || i >= N || j >= N) {
      throw "Index out of bounds!";
    }

    return _data.at(
      UpperTriangularMatrixImpl::index_conversion::toSingleIndex<N>(i, j)
    );
  }

/* Information */
  constexpr const ValueType& at(
    const unsigned& i,
    const unsigned& j
  ) const {
    if(i >= j || i >= N || j >= N) {
      throw "Index out of bounds!";
    }

    return _data.at(
      UpperTriangularMatrixImpl::index_conversion::toSingleIndex<N>(i, j)
    );
  }

  constexpr const decltype(_data)& getData() const {
    return _data;
  }
};

//! Helper constructing function that deduces the required type signature
template<
  template<typename, size_t> class ArrayType,
  typename ValueType,
  size_t size
> constexpr UpperTriangularMatrix<ValueType, size> makeUpperTriangularMatrix(
  const ArrayType<ValueType, size>& data
) {
  return UpperTriangularMatrix<ValueType, size>(data);
}

} // namespace constable

#endif
