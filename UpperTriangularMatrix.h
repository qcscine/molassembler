#ifndef INCLUDE_CONSTEXPR_MAGIC_UPPER_TRIANGULAR_MATRIX_H
#define INCLUDE_CONSTEXPR_MAGIC_UPPER_TRIANGULAR_MATRIX_H

/* TODO
 * - remove dependency
 * - current constructors of UpperTriangular force you to specify the size of
 *   the resulting matrix, but it could be inferred from the size of the 
 *   array used for construction (and checked for correctness!)
 */

#include "Vector.h"

namespace ConstexprMagic {

namespace UpperTriangularMatrixImpl {

namespace index_conversion {
  template<unsigned long size, typename UnsignedType>
  constexpr unsigned toSingleIndex(const UnsignedType& i, const UnsignedType& j) {
    return (
      size * (size - 1) / 2
      - (size - i) * ((size - i) - 1) / 2 
      + j - i - 1
    );
  }

  template<unsigned long size, typename UnsignedType>
  constexpr std::pair<UnsignedType, UnsignedType> toDoubleIndex(const UnsignedType& k) {
    UnsignedType i {
      size - 2 
      - static_cast<unsigned>(
        Math::floor(
          Math::sqrt(
            0.0
            + 4 * size * (size - 1) 
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
        - size * (size - 1) / 2 
        + (size - i) * ((size - i) - 1) / 2
      )
    };
  }
} // namespace index_conversion

} // namespace UpperTriangularMatrixImpl

template<typename ValueType, unsigned long size>
class UpperTriangularMatrix {
private:
/* State */
  std::array<ValueType, size * (size - 1) / 2> _data;

public:
/* Public information */
  const unsigned N = size;

/* Constructors */
  constexpr UpperTriangularMatrix(
    const std::array<ValueType, size * (size - 1) / 2>& data
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
      UpperTriangularMatrixImpl::index_conversion::toSingleIndex<size>(i, j)
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
      UpperTriangularMatrixImpl::index_conversion::toSingleIndex<size>(i, j)
    ];
  }
};

template<unsigned long size, typename ValueType>
constexpr UpperTriangularMatrix<ValueType, size> makeUpperTriangularMatrix(
  const std::array<ValueType, size * (size - 1) / 2>& data
) {
  return UpperTriangularMatrix<ValueType, size>(data);
}

} // namespace ConstexprMagic

#endif
