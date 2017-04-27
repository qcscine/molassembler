#include "Vector.h"
#include "UpperTriangularMatrix.h"

#include <iostream>

/* True client code ------------------------------- */

namespace make_constexpr_array {
  template<unsigned long size>
  constexpr double makeElement(
    const std::array<ConstexprMagic::Vector, size>& positions,
    const size_t& i
  ) {
    const auto indexPair = ConstexprMagic::UpperTriangularMatrixImpl::index_conversion::toDoubleIndex<size * (size - 1) /2>(i);

    return ConstexprMagic::toDegrees(
      ConstexprMagic::angle(
        positions[indexPair.first],
        positions[indexPair.second]
      )
    );
  }

  template<unsigned long size, size_t... Inds>
  constexpr std::array<double, size * (size - 1) / 2> makeArrayImpl(
    const std::array<ConstexprMagic::Vector, size>& positions,
    std::integer_sequence<size_t, Inds...>
  ) {
    return { makeElement(positions, Inds)... };
  }

  template<unsigned long size>
  constexpr std::array<double, size * (size - 1) / 2> makeArray(
    const std::array<ConstexprMagic::Vector, size>& positions
  ) {
    return makeArrayImpl(
      positions,
      std::make_index_sequence<size * (size - 1) / 2>{}
    );
  }

} // namespace make_constexpr_array


// Say we have three positions
constexpr std::array<
  ConstexprMagic::Vector,
  3
> positions {{
  {{-0.00928803, 0.61156848, 0.79113698}},
  {{0.79562737, 0.60564101, -0.01326839}},
  {{0.79562737, -0.60564101, -0.01326839}}
}};

constexpr auto matrixData = make_constexpr_array::makeArray<3>(positions);

constexpr auto angleMatrix = ConstexprMagic::makeUpperTriangularMatrix<3>(
  matrixData
);

int main() {
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = i + 1; j < 3; j++) {
      std::cout << "(" << i << ", " << j << ") = " << angleMatrix.at(i, j) << std::endl;
    }
  }

}
