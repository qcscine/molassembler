#include "SymmetryInformation.h"

namespace PermSymmetry {

constexpr unsigned Linear::size;
constexpr unsigned TrigonalPlanar::size;
constexpr unsigned Tetrahedral::size;
constexpr unsigned SquarePlanar::size;
constexpr unsigned SquarePyramidal::size;
constexpr unsigned TrigonalBiPyramidal::size;
constexpr unsigned Octahedral::size;

constexpr std::array<
  std::array<unsigned, 2>,
  1
> Linear::rotations;

constexpr std::array<
  std::array<unsigned, 3>,
  2
> TrigonalPlanar::rotations;

constexpr std::array<
  std::array<unsigned, 4>,
  4
> Tetrahedral::rotations;

constexpr std::array<
  std::array<unsigned, 4>,
  3
> SquarePlanar::rotations;

constexpr std::array<
  std::array<unsigned, 5>,
  1
> SquarePyramidal::rotations;

constexpr std::array<
  std::array<unsigned, 5>,
  4
> TrigonalBiPyramidal::rotations;

constexpr std::array<
  std::array<unsigned, 6>,
  3
> Octahedral::rotations;

}
