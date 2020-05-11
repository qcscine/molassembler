/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "Molassembler/Prng.h"

void init_random_engine(pybind11::module& m) {
  using namespace Scine::Molassembler;

  pybind11::class_<Random::Engine> engine(
    m,
    "PRNG",
    R"delim(
      Pseudo-random number generator

      Central source of pseudo-randomness for the library.
    )delim"
  );

  engine.def(
    "seed",
    pybind11::overload_cast<int>(&Random::Engine::seed),
    pybind11::arg("seed_number"),
    "Seed the PRNG with state"
  );
}
