/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "molassembler/PRNG.h"

void init_random_engine(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<random::Engine> engine(m, "PRNG", "Pseudo-randomness source");
  engine.def(
    "seed",
    pybind11::overload_cast<int>(&random::Engine::seed),
    pybind11::arg("seed_number"),
    "Seed the PRNG with state"
  );
}
