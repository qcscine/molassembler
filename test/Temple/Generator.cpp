/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/constexpr/Jsf.h"

using namespace Scine::Molassembler;

Temple::Generator<> generator;

BOOST_AUTO_TEST_CASE(JsfRandomness, *boost::unit_test::label("Temple")) {
  const int fixedSeed = 1042;

  Temple::JSF64 engine;
  engine.seed(fixedSeed);
  auto engineStateCopy = engine;
  for(unsigned i = 0; i < 10; ++i) {
    engine();
  }
  engine.seed(fixedSeed);
  BOOST_CHECK(engine == engineStateCopy);

  Temple::JSF64 directlySeededEngine {fixedSeed};
  BOOST_CHECK(engine == directlySeededEngine);

  std::seed_seq sequence {fixedSeed};
  directlySeededEngine = Temple::JSF64 {sequence};
  BOOST_CHECK(engine == directlySeededEngine);
}
