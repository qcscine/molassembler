/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/ContainerTraits.h"
#include "Molassembler/Temple/StrongIndex.h"
#include "Molassembler/Temple/StrongIndexPermutation.h"
#include <vector>

using namespace Scine::Molassembler;

BOOST_AUTO_TEST_CASE(PermutationsTests, *boost::unit_test::label("Temple")) {
  Temple::Permutation p {0, 1, 3, 2};
  BOOST_CHECK_EQUAL(p.index(), 1);

  // Applying it to a vector of other objects yields properly permuted elements
  const std::vector<int> nums {{10, -3, 4, -1}};
  const std::vector<int> permutedNums = p.apply(nums);
  const std::vector<int> expected {{10, -3, -1, 4}};
  BOOST_CHECK(permutedNums == expected);

  // Applying an ordering permutation yields an ordered vector
  const auto ordered = Temple::Permutation::ordering(nums).apply(nums);
  BOOST_CHECK(std::is_sorted(std::begin(ordered), std::end(ordered)));

  // Move forward and backward through permutations
  BOOST_CHECK(p.prev());
  BOOST_CHECK_EQUAL(p.index(), 0);
  BOOST_CHECK(p.next());
  BOOST_CHECK(p.next());
  BOOST_CHECK_EQUAL(p.index(), 2);

  // Composing with the inverse yields an identity permutation
  BOOST_CHECK_EQUAL(p.compose(p.inverse()).index(), 0);
  BOOST_CHECK_EQUAL(p.inverse().compose(p).index(), 0);

  auto q = Temple::Permutation::identity(4);
  BOOST_CHECK_EQUAL(q.index(), 0);
  BOOST_CHECK(q.next());
  BOOST_CHECK_EQUAL(q.index(), 1);

  BOOST_CHECK_EQUAL(q.compose(q.inverse()).index(), 0);

  Temple::Permutation r(4, 2);
  BOOST_CHECK_EQUAL(r.index(), 2);

  Temple::Permutation s(10, 1001);
  BOOST_CHECK_EQUAL(s.index(), 1001);

  // Composition order
  const Temple::Permutation t {2, 0, 3, 1};
  const Temple::Permutation u {2, 3, 1, 0};
  const Temple::Permutation tu {1, 2, 0, 3};
  const Temple::Permutation ut {3, 1, 0, 2};

  BOOST_CHECK(t.compose(u) == tu);
  BOOST_CHECK(u.compose(t) == ut);

  // Application
  const std::vector<int> forward {{-1, 4, 10, -3}};
  for(unsigned i = 0; i < u.size(); ++i) {
    BOOST_CHECK(forward.at(u(i)) == nums.at(i));
  }
  BOOST_CHECK(u.apply(nums) == forward);
  BOOST_CHECK(u.inverse().apply(forward) == nums);
}

BOOST_AUTO_TEST_CASE(StrongIndices, *boost::unit_test::label("Temple")) {
  struct FooTag;
  using Foo = Temple::StrongIndex<FooTag, unsigned>;

  struct BarTag;
  using Bar = Temple::StrongIndex<BarTag, unsigned>;

  Foo f(4);
  ++f;
  BOOST_CHECK(f == 5);

  Foo g(5);
  BOOST_CHECK(f == g);

  // This is an eyesore, possible because of non-explicit conversion
  Bar h(5);
  BOOST_CHECK(g == h);

  // This is even worse
  g = h;
}

template<class Permutation>
using TestCanIndexWithSizeT = decltype(std::declval<Permutation>().at(std::size_t {0}));

template<class Permutation>
struct CanIndexWithSizeT : std::integral_constant<bool, Temple::Traits::Detail::is_detected_v<TestCanIndexWithSizeT, Permutation>> {};

BOOST_AUTO_TEST_CASE(StrongPermutations, *boost::unit_test::label("Temple")) {
  struct FooTag;
  using Foo = Temple::StrongIndex<FooTag, unsigned>;

  struct BarTag;
  using Bar = Temple::StrongIndex<BarTag, unsigned>;

  using StrongPermutation = Temple::StrongIndexPermutation<Foo, Bar>;
  const StrongPermutation p {0, 2, 1};

  // Indexing with strong indices only
  static_assert(
    !CanIndexWithSizeT<StrongPermutation>::value,
    "Oops, can index StrongPermutation with weak key"
  );
  BOOST_CHECK(p.at(Foo {1}) == Bar {2});

  // Can invert to reverse template argument order
  const auto inverse = p.inverse();
  static_assert(
    std::is_same<std::decay_t<decltype(inverse)>, Temple::StrongIndexPermutation<Bar, Foo>>::value,
    "Bad type returned from inverse"
  );

  // Composition with strongly typed pre and post-permutations
  const Temple::StrongIndexPermutation<Bar, Bar> q {1, 0, 2};
  const Temple::StrongIndexPermutation<Foo, Bar> pq {1, 2, 0};
  BOOST_CHECK(p.compose(q) == pq);
  // auto x = q.compose(p); // doesn't compile!

  const Temple::StrongIndexPermutation<Foo, Foo> r {2, 1, 0};
  const Temple::StrongIndexPermutation<Foo, Bar> rp {1, 2, 0};
  BOOST_CHECK(r.compose(p) == rp);
  // auto x = p.compose(r);  // doesn't compile!
}
