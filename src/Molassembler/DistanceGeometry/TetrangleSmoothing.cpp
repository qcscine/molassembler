/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/DistanceGeometry/TetrangleSmoothing.h"
#include "Molassembler/DistanceGeometry/ValueBounds.h"
#include "Molassembler/Temple/constexpr/Array.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/Invoke.h"

#include <Eigen/Dense>
#include <cfenv>

namespace Scine {
namespace Molassembler {
namespace DistanceGeometry {

struct LU {
  Eigen::Matrix4d L = Eigen::Matrix4d::Zero();
  Eigen::Matrix4d U = Eigen::Matrix4d::Zero();

  LU(
    const Eigen::MatrixXd& bounds,
    const std::array<unsigned, 4>& indices
  ) {
    for(unsigned i = 0; i < 4; ++i) {
      for(unsigned j = i + 1; j < 4; ++j) {
        unsigned a = indices[i];
        unsigned b = indices[j];
        if(a > b) {
          std::swap(a, b);
        }

        U(i, j) = bounds(a, b);
        U(j, i) = bounds(a, b);
        L(i, j) = bounds(b, a);
        L(j, i) = bounds(b, a);
      }
    }
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Named parameters to indicate when indices passed to lower and upper are ordered
struct OrderedIndicesTag {};
constexpr OrderedIndicesTag orderedIndicesTag;

constexpr unsigned pickMissingInZeroToThree(unsigned a, unsigned b, unsigned c) {
  if(
    !(a < 4 && b < 4 && c < 4)
    || a == b
    || a == c
    || b == c
  ) {
    throw std::logic_error("Unsafe arguments!");
  }

  unsigned x = 0;

  for(unsigned i = 0; i < 4; ++i) {
    if(a != i && b != i && c != i) {
      x = i;
      break;
    }
  }

  return x;
}

struct TriCheckResult {
  /* Lower and upper bounds on the 3-4/k-l distance as per triangle
   * inequalities
   */
  double klLowerBound; // l3_lim
  double klUpperBound; // u3_lim

  /* Whether or not the various triangle limits can be attained without
   * violating any of the passed bounds.
   */
  bool l_col = false;
  bool u_col = false;
};

template<typename MatrixType, typename ... AssumeISmallerJPack>
std::conditional_t<
  std::is_const<std::remove_reference_t<MatrixType>>::value,
  double,
  double&
> upper(MatrixType&& matrix, const unsigned i, const unsigned j, AssumeISmallerJPack ... pack) {
  if(sizeof...(pack) > 0) {
    return matrix(i, j);
  }

  if(i < j) {
    return matrix(i, j);
  }

  return matrix(j, i);
}

template<typename MatrixType, typename ... AssumeISmallerJPack>
std::conditional_t<
  std::is_const<std::remove_reference_t<MatrixType>>::value,
  double,
  double&
> lower(MatrixType&& matrix, const unsigned i, const unsigned j, AssumeISmallerJPack ... pack) {
  if(sizeof...(pack) > 0) {
    return matrix(j, i);
  }

  if(i < j) {
    return matrix(j, i);
  }

  return matrix(i, j);
}

double T(
  const double d_pr,
  const double d_qr,
  const double d_ps,
  const double d_qs
) {
  if(d_pr + d_qr > d_ps + d_qs || d_pr + d_qr < std::fabs(d_ps - d_qs)) {
    throw std::logic_error("T preconditions violated!");
  }

  return (
    d_pr * (std::pow(d_qs, 2) - std::pow(d_qr, 2))
    + d_qr * (std::pow(d_ps, 2) - std::pow(d_pr, 2))
  ) / (d_pr + d_qr);
}

bool collinear(
  const double d_pr,
  const double d_qr,
  const double l_ps,
  const double l_qs,
  const double l_rs,
  const double u_ps,
  const double u_qs,
  const double u_rs
) {
  auto U = [&]() -> double {
    if(d_pr + d_qr < u_ps - u_qs) {
      return std::pow(d_qr + u_qs, 2);
    }

    if(d_pr + d_qr < u_qs - u_ps) {
      return std::pow(d_pr + u_ps, 2);
    }

    return T(d_pr, d_qr, u_ps, u_qs);
  };

  auto L = [&]() -> double {
    if(d_pr + d_qr > l_ps + l_qs) {
      return std::pow(
        std::max({
          0.0,
          d_qr - u_qs,
          d_pr - u_ps,
          l_qs - d_qr,
          l_ps - d_pr
        }),
        2
      );
    }

    if(d_pr + d_qr < l_ps - l_qs) {
      return std::pow(l_ps - d_pr, 2);
    }

    if(d_pr + d_qr < l_qs - l_ps) {
      return std::pow(l_qs - d_qr, 2);
    }

    return T(d_pr, d_qr, l_ps, l_qs);
  };

  return (
    d_pr + d_qr <= u_ps + u_qs
    && u_rs >= L()
    && l_rs <= U()
  );
}

[[gnu::pure]] constexpr unsigned decrement(unsigned i) {
  if(i == 0) {
    throw std::logic_error("Underflow");
  }

  return i - 1;
}

/* limitTest overload set definitions section */

template<bool isUpper, unsigned i, unsigned j, unsigned k>
std::enable_if_t<isUpper, bool> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double minimalUpperLimit,
  const double tripletLimit
) {
  // Triangle upper limit test
  static_assert(i == 2 && k == 3, "Unexpected instantiation!");
  constexpr unsigned l = pickMissingInZeroToThree(i, j, k);
  static_assert(j < 4 && l < 4, "Zero-based indices");

  return (
    tripletLimit == minimalUpperLimit
    && collinear(
      upper(i, decrement(3)),
      upper(i, decrement(4)),
      lower(j, decrement(3)),
      lower(j, decrement(4)),
      lower(i, j),
      upper(j, decrement(3)),
      upper(j, decrement(4)),
      upper(i, j)
    )
  );
}

template<bool isUpper, unsigned i, unsigned k, unsigned l>
std::enable_if_t<!isUpper, bool> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double tripletLimit
) {
  // Triangle lower limit test
  constexpr unsigned j = pickMissingInZeroToThree(i, k, l);
  return(
    tripletLimit == maximalLowerLimit
    && collinear(
      upper(i, k),
      /* Paper has l_3(b_k, b_l) below, assuming that's a typo, as it would
       * mean the lower triangle inequality bound between k and l instead
       * of the current lower bound on that quantity. There is no such
       * subscripted quantity in the place of the paper where triangle upper
       * limit tests are performed.
       */
      lower(k, l),
      lower(i, j),
      lower(j, l),
      lower(j, k),
      upper(i, j),
      upper(j, l),
      upper(j, k)
    )
  );
}

template<bool isUpper, unsigned k, unsigned i, unsigned j, unsigned l>
std::enable_if_t<isUpper, bool> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double minimalUpperLimit,
  const double quadrupletLimit
) {
  // Tetrangle upper limit test
  if(quadrupletLimit == minimalUpperLimit) {
    const double firstConditionA = upper(j, k);
    const double firstConditionB = (
      upper(i, k) + upper(i, j)
    );
    const double firstConditionC = lower(j, k);

    // ujk >= uik + uij >= ljk
    if(firstConditionA >= firstConditionB && firstConditionB >= firstConditionC) {
      const double secondConditionA = upper(i, l);
      const double secondConditionB = (
        upper(j, l) + upper(i, l)
      );
      const double secondConditionC = lower(i, l);

      // uil >= ujl + uil >= lil
      return (secondConditionA >= secondConditionB && secondConditionB >= secondConditionC);
    }

    return false;
  }

  return false;
}

// First set of 4-atom lower limits with L[i, k, l, j] with {k, l} = {3, 4}
template<bool isUpper, unsigned i, unsigned k, unsigned l, unsigned j>
std::enable_if_t<
  !isUpper && (
    (k == 2 && l == 3) || (k == 3 && l == 2)
  ),
  bool
> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double quadrupletLimit
) {
  if(quadrupletLimit == maximalLowerLimit) {
    const double firstConditionA = upper(i, k);
    const double firstConditionB = lower(i, j) - upper(i, k);
    const double firstConditionC = lower(j, k);

    // uik >= lij - uik >= ljk
    if(firstConditionA >= firstConditionB && firstConditionB >= firstConditionC) {
      const double secondConditionA = upper(i, l);
      const double secondConditionB = lower(i, j) - upper(j, l);
      const double secondConditionC = lower(i, l);

      // uil >= lij - ujl >= lil
      return (secondConditionA >= secondConditionB && secondConditionB >= secondConditionC);
    }

    return false;
  }

  return false;
}

/* Second set of 4-atom lower limits with L[i, j, k, l] with {k, l} = {3, 4}
 * Here, k and l are the last two indices instead of in the middle.
 */
template<bool isUpper, unsigned i, unsigned j, unsigned k, unsigned l>
std::enable_if_t<
  !isUpper && (
    (k == 2 && l == 3) || (k == 3 && l == 2)
  ),
  bool
> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double quadrupletLimit
) {
  if(quadrupletLimit == maximalLowerLimit) {
    const double firstConditionA = upper(j, l);
    const double firstConditionB = lower(i, l) - upper(i, j);
    const double firstConditionC = lower(j, l);

    // ujl >= lil - uij >= ljl
    if(firstConditionA >= firstConditionB && firstConditionB >= firstConditionC) {
      const double secondConditionA = upper(i, k);
      const double secondConditionB = upper(i, j) + upper(j, k);
      const double secondConditionC = lower(i, k);

      // uik >= uij + ujk >= lik
      return (secondConditionA >= secondConditionB && secondConditionB >= secondConditionC);
    }

    return false;
  }

  return false;
}

bool triangleInequalitiesHold(const Eigen::Matrix3d& matrix) {
  // Ensure triangle inequalities are satisified in this matrix
  for(unsigned k = 0; k < 3; ++k) {
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = i + 1; j < 3; ++j) {
        if(upper(matrix, i, j, orderedIndicesTag) > upper(matrix, i, k) + upper(matrix, k, j)) {
          return false;
        }

        if(lower(matrix, i, j, orderedIndicesTag) < lower(matrix, i, k) - upper(matrix, k, j)) {
          return false;
        }
      }
    }
  }

  return true;
}

bool zeroBoundTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit
) {
  auto l = [&lower](unsigned i, unsigned j) {
    return lower(decrement(i), decrement(j));
  };
  auto u = [&upper](unsigned i, unsigned j) {
    return upper(decrement(i), decrement(j));
  };

  if(maximalLowerLimit == 0.0) {
    Eigen::Matrix3d c;

    // if i < j, then i,j an upper bound, i > j, a lower bound
    c(1, 0) = l(1, 2);
    c(0, 1) = u(1, 2);
    c(2, 0) = std::max(l(1, 3), l(1, 4));
    c(0, 2) = std::min(u(1, 3), u(1, 4));
    c(2, 1) = std::max(l(2, 3), l(2, 4));
    c(1, 2) = std::min(u(2, 3), u(2, 4));

    return triangleInequalitiesHold(c);
  }

  return false;
}

/* This function is a catch-all for the various bounds calculated in the
 * beginning of triCheck based on an index_sequence type argument.
 *
 * Although it may seem ridiculous to instantiate this for each index sequence
 * if the matrix element access is resolved at runtime anyway, using
 * index_sequence to store the various sequences allows for code reuse in the
 * checking stage of triCheck. It also allows us to group the expressions
 * for the various bounds nicely.
 *
 * Note that this takes zero-based indices! Make sure to use decrement() for
 * the index_sequence with which this is instantiated.
 */
constexpr unsigned NAngleBoundNoneValue = 100;
template<bool isUpper, bool firstPattern, unsigned i, unsigned j, unsigned k, unsigned l = NAngleBoundNoneValue>
double NAngleBound(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper) {
  static_assert(i < 4 && j < 4 && k < 4, "Zero-based indices!");
  static_assert(l == NAngleBoundNoneValue || l < 4, "Optional d argument is not zero-based!");

  // NOTE all ifs here could be if constexpr in C++17
  if(l == NAngleBoundNoneValue) {
    // Triangle algorithms
    if(isUpper) {
      // Upper triangle bound
      return upper(i, j) + upper(j, k);
    }

    // Lower triangle bound
    return lower(i, k) - upper(i, j);
  }

  // Tetrangle algorithms: l != NAngleBoundNoneValue
  if(isUpper) {
    // Upper tetrangle bound
    return upper(i, j) + upper(j, k) + upper(k, l);
  }

  // Lower tetrangle bound (THESE ARE IRREGULAR!)
  if(firstPattern) {
    return lower(i, l) - upper(i, j) - upper(k, l); // for the first two
  }

  return lower(i, l) - upper(i, j) - upper(j, k); // for the last four
}

template<bool isUpper, bool firstPattern, std::size_t ... Inds>
double NAngleBoundsMap(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper, std::index_sequence<Inds ...> /* inds */) {
  // Raise the index sequence from function argument to template argument
  return NAngleBound<isUpper, firstPattern, decrement(Inds)...>(lower, upper);
}

template<typename TupleType, bool isUpper, bool firstPattern, std::size_t ... TupleInds>
auto mapIndexSetsHelper(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper, std::index_sequence<TupleInds ...> /* inds */) {
  // Create an array with the results of bounds evaluations
  return std::array<double, sizeof...(TupleInds)> {
    NAngleBoundsMap<isUpper, firstPattern>(
      lower,
      upper,
      std::tuple_element_t<TupleInds, TupleType> {}
    )...
  };
}

/* isUpper denotes whether we are calculating an upper or lower bound
 * firstPattern is only for the tetrangle lower limits, where there are two
 * index patterns.
 */
template<typename TupleType, bool isUpper, bool firstPattern>
auto mapIndexSets(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper) {
  // Enumerate the index sequences in TupleType
  return mapIndexSetsHelper<TupleType, isUpper, firstPattern>(
    lower,
    upper,
    std::make_index_sequence<std::tuple_size<TupleType>::value> {}
  );
}

template<typename ArgsTuple, typename CallableArgPair>
bool foldCallableArgsPairLogicalOr(
  const ArgsTuple& args,
  CallableArgPair&& a
) {
  /* Piece together the arguments in order to call LimitTester::operator(),
   * which forwards its arguments to the limitTest overload set
   *
   * The CallableArgPair is created in limitAnyOfHelper and consists of a
   * LimitTester and a part of the arguments needed to call it.
   */
  return Temple::invoke(a.first, std::tuple_cat(args, std::tie(a.second)));
}

// C++17 fold
template<typename ArgsTuple, typename CallableArgPair, typename OtherCallableArgPair, typename ... Conditionals>
bool foldCallableArgsPairLogicalOr(
  const ArgsTuple& args,
  CallableArgPair&& a,
  OtherCallableArgPair&& b,
  Conditionals ... conditionals
) {
  // Fold the results of calls hopefully with short-circuiting
  return (
    foldCallableArgsPairLogicalOr(args, a)
    || foldCallableArgsPairLogicalOr(args, b, conditionals...)
  );
}

template<bool isUpper, std::size_t ... Inds>
struct LimitTester {
  template<typename ... Args>
  bool operator() (Args&& ... args) {
    return limitTest<isUpper, decrement(Inds)...>(std::forward<Args>(args)...);
  }
};

template<bool isUpper, std::size_t ... Inds>
auto makeLimitTester(std::index_sequence<Inds ...> /* inds */) {
  // Raise index sequence from function argument to template argument
  return LimitTester<isUpper, Inds ...> {};
}

template<typename TupleType, bool isUpper, typename ArgsTuple, std::size_t ... Inds>
auto limitAnyOfHelper(
  const ArgsTuple& args,
  const std::array<double, std::tuple_size<TupleType>::value>& values,
  std::index_sequence<Inds ...> /* inds */
) {
  /* Create a long list of pairs of test functors corresponding to index
   * sequences in TupleType that we can fold with logical or
   */
  return foldCallableArgsPairLogicalOr(
    args,
    std::make_pair(
      makeLimitTester<isUpper>(std::tuple_element_t<Inds, TupleType> {}),
      values.at(Inds)
    )...
  );
}

template<typename TupleType, bool isUpper, typename ArgsTuple>
auto limitAnyOf(
  const ArgsTuple& args,
  const std::array<double, std::tuple_size<TupleType>::value>& values
) {
  // Enumerate the indices into TupleType
  return limitAnyOfHelper<TupleType, isUpper>(
    args,
    values,
    std::make_index_sequence<std::tuple_size<TupleType>::value> {}
  );
}

TriCheckResult triCheck(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper
) {
  TriCheckResult result;
  /* The comments here sometimes refer to the parts of TRI_CHECK that they
   * implement from the original algorithm.
   *
   * Here we first define all the index sequences we will need as types so we
   * can freely mess with them. Note that these are one-based purely to match
   * the definitions from the original algorithm description.
   */
  // Definitions of upper triangle index sets that `U` is called with
  using UpperTriangleIndexSets = std::tuple<
    std::index_sequence<3, 1, 4>,
    std::index_sequence<3, 2, 4>
  >;

  // Definitions of upper tetrangle index sets that `U` is called with
  using UpperTetrangleIndexSets = std::tuple<
    std::index_sequence<3, 1, 2, 4>,
    std::index_sequence<4, 1, 2, 3>
  >;

  // Definitions of lower triangle index sets that `L` is called with
  using LowerTriangleIndexSets = std::tuple<
    std::index_sequence<1, 3, 4>,
    std::index_sequence<2, 3, 4>,
    std::index_sequence<1, 4, 3>,
    std::index_sequence<2, 4, 3>
  >;

  // Lower tetrangle index sets following the pattern l_ik - u_ij - u_kl
  using FirstLowerTetrangleIndexSets = std::tuple<
    std::index_sequence<1, 3, 4, 2>,
    std::index_sequence<1, 4, 3, 2>
  >;

  // Lower tetrangle index sets following the pattern l_ik - u_ij - u_jk
  using SecondLowerTetrangleIndexSets = std::tuple<
    std::index_sequence<1, 2, 3, 4>,
    std::index_sequence<2, 1, 3, 4>,
    std::index_sequence<1, 2, 4, 3>,
    std::index_sequence<2, 1, 4, 3>
  >;

  /* Now we actually calculate all of those values by calling `U` and `L` with
   * them. mapIndexSets takes care of forwarding the index sequences to the
   * right functions and grouping the results into arrays. The first boolean
   * template argument denotes whether we are calculating an upper bound, which
   * we are.
   *
   * The second is relevant only to the lower tetrangle bounds, and its value
   * is just set false here.
   *
   * The three- and four-argument `L` and `U` functions are represented by
   * the function NAngleBound.
   */
  const auto upperTriangleBounds = mapIndexSets<UpperTriangleIndexSets, true, false>(lower, upper);
  const auto upperTetrangleBounds = mapIndexSets<UpperTetrangleIndexSets, true, false>(lower, upper);

  /* Now for the lower bounds. Regarding the tetrangle bounds, you saw in the
   * index sequence definition that there are two sets, each following a
   * different pattern. If the second boolean template argument to mapIndexSets
   * is true, that means the current set follows the first pattern.
   */
  const auto lowerTriangleBounds = mapIndexSets<LowerTriangleIndexSets, false, false>(lower, upper);
  const auto firstLowerTetrangleBounds = mapIndexSets<FirstLowerTetrangleIndexSets, false, true>(lower, upper);
  const auto secondLowerTetrangleBounds = mapIndexSets<SecondLowerTetrangleIndexSets, false, false>(lower, upper);

  /* We set u3 and l3 from the minimum and maximum of the calculated bounds */
  result.klUpperBound = std::min(
    Temple::min(upperTriangleBounds), // min(U[i, j, k])
    Temple::min(upperTetrangleBounds) // min(U[i, j, k, l])
  );

  result.klLowerBound = std::max({
    0.0, // L[0] := 0
    Temple::max(lowerTriangleBounds), // min(L[i, j, k])
    Temple::max(firstLowerTetrangleBounds), // min(L[i, j, k, l]) of first pattern
    Temple::max(secondLowerTetrangleBounds) // min(L[i, j, k, l]) of second pattern
  });

  /* Now we determine u_col and l_col. These are "true whenever the triangle
   * inequality limits on the (3, 4)-distance are attainable without violating
   * the given bounds or the tetrangle inequality.
   *
   * The algorithm itself has multiple "For each 3-atom upper limit", "For each
   * 4-atom upper limit" parts that we abstract over here. If any of those
   * trigger, then the rest need not be evaluated, so we short-circuit with
   * logical ors.
   */

  const auto upperArgs = std::tie(lower, upper, result.klUpperBound);
  result.u_col = (
    limitAnyOf<UpperTriangleIndexSets, true>(upperArgs, upperTriangleBounds)
    || limitAnyOf<UpperTetrangleIndexSets, true>(upperArgs, upperTetrangleBounds)
  );

  const auto lowerArgs = std::tie(lower, upper, result.klLowerBound);
  result.l_col = (
    Temple::invoke(zeroBoundTest, lowerArgs)
    || limitAnyOf<LowerTriangleIndexSets, false>(lowerArgs, lowerTriangleBounds)
    || limitAnyOf<FirstLowerTetrangleIndexSets, false>(lowerArgs, firstLowerTetrangleBounds)
    || limitAnyOf<SecondLowerTetrangleIndexSets, false>(lowerArgs, secondLowerTetrangleBounds)
  );

  return result;
}

double CMUpper(
  const double d12,
  const double d13,
  const double d14,
  const double d23,
  const double d24
) {
  const double a = std::pow(d12, 2);
  const double b = std::pow(d13, 2);
  const double c = std::pow(d14, 2);
  const double d = std::pow(d23, 2);
  const double e = std::pow(d24, 2);

  const double A = -a / 4;
  const double B = (
    a * (-a + b + c + d + e)
    - b * c + b * e + c * d - d * e
  ) / 4;
  const double C = (
    d * (- a * b + a * c + b * c - c * c - c * d)
    + e * (a * b - b * b - a * c + b * c + b * d + c * d - b * e)
  ) / 4;

  const double discriminant = B * B + a * C; // A = -a / 4 -> -4 A = a

  if(discriminant >= 0) {
    return (-B - std::sqrt(discriminant)) / (2.0 * A);
  }

  return std::numeric_limits<double>::lowest();
}

double CMLower(
  const double d12,
  const double d13,
  const double d14,
  const double d23,
  const double d24
) {
  // NOTE: Nearly identical to CMUpper
  const double a = std::pow(d12, 2);
  const double b = std::pow(d13, 2);
  const double c = std::pow(d14, 2);
  const double d = std::pow(d23, 2);
  const double e = std::pow(d24, 2);

  const double A = -a / 4;
  const double B = (
    a * (-a + b + c + d + e)
    - b * c + b * e + c * d - d * e
  ) / 4;
  const double C = (
    d * (- a * b + a * c + b * c - c * c - c * d)
    + e * (a * b - b * b - a * c + b * c + b * d + c * d - b * e)
  ) / 4;

  const double discriminant = B * B + a * C; // A = -a / 4 -> -4 A = a

  if(discriminant >= 0) {
    return (-B + std::sqrt(discriminant)) / (2.0 * A);
  }

  return std::numeric_limits<double>::max();
}

/* Namespace with non-templated lower and upper functions so they can be
 * template parameters for other functions
 */
namespace m {

template<typename FnTuple, std::size_t ... Inds>
auto makeTetrangleCallTupleHelper(
  const Eigen::MatrixXd& bounds,
  const std::array<unsigned, 4>& b,
  FnTuple fnTuple,
  std::index_sequence<Inds...> /* inds */
) {
  constexpr auto indexTuple = std::make_tuple(
    std::pair<unsigned, unsigned> {0, 1},
    std::pair<unsigned, unsigned> {0, 2},
    std::pair<unsigned, unsigned> {0, 3},
    std::pair<unsigned, unsigned> {1, 2},
    std::pair<unsigned, unsigned> {1, 3}
  );

  return std::make_tuple(
    std::get<Inds>(fnTuple)(
      bounds,
      b[std::get<Inds>(indexTuple).first],
      b[std::get<Inds>(indexTuple).second]
    )...
  );
}

template<typename ... Fns>
auto fetchAppropriateBounds(
  const Eigen::MatrixXd& bounds,
  const std::array<unsigned, 4>& b,
  Fns ... fns
) {
  return makeTetrangleCallTupleHelper(
    bounds,
    b,
    std::forward_as_tuple(fns...),
    std::make_index_sequence<sizeof...(fns)>()
  );
}


double l(const Eigen::MatrixXd& bounds, const unsigned i, const unsigned j) {
  if(i < j) {
    return bounds(j, i);
  }

  return bounds(i, j);
}

double u(const Eigen::MatrixXd& bounds, const unsigned i, const unsigned j) {
  if(i < j) {
    return bounds(i, j);
  }

  return bounds(j, i);
}

} // namespace m

double upperTetrangleLimit(
  const Eigen::MatrixXd& bounds,
  const std::array<unsigned, 4>& b
) {
  /* Indices are always the following sequence: 12, 13, 14, 23, 24.
   * The only thing different among the tetrangle limits is whether we pass
   * the lower or upper bound. To make this readable, we pass merely functions
   * fetching the lower or upper bound from bounds and let
   * fetchAppropriateBounds take care of subindexing the indices stored in b
   * appropriately for each passed function.
   *
   * m::l is lower
   * m::u is upper
   *
   * These functions had to be duplicated from the original lower and upper
   * functions and fully qualified since we cannot pass templated functions as
   * an overload set.
   */
  return std::max({
    Temple::invoke(CMUpper, m::fetchAppropriateBounds(bounds, b, m::l, m::u, m::u, m::u, m::u)),
    Temple::invoke(CMUpper, m::fetchAppropriateBounds(bounds, b, m::u, m::l, m::l, m::u, m::u)),
    Temple::invoke(CMUpper, m::fetchAppropriateBounds(bounds, b, m::u, m::u, m::u, m::l, m::l))
  });
}

double lowerTetrangleLimit(
  const Eigen::MatrixXd& bounds,
  const std::array<unsigned, 4>& b
) {
  // See upperTetrangleLimit to explain the matrix below
  return std::min({
    Temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::u, m::u, m::l, m::l, m::u)),
    Temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::u, m::l, m::u, m::u, m::l)),
    Temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::l, m::l, m::u, m::l, m::u)),
    Temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::l, m::u, m::l, m::u, m::l))
  });
}

struct TetrangleLimits {
  DistanceGeometry::ValueBounds klLimits;
  bool boundViolation;

  TetrangleLimits(
    const Eigen::MatrixXd& bounds,
    const std::array<unsigned, 4>& b
  ) {
    LU matrixPair {bounds, b};
    // Can try the original triCheck here too
    TriCheckResult check = triCheck(matrixPair.L, matrixPair.U);

    /* NOTE: upperTetrangleLimit and lowerTetrangleLimit could use the LU
     * matrices for better data locality if this ever needs optimization
     */
    if(check.u_col) {
      klLimits.upper = check.klUpperBound;
    } else {
      const double limit = upperTetrangleLimit(bounds, b);
      klLimits.upper = std::sqrt(limit);
    }

    if(check.l_col) {
      klLimits.lower = check.klLowerBound;
    } else {
      const double limit = lowerTetrangleLimit(bounds, b);
      klLimits.lower = std::sqrt(limit);
    }

    boundViolation = klLimits.upper < klLimits.lower;
  }
};

unsigned tetrangleSmooth(Eigen::Ref<Eigen::MatrixXd> bounds) {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  const unsigned N = bounds.cols();

  // Minimal change in the bounds required to consider something has changed
  constexpr double epsilon = 0.01;

  bool changedSomething;
  unsigned iterations = 0;
  do {
    changedSomething = false;

    for(unsigned i = 0; i < N - 1; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        // NOTE (i,j) and (k,l) are not mutually disjoint
        for(unsigned k = 0; k < N - 1; ++k) {
          for(unsigned l = k + 1; l < N; ++l) {
            // Equal index pairs are trouble
            if(i == k && j == l) {
              continue;
            }

            const TetrangleLimits limits {bounds, {i, j, k, l}};

            if(limits.boundViolation) {
              throw std::runtime_error("Bound violation found!");
            }

            // k < l, so bounds(k, l) is the upper bound, bounds(l, k) the lower
            double& klLowerBound = bounds(l, k);
            double& klUpperBound = bounds(k, l);

            assert(klLowerBound <= klUpperBound);

            if(
              limits.klLimits.lower > klLowerBound
              && std::fabs(limits.klLimits.lower - klLowerBound) / klLowerBound > epsilon
            ) {
              if(limits.klLimits.lower > klUpperBound) {
                throw std::runtime_error("Bound violation found!");
              }

              klLowerBound = limits.klLimits.lower;
              changedSomething = true;
            }

            if(
              limits.klLimits.upper < klUpperBound
              && std::fabs(klUpperBound - limits.klLimits.upper) / klUpperBound > epsilon
            ) {
              if(limits.klLimits.upper < klLowerBound) {
                throw std::runtime_error("Bound violation found!");
              }
              klUpperBound = limits.klLimits.upper;
              changedSomething = true;
            }
          }
        }
      }
    }

    ++iterations;
  } while(changedSomething);

  fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  return iterations;
}

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine
