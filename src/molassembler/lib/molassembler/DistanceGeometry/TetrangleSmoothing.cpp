#include "molassembler/DistanceGeometry/TetrangleSmoothing.h"
#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "temple/constexpr/Array.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Invoke.h"

#include <Eigen/Dense>

/* TODO
 * - test this with some examples to see if it blows up spectacularly or does
 *   its job on the first attempt. Try both variants of triCheck, too. The
 *   template trickery might mess with correctness in some places.
 */

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {

std::pair<
  Eigen::Matrix4d,
  Eigen::Matrix4d
> makeLU(
  const Eigen::MatrixXd& bounds,
  const std::array<unsigned, 4>& indices
) {
  auto result = std::make_pair<Eigen::Matrix4d, Eigen::Matrix4d>(
    Eigen::Matrix4d::Zero(),
    Eigen::Matrix4d::Zero()
  );

  Eigen::Matrix4d& lower = result.first;
  Eigen::Matrix4d& upper = result.first;

  for(unsigned i = 0; i < 4; ++i) {
    for(unsigned j = i + 1; j < 4; ++j) {
      unsigned a = indices[i];
      unsigned b = indices[j];
      if(b < a) {
        std::swap(a, b);
      }

      upper(i, j) = bounds(a, b);
      upper(j, i) = bounds(a, b);

      lower(i, j) = bounds(b, a);
      lower(j, i) = bounds(b, a);
    }
  }

  return result;
}

// Named parameters to indicate when indices passed to lower and upper are ordered
struct OrderedIndicesTag {};
constexpr OrderedIndicesTag orderedIndicesTag;

constexpr unsigned pickMissingInOneToFour(unsigned a, unsigned b, unsigned c) {
  if(
    !(1 <= a && a <= 4 && 1 <= b && b <= 4 && 1 <= c && c <= 4)
    || a == b
    || a == c
    || b == c
  ) {
    throw std::logic_error("Unsafe arguments!");
  }

  unsigned x = 0;

  for(unsigned i = 1; i < 5; ++i) {
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

  return (i - 1);
}

/* limitTest overload set definitions section */

template<bool isUpper, unsigned a, unsigned b, unsigned c>
std::enable_if_t<isUpper, bool> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double minimalUpperLimit,
  const double tripletLimit
) {
  // Triangle upper limit test
  static_assert(a == 3 && c == 4, "Unexpected instantiation!");
  constexpr unsigned d = pickMissingInOneToFour(a, b, c);
  static_assert(0 < b && b < 5 && 0 < d && d < 5, "One-based indices!");
  constexpr unsigned i = decrement(b);
  constexpr unsigned j = decrement(d);

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

template<bool isUpper, unsigned a, unsigned c, unsigned d>
std::enable_if_t<!isUpper, bool> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double tripletLimit
) {
  // Triangle lower limit test
  constexpr unsigned b = pickMissingInOneToFour(a, c, d);

  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);
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

template<bool isUpper, unsigned c, unsigned a, unsigned b, unsigned d>
std::enable_if_t<isUpper, bool> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double minimalUpperLimit,
  const double quadrupletLimit
) {
  // Tetrangle upper limit test
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);

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

// First set of 4-atom lower limits with L[a, c, d, b] with {c, d} = {3, 4}
template<bool isUpper, unsigned a, unsigned c, unsigned d, unsigned b>
std::enable_if_t<
  !isUpper && (
    (c == 3 && d == 4) || (c == 4 && d == 3)
  ),
  bool
> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double quadrupletLimit
) {
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);

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

/* Second set of 4-atom lower limits with L[a, b, c, d] with {c, d} = {3, 4}
 * Here, c and d are the last two indices instead of in the middle.
 */
template<bool isUpper, unsigned a, unsigned b, unsigned c, unsigned d>
std::enable_if_t<
  !isUpper && (
    (c == 3 && d == 4) || (c == 4 && d == 3)
  ),
  bool
> limitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double quadrupletLimit
) {
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);

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

/* NOTE: from here to MARKERA is needed only for the original triCheck impl */

template<unsigned a, unsigned b, unsigned c>
bool triangleUpperLimitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double minimalUpperLimit,
  const double tripletLimit
) {
  static_assert(a == 3 && c == 4, "Unexpected instantiation!");
  constexpr unsigned d = pickMissingInOneToFour(a, b, c);
  static_assert(0 < b && b < 5 && 0 < d && d < 5, "One-based indices!");
  constexpr unsigned i = decrement(b);
  constexpr unsigned j = decrement(d);

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

template<unsigned c, unsigned a, unsigned b, unsigned d>
bool tetrangleUpperLimitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double minimalUpperLimit,
  const double quadrupletLimit
) {
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);

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

template<unsigned a, unsigned c, unsigned d>
bool triangleLowerLimitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double triangleLimit
) {
  constexpr unsigned b = pickMissingInOneToFour(a, c, d);

  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);
  return(
    triangleLimit == maximalLowerLimit
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

// First set of 4-atom lower limits with L[a, c, d, b] with {c, d} = {3, 4}
template<unsigned a, unsigned c, unsigned d, unsigned b>
std::enable_if_t<
  (c == 3 && d == 4) || (c == 4 && d == 3),
  bool
> tetrangleLowerLimitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double quadrupletLimit
) {
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);

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

/* Second set of 4-atom lower limits with L[a, b, c, d] with {c, d} = {3, 4}
 * Here, c and d are the last two indices instead of in the middle.
 */
template<unsigned a, unsigned b, unsigned c, unsigned d>
std::enable_if_t<
  (c == 3 && d == 4) || (c == 4 && d == 3),
  bool
> tetrangleLowerLimitTest(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper,
  const double maximalLowerLimit,
  const double quadrupletLimit
) {
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);
  constexpr unsigned l = decrement(d);

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

template<unsigned a, unsigned b, unsigned c>
double upperNAngleBound(const Eigen::Matrix4d& upper) {
  static_assert(0 < a && a < 5 && 0 < b && b < 5 && 0 < c && c < 5, "One-based indices!");
  constexpr unsigned i = (a - 1);
  constexpr unsigned j = (b - 1);
  constexpr unsigned k = (c - 1);
  return upper(i, j) + upper(j, k);
}

template<unsigned a, unsigned b, unsigned c, unsigned d>
double upperNAngleBound(const Eigen::Matrix4d& upper) {
  static_assert(0 < a && a < 5 && 0 < b && b < 5 && 0 < c && c < 5, "One-based indices!");
  constexpr unsigned i = (a - 1);
  constexpr unsigned j = (b - 1);
  constexpr unsigned k = (c - 1);
  constexpr unsigned l = (d - 1);
  return upper(i, j) + upper(j, k) + upper(k, l);
}

template<unsigned a, unsigned b, unsigned c>
double lowerNAngleBound(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper) {
  static_assert(0 < a && a < 5 && 0 < b && b < 5 && 0 < c && c < 5, "One-based indices!");
  constexpr unsigned i = (a - 1);
  constexpr unsigned j = (b - 1);
  constexpr unsigned k = (c - 1);
  return lower(i, k) - upper(i, j);
}

template<unsigned a, unsigned b, unsigned c, unsigned d>
double lowerNAngleBound(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper) {
  static_assert(0 < a && a < 5 && 0 < b && b < 5 && 0 < c && c < 5 && 0 < d && d < 5, "One-based indices!");
  constexpr unsigned i = (a - 1);
  constexpr unsigned j = (b - 1);
  constexpr unsigned k = (c - 1);
  constexpr unsigned l = (d - 1);
  return lower(i, l) - upper(i, j) - upper(j, k);
}

TriCheckResult triCheck(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper
) {
  TriCheckResult result;

  /* Check to see if the triangle inequality and not the tetrangle inequality
   * is the limiting lower and/or upper constraint on the k-l distance.
   *
   * Determine whether or not the various triangle limits
   * (which are possible among four atoms)
   * can be attained without violating any of the given bounds.
   */
  const double U314 = upperNAngleBound<3, 1, 4>(upper);
  const double U324 = upperNAngleBound<3, 2, 4>(upper);
  const double U3124 = upperNAngleBound<3, 1, 2, 4>(upper);
  const double U4123 = upperNAngleBound<4, 1, 2, 3>(upper);

  result.klUpperBound = std::min({
    U314, U324,
    U3124, U4123
  });

  const double L134 = lowerNAngleBound<1, 3, 4>(lower, upper);
  const double L234 = lowerNAngleBound<2, 3, 4>(lower, upper);
  const double L143 = lowerNAngleBound<1, 4, 3>(lower, upper);
  const double L243 = lowerNAngleBound<2, 4, 3>(lower, upper);
  const double L1342 = lowerNAngleBound<1, 3, 4, 2>(lower, upper);
  const double L1432 = lowerNAngleBound<1, 4, 3, 2>(lower, upper);
  const double L1234 = lowerNAngleBound<1, 2, 3, 4>(lower, upper);
  const double L2134 = lowerNAngleBound<2, 1, 3, 4>(lower, upper);
  const double L1243 = lowerNAngleBound<1, 2, 4, 3>(lower, upper);
  const double L2143 = lowerNAngleBound<2, 1, 4, 3>(lower, upper);

  result.klLowerBound = std::max({
    0.0,
    L134, L234, L143, L243,
    L1342, L1432, L1234, L2134, L1243, L2143
  });

  result.u_col = (
    triangleUpperLimitTest<3, 1, 4>(lower, upper, result.klUpperBound, U314)
    || triangleUpperLimitTest<3, 2, 4>(lower, upper, result.klUpperBound, U324)
    || tetrangleUpperLimitTest<3, 1, 2, 4>(lower, upper, result.klUpperBound, U3124)
    || tetrangleUpperLimitTest<4, 1, 2, 3>(lower, upper, result.klUpperBound, U4123)
  );

  result.l_col = (
    zeroBoundTest(lower, upper, result.klLowerBound)
    || triangleLowerLimitTest<1, 3, 4>(lower, upper, result.klLowerBound, L134)
    || triangleLowerLimitTest<2, 3, 4>(lower, upper, result.klLowerBound, L234)
    || triangleLowerLimitTest<1, 4, 3>(lower, upper, result.klLowerBound, L143)
    || triangleLowerLimitTest<2, 4, 3>(lower, upper, result.klLowerBound, L243)
    || tetrangleLowerLimitTest<1, 3, 4, 2>(lower, upper, result.klLowerBound, L1342)
    || tetrangleLowerLimitTest<1, 4, 3, 2>(lower, upper, result.klLowerBound, L1432)
    || tetrangleLowerLimitTest<1, 2, 3, 4>(lower, upper, result.klLowerBound, L1234)
    || tetrangleLowerLimitTest<2, 1, 3, 4>(lower, upper, result.klLowerBound, L2134)
    || tetrangleLowerLimitTest<1, 2, 4, 3>(lower, upper, result.klLowerBound, L1243)
    || tetrangleLowerLimitTest<2, 1, 4, 3>(lower, upper, result.klLowerBound, L2143)
  );

  return result;
}

/* NOTE: MARKERA, below is relevant only to new implementation */

constexpr unsigned NAngleBoundNoneValue = 100;
template<bool isUpper, unsigned a, unsigned b, unsigned c, unsigned d = NAngleBoundNoneValue>
double NAngleBound(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper) {
  static_assert(0 < a && a < 5 && 0 < b && b < 5 && 0 < c && c < 5, "One-based indices!");
  static_assert(d == NAngleBoundNoneValue || (0 < d && d < 5), "Optional d argument is not one-based!");
  constexpr unsigned i = decrement(a);
  constexpr unsigned j = decrement(b);
  constexpr unsigned k = decrement(c);

  // NOTE all ifs here could be if constexpr in C++17

  if(d == NAngleBoundNoneValue) {
    // Triangle algorithms
    if(isUpper) {
      // Upper triangle bound
      return upper(i, j) + upper(j, k);
    }

    // Lower triangle bound
    return lower(i, k) - upper(i, j);
  }

  // Tetrangle algorithms: d != NAngleBoundNoneValue
  constexpr unsigned l = decrement(d);
  if(isUpper) {
    // Upper tetrangle bound
    return upper(i, j) + upper(j, k) + upper(k, l);
  }

  // Lower tetrangle bound
  return lower(i, l) - upper(i, j) - upper(j, k);
}

template<bool isUpper, std::size_t ... Inds>
double NAngleBoundsMap(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper, std::index_sequence<Inds ...> /* inds */) {
  return NAngleBound<isUpper, Inds...>(lower, upper);
}

template<typename TupleType, bool isUpper, std::size_t ... TupleInds>
auto mapIndexSetsHelper(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper, std::index_sequence<TupleInds ...> /* inds */) {
  return std::array<double, sizeof...(TupleInds)> {
    NAngleBoundsMap<isUpper>(
      lower,
      upper,
      std::tuple_element_t<TupleInds, TupleType> {}
    )...
  };
}

template<typename TupleType, bool isUpper>
auto mapIndexSets(const Eigen::Matrix4d& lower, const Eigen::Matrix4d& upper) {
  return mapIndexSetsHelper<TupleType, isUpper>(
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
  return temple::invoke(a.first, std::tuple_cat(args, std::tie(a.second)));
}

template<typename ArgsTuple, typename CallableArgPair, typename OtherCallableArgPair, typename ... Conditionals>
bool foldCallableArgsPairLogicalOr(
  const ArgsTuple& args,
  CallableArgPair&& a,
  OtherCallableArgPair&& b,
  Conditionals ... conditionals
) {
  return (
    foldCallableArgsPairLogicalOr(args, a)
    || foldCallableArgsPairLogicalOr(args, b, conditionals...)
  );
}

template<bool isUpper, std::size_t ... Inds>
struct LimitTester {
  template<typename ... Args>
  bool operator() (Args ... args) {
    return limitTest<isUpper, Inds ...>(args...);
  }
};

template<bool isUpper, std::size_t ... Inds>
auto makeLimitTester(std::index_sequence<Inds ...> /* inds */) {
  return LimitTester<isUpper, Inds ...> {};
}

template<typename TupleType, bool isUpper, typename ArgsTuple, std::size_t ... Inds>
auto limitAnyOfHelper(
  const ArgsTuple& args,
  const std::array<double, std::tuple_size<TupleType>::value>& values,
  std::index_sequence<Inds ...> /* inds */
) {
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
  return limitAnyOfHelper<TupleType, isUpper>(
    args,
    values,
    std::make_index_sequence<std::tuple_size<TupleType>::value> {}
  );
}

TriCheckResult triCheckRetry(
  const Eigen::Matrix4d& lower,
  const Eigen::Matrix4d& upper
) {
  TriCheckResult result;

  using UpperTriangleIndexSets = std::tuple<
    std::index_sequence<3, 1, 4>,
    std::index_sequence<3, 2, 4>
  >;

  using UpperTetrangleIndexSets = std::tuple<
    std::index_sequence<3, 1, 2, 4>,
    std::index_sequence<4, 1, 2, 3>
  >;

  using LowerTriangleIndexSets = std::tuple<
    std::index_sequence<1, 3, 4>,
    std::index_sequence<2, 3, 4>,
    std::index_sequence<1, 4, 3>,
    std::index_sequence<2, 4, 3>
  >;

  using LowerTetrangleIndexSets = std::tuple<
    std::index_sequence<1, 3, 4, 2>,
    std::index_sequence<1, 4, 3, 2>,
    std::index_sequence<1, 2, 3, 4>,
    std::index_sequence<2, 1, 3, 4>,
    std::index_sequence<1, 2, 4, 3>,
    std::index_sequence<2, 1, 4, 3>
  >;

  const auto upperTriangleBounds = mapIndexSets<UpperTriangleIndexSets, true>(lower, upper);
  const auto upperTetrangleBounds = mapIndexSets<UpperTetrangleIndexSets, true>(lower, upper);
  const auto lowerTriangleBounds = mapIndexSets<LowerTriangleIndexSets, false>(lower, upper);
  const auto lowerTetrangleBounds = mapIndexSets<LowerTetrangleIndexSets, false>(lower, upper);

  result.klUpperBound = std::min(
    temple::min(upperTriangleBounds),
    temple::min(upperTetrangleBounds)
  );

  result.klLowerBound = std::max({
    0.0,
    temple::max(lowerTriangleBounds),
    temple::max(lowerTetrangleBounds)
  });

  const auto args = std::tie(lower, upper, result.klUpperBound);

  result.u_col = (
    limitAnyOf<UpperTriangleIndexSets, true>(args, upperTriangleBounds)
    || limitAnyOf<UpperTetrangleIndexSets, true>(args, upperTetrangleBounds)
  );

  result.l_col = (
    temple::invoke(zeroBoundTest, args)
    || limitAnyOf<LowerTriangleIndexSets, false>(args, lowerTriangleBounds)
    || limitAnyOf<LowerTetrangleIndexSets, false>(args, lowerTetrangleBounds)
  );

  return result;
}

// Please inline and optimize repeated expressions of me my dear compiler
[[gnu::const]] constexpr inline double sq(const double a) noexcept {
  return a * a;
}

// Calculates D(1, 2, 3; 1, 2, 4)|D_34 = 0
double CMMixed(
  const double d12,
  const double d13,
  const double d14,
  const double d23,
  const double d24
) {
  Eigen::Matrix4d cayleyMengerMatrix;

  /* Zero entries in this matrix:
   * d_ii = 0 (distance to itself is zero)
   * d_34 = 0 (additional condition)
   */
  // Indices colwise             1        2        3
  cayleyMengerMatrix << 0,       1,       1,       1,
                        1,       0, sq(d12), sq(d13),  // 1
                        1, sq(d12),       0, sq(d23),  // 2
                        1, sq(d14), sq(d24),       0;  // 4

  return -cayleyMengerMatrix.determinant() / 4;
}

double CMWith34Zeroed(
  const double d12,
  const double d13,
  const double d14,
  const double d23,
  const double d24
) {
  Eigen::Matrix<double, 5, 5> cayleyMengerMatrix;
  cayleyMengerMatrix << 0,       1,       1,       1,       1,
                        1,       0, sq(d12), sq(d13), sq(d14),
                        1, sq(d12),       0, sq(d23), sq(d24),
                        1, sq(d13), sq(d23),       0,       0,
                        1, sq(d14), sq(d24),       0,       0;

  return cayleyMengerMatrix.determinant() / 8;
}

double CMUpper(
  const double d12,
  const double d13,
  const double d14,
  const double d23,
  const double d24
) {
  const double ATimesMinusFour = std::pow(d12, 2);
  const double B = CMMixed(d12, d13, d14, d23, d24);
  const double C = CMWith34Zeroed(d12, d13, d14, d23, d24);
  const double discriminant = std::pow(B, 2) + ATimesMinusFour * C;

  if(discriminant >= 0) {
    return (-B - std::sqrt(discriminant)) / ((2.0 / -4) * ATimesMinusFour);
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
  const double ATimesMinusFour = std::pow(d12, 2);
  const double B = CMMixed(d12, d13, d14, d23, d24);
  const double C = CMWith34Zeroed(d12, d13, d14, d23, d24);
  const double discriminant = std::pow(B, 2) + ATimesMinusFour * C;

  if(discriminant >= 0) {
    return (-B + std::sqrt(discriminant)) / ((2.0 / -4) * ATimesMinusFour);
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
   * These functions had to be duplicated and fully qualified since we cannot
   * pass templated functions as an overload set.
   */
  return std::max({
    temple::invoke(CMUpper, m::fetchAppropriateBounds(bounds, b, m::l, m::u, m::u, m::u, m::u)),
    temple::invoke(CMUpper, m::fetchAppropriateBounds(bounds, b, m::u, m::l, m::l, m::u, m::u)),
    temple::invoke(CMUpper, m::fetchAppropriateBounds(bounds, b, m::u, m::u, m::u, m::l, m::l))
  });
}

double lowerTetrangleLimit(
  const Eigen::MatrixXd& bounds,
  const std::array<unsigned, 4>& b
) {
  // See upperTetrangleLimit to explain the matrix below
  return std::min({
    temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::u, m::u, m::l, m::l, m::u)),
    temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::u, m::l, m::u, m::u, m::l)),
    temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::l, m::l, m::u, m::l, m::u)),
    temple::invoke(CMLower, m::fetchAppropriateBounds(bounds, b, m::l, m::u, m::l, m::u, m::l))
  });
}

struct TetrangleLimits {
  DistanceGeometry::ValueBounds klLimits;
  bool boundViolation;

  TetrangleLimits(
    const Eigen::MatrixXd& bounds,
    const std::array<unsigned, 4>& b
  ) {
    auto matrixPair = makeLU(bounds, b);
    TriCheckResult check = triCheckRetry(matrixPair.first, matrixPair.second);

    if(check.u_col) {
      klLimits.upper = check.klUpperBound;
    } else {
      klLimits.upper = std::sqrt(upperTetrangleLimit(bounds, b));
    }

    if(check.l_col) {
      klLimits.lower = check.klLowerBound;
    } else {
      klLimits.lower = std::sqrt(lowerTetrangleLimit(bounds, b));
    }

    boundViolation = klLimits.upper < klLimits.upper;
  }
};

Eigen::MatrixXd tetrangleSmooth(Eigen::MatrixXd bounds) {
  const unsigned N = bounds.cols();

  bool changedSomething;
  do {
    changedSomething = false;

    for(unsigned i = 0; i < N - 1; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        // NOTE (i,j) and (k,l) are not mutually disjoint
        for(unsigned k = 0; k < N - 1; ++k) {
          for(unsigned l = k + 1; l < N; ++l) {
            TetrangleLimits limits {bounds, {i, j, k, l}};

            if(limits.boundViolation) {
              throw std::logic_error("Bound violation found!");
            }

            // k < l, so bounds(k, l) is the upper bound, bounds(l, k) the lower
            double& klLowerBound = bounds(l, k);
            double& klUpperBound = bounds(k, l);

            if(
              limits.klLimits.lower > klLowerBound
              && std::fabs(limits.klLimits.lower - klLowerBound) / klLowerBound > 0.01
            ) {
              klLowerBound = limits.klLimits.lower;
              changedSomething = true;
            }

            if(
              limits.klLimits.upper < klUpperBound
              && std::fabs(klUpperBound - limits.klLimits.upper) / klUpperBound > 0.01
            ) {
              klUpperBound = limits.klLimits.upper;
              changedSomething = true;
            }
          }
        }
      }
    }
  } while(changedSomething);

  return bounds;
}

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine
