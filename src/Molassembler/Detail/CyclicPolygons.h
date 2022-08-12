/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Calculate internal angles of cyclic polygons from edge lengths
 *
 * Provides functionality to calculate the circumradii and internal angles of
 * cyclic polygons of any size.  Contains shortcut calculations of the internal angles
 * for triangles and quadrilaterals, which are easier to treat.
 *
 * The general implementation is split into two variants for when the
 * circumcenter is inside the convex cyclic polygon, and when it is outside.
 * The overall strategy is to try out first the variant if it is inside, and
 * if finding a result fails, try if it is outside. One of both must work, but
 * there is no obvious algorithm to determine the position of the circumcenter
 * a priori.
 */

#ifndef INCLUDE_CYCLIC_POLYGONS_LIB_H
#define INCLUDE_CYCLIC_POLYGONS_LIB_H

#include "boost/math/tools/roots.hpp"

#include "Molassembler/Temple/Adaptors/Transform.h"
#include "Molassembler/Temple/Adaptors/SequentialPairs.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include <array>

namespace Scine {
namespace Molassembler {
namespace CyclicPolygons {
namespace Detail {

/*!
 * @brief Helper function to calculate an angle from three side lengths in a
 *   triangle using the law of cosines
 */
template<typename FloatType>
inline FloatType inverseLawOfCosines(
  const FloatType opposingSideLength,
  const FloatType adjacentSideLengthA,
  const FloatType adjacentSideLengthB
) {
  return std::acos(
    (
      adjacentSideLengthA * adjacentSideLengthA
      + adjacentSideLengthB * adjacentSideLengthB
      - opposingSideLength * opposingSideLength
    ) / (
      2 * adjacentSideLengthA * adjacentSideLengthB
    )
  );
}

/*!
 * For triangles, there are no degrees of freedom in arranging the edges at all,
 * and the triangle's inner angles can be calculated directly.
 */
template<typename FloatType>
std::vector<FloatType> triangleShortcut(const std::vector<FloatType>& edgeLengths) {
  return {
    inverseLawOfCosines(edgeLengths[2], edgeLengths[0], edgeLengths[1]),
    inverseLawOfCosines(edgeLengths[0], edgeLengths[1], edgeLengths[2]),
    inverseLawOfCosines(edgeLengths[1], edgeLengths[0], edgeLengths[2])
  };
}

//! Readability of implementation helper
template<typename T>
inline T square(const T& value) {
  return value * value;
}

/* For a cyclic quadrilateral, the internal angle between adjacent edges a and b
 * is given as
 *
 *               a² + b² - c² - d²
 *    cos(phi) = -----------------
 *                 2 (ab + cd)
 *
 * This general structure from adjacent and non-adjacent edge lengths is
 * calculated below.
 */
template<typename FloatType>
FloatType quadrilateralInternalAngle(
  const std::vector<FloatType>& edgeLengths,
  const std::array<unsigned, 2>& adjacentIndices,
  const std::array<unsigned, 2>& nonAdjacentIndices
) {

  return std::acos(
    (
      square(edgeLengths.at(adjacentIndices[0]))
      + square(edgeLengths.at(adjacentIndices[1]))
      - square(edgeLengths.at(nonAdjacentIndices[0]))
      - square(edgeLengths.at(nonAdjacentIndices[1]))
    ) / (
      2 * (
        edgeLengths.at(adjacentIndices[0]) * edgeLengths.at(adjacentIndices[1])
        + edgeLengths.at(nonAdjacentIndices[0]) * edgeLengths.at(nonAdjacentIndices[1])
      )
    )
  );
}

//! The inner angles of a cyclic quadrilateral is also explicitly calculable
template<typename FloatType>
std::vector<FloatType> quadrilateralShortcut(const std::vector<FloatType>& edgeLengths) {
  return {
    quadrilateralInternalAngle(edgeLengths, {{0, 1}}, {{2, 3}}),
    quadrilateralInternalAngle(edgeLengths, {{1, 2}}, {{3, 0}}),
    quadrilateralInternalAngle(edgeLengths, {{2, 3}}, {{0, 1}}),
    quadrilateralInternalAngle(edgeLengths, {{3, 0}}, {{1, 2}}),
  };
}

/*!
 * @brief Calculates the angle for the isosceles triangle spanned from the
 *   circumcenter to an edge
 */
template<typename FloatType>
FloatType centralAnglesHelper(const FloatType edgeLength, const FloatType circumradius) {
  return std::acos(
    1 - (edgeLength * edgeLength) / (
      2 * circumradius * circumradius
    )
  );
}

namespace circumcenterInside {

//! Calculates the angles at the circumcenter for each edge
template<typename FloatType>
std::vector<FloatType> centralAngles(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  return Temple::map(
    edgeLengths,
    [&](const FloatType edgeLength) -> FloatType {
      return centralAnglesHelper(edgeLength, circumradius);
    }
  );
}

/*! Calculates how close to closed a hypothetical convex cyclic polygon is
 *
 * Calculates the deviation of the sum of central angles of the hypothetic
 * cyclic polygon to 2π, which would signify that it is a closed polygon.
 *
 * @param circumradius The circumradius of the hypothetical cyclic polygon
 * @param edgeLengths The adjacent edge lengths composing the polygon
 * @param longestEdge the length of the maximum edge in the edge lengths
 */
template<typename FloatType>
FloatType centralAnglesDeviation(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  return Temple::sum(
    centralAngles(circumradius, edgeLengths)
  ) - 2 * M_PI;
}

//! The first derivative of the central angles deviation
template<typename FloatType>
FloatType centralAnglesDeviationDerivative(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  return Temple::sum(
    Temple::Adaptors::transform(
      edgeLengths,
      [&](const FloatType a) -> FloatType {
        return -2 * a / (
          circumradius * std::sqrt(
            4 * circumradius * circumradius - a * a
          )
        );
      }
    )
  );
}

//! The second derivative of the central angles deviation
template<typename FloatType>
FloatType centralAnglesDeviationSecondDerivative(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  const FloatType squareCircumradius = circumradius * circumradius;

  return Temple::sum(
    Temple::Adaptors::transform(
      edgeLengths,
      [&](const FloatType a) -> FloatType {
        const auto temp = 4 * squareCircumradius - a * a;
        return -2 * a * (
          - 4 * std::pow(temp, -1.5)
          - std::pow(temp, -0.5) / (
            squareCircumradius
          )
        );
      }
    )
  );
}

} // namespace circumcenterInside


namespace circumcenterOutside {

//! Calculates the angles at the circumcenter for each edge
template<typename FloatType>
std::vector<FloatType> centralAngles(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType longestEdge
) {
  return Temple::map(
    edgeLengths,
    [&](const FloatType edgeLength) -> FloatType {
      if(edgeLength == longestEdge) {
        return 2 * M_PI - centralAnglesHelper(edgeLength, circumradius);
      }

      return centralAnglesHelper(edgeLength, circumradius);
    }
  );
}

/*! Calculates how close to closed a hypothetical convex cyclic polygon is
 *
 * Calculates the deviation of the sum of central angles of the hypothetic
 * cyclic polygon to 2π, which would signify that it is a closed polygon.
 *
 * @param circumradius The circumradius of the hypothetical cyclic polygon
 * @param edgeLengths The adjacent edge lengths composing the polygon
 * @param longestEdge the length of the maximum edge in the edge lengths
 */
template<typename FloatType>
FloatType centralAnglesDeviation(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType longestEdge
) {
  assert(circumradius > longestEdge / 2);

  return Temple::sum(
    centralAngles(circumradius, edgeLengths, longestEdge)
  ) - 2 * M_PI;
}

//! The first derivative of the central angles deviation
template<typename FloatType>
FloatType centralAnglesDeviationDerivative(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType longestEdge
) {
  return Temple::sum(
    Temple::Adaptors::transform(
      edgeLengths,
      [&](const FloatType a) -> FloatType {
        FloatType value = -2 * a / (
          circumradius * std::sqrt(
            4 * circumradius * circumradius - a * a
          )
        );

        if(a == longestEdge) {
          return -value;
        }

        return value;
      }
    )
  );
}

//! The second derivative of the central angles deviation
template<typename FloatType>
FloatType centralAnglesDeviationSecondDerivative(
  const FloatType circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType longestEdge
) {
  const FloatType squareCircumradius = circumradius * circumradius;

  return Temple::sum(
    Temple::Adaptors::transform(
      edgeLengths,
      [&](const FloatType a) -> FloatType {
        const auto temp = 4 * squareCircumradius - a * a;
        FloatType value = -2 * a * (
          - 4 * std::pow(temp, -1.5)
          - std::pow(temp, -0.5) / (
            squareCircumradius
          )
        );

        if(a == longestEdge) {
          return -value;
        }

        return value;
      }
    )
  );
}

} // namespace circumcenterOutside

//! The circumradius of a regular polygon composed of a single edge length
template<typename FloatType>
FloatType regularCircumradius(
  const unsigned nSides,
  const FloatType a
) {
  return 0.5 * a / std::sin(M_PI / nSides);
}

/*!
 * @brief Calculates the circumradius of a convex cyclic polygon
 *
 * Calculates the circumradius and whether the circumcenter is inside of the
 * cyclic polygon constructed using the passed adjacent edge lengths.
 *
 * @param edgeLengths adjacent edge lengths composing the polygon
 *
 * @note Precise only to 6 * sizeof(FloatType) bits
 */
template<typename FloatType>
std::pair<FloatType, bool> convexCircumradius(const std::vector<FloatType>& edgeLengths) {
  const unsigned bitPrecision = 6 * sizeof(FloatType);

  const FloatType longestEdge = Temple::max(edgeLengths);

  const FloatType lowerBound = longestEdge / 2 + 1e-10;

  const FloatType upperBound = std::numeric_limits<FloatType>::max();

  FloatType rootGuess = regularCircumradius(
    edgeLengths.size(),
    std::max(
      Temple::average(edgeLengths),
      lowerBound
    )
  );

  if(rootGuess < lowerBound) {
    rootGuess = lowerBound;
  }

  assert(lowerBound <= rootGuess && rootGuess <= upperBound);

  auto rootSearchLambda = [&](const FloatType circumradius) -> std::tuple<FloatType, FloatType, FloatType> {
    return std::make_tuple<FloatType, FloatType, FloatType>(
      circumcenterInside::centralAnglesDeviation(circumradius, edgeLengths),
      circumcenterInside::centralAnglesDeviationDerivative(circumradius, edgeLengths),
      circumcenterInside::centralAnglesDeviationSecondDerivative(circumradius, edgeLengths)
    );
  };

  const unsigned maxIter = 1000;
  boost::uintmax_t iterCount = maxIter;

  auto root = boost::math::tools::schroder_iterate(
    rootSearchLambda,
    rootGuess,
    lowerBound,
    upperBound,
    bitPrecision,
    iterCount
  );

  if(iterCount == static_cast<boost::uintmax_t>(maxIter)) {
    throw std::logic_error("Could not find polygon circumradius!");
  }

  if(std::fabs(circumcenterInside::centralAnglesDeviation(root, edgeLengths)) >= 1e-6) {
    // Perhaps the circumcenter is outside the polygon?

    auto alternateSearchLambda = [&](const FloatType circumradius) -> std::tuple<FloatType, FloatType, FloatType> {
      return std::make_tuple<FloatType, FloatType, FloatType>(
        circumcenterOutside::centralAnglesDeviation(circumradius, edgeLengths, longestEdge),
        circumcenterOutside::centralAnglesDeviationDerivative(circumradius, edgeLengths, longestEdge),
        circumcenterOutside::centralAnglesDeviationSecondDerivative(circumradius, edgeLengths, longestEdge)
      );
    };

    iterCount = maxIter;

    root = boost::math::tools::schroder_iterate(
      alternateSearchLambda,
      rootGuess,
      lowerBound,
      upperBound,
      bitPrecision,
      iterCount
    );

    return {root, false};
  }

  return {root, true};
}

/*!
 * @brief Calculates the internal angles of a polygon when the circumradius and
 *   position of the circumcenter are known
 *
 * @param edgeLengths Consecutive adjacent edge lengths composing the polygon
 * @param circumradius The discovered circumradius on which all vertices can
 *   be placed
 * @param circumcenterInside Whether the circumcenter is inside the final
 *   polygon or outside
 */
template<typename FloatType>
std::vector<FloatType> generalizedInternalAngles(
  const std::vector<FloatType>& edgeLengths,
  const FloatType circumradius,
  bool circumcenterInside
) {
  // Immediately multiply with 2 to avoid doing so in every calculation
  const FloatType doubleR = 2 * circumradius;

  /* Add the first edge length onto the list so we can use pairwiseMap to get
   * all subsequent pairs
   */
  auto lengthsCopy = edgeLengths;
  lengthsCopy.emplace_back(lengthsCopy.front());

  double longestEdge = Temple::max(edgeLengths);

  if(circumcenterInside) {
    return Temple::map(
      Temple::Adaptors::sequentialPairs(lengthsCopy),
      [&](const FloatType a, const FloatType b) -> FloatType {
        return std::acos(a / doubleR) + std::acos(b / doubleR);
      }
    );
  }

  auto delta = [&](const FloatType a) -> FloatType {
    FloatType value = std::acos(a / doubleR);

    if(a == longestEdge) {
      return -value;
    }

    return value;
  };

  return Temple::map(
    Temple::Adaptors::sequentialPairs(lengthsCopy),
    [&](const FloatType a, const FloatType b) -> FloatType {
      return delta(a) + delta(b);
    }
  );
}

} // namespace Detail


/*!
 * @brief Check whether a cyclic polygon exists for a specified sequence of edge
 *   lengths.
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param edgeLengths The lengths of the edges that will compose the polygon.
 *   Positionally adjacent edge lengths will be adjacent in the polygon.
 *
 * @pre Passed edge lengths must be positive.
 *
 * @return Whether the cyclic polygon exists. It does if the longest edge is
 * smaller than the sum of all other edge lengths.
 */
template<typename FloatType>
bool exists(const std::vector<FloatType>& edgeLengths) {
  /* If a1, a2, ..., aN satisfy: Each edge length smaller than sum of others
   * -> There exists a convex cyclic polygon (Iosif Pinelis, 2005)
   *
   * Equivalent to saying largest value in set of edge lengths smaller than
   * sum of remainder, no need to check all of them.
   */

  const FloatType maxValue = Temple::max(edgeLengths);

  return maxValue < Temple::sum(edgeLengths) - maxValue - 0.01 * maxValue;
}

/*!
 * @brief Calculates the internal angles of a convex cyclic polygon
 *
 * @complexity{@math{\Theta(1)} if edgeLengths is three or four, otherwise
 * @math{\Theta(N)}}
 *
 * @param edgeLengths The lengths of the edges that will compose the polygon.
 *   Positionally adjacent edge lengths will be adjacent in the polygon.
 *
 * @pre You must pass at least three positive edge lengths. The cyclic polygon
 *   must exist (use exists to check beforehand) and no edge length may be zero
 *   (assumes logical error in calling code).
 *
 * @return Internal angles of the convex cyclic polygon specified by the passed
 *   edge lengths. Angles are returned in the following sequence:
 *   - passed edge lengths: a1, a2, ..., aN
 *   - angles: a1 ∡ a2, a2 ∡ a3, ..., a(N-1) ∡ aN, aN ∡ a1
 */
template<typename FloatType>
std::enable_if_t<
  std::is_floating_point<FloatType>::value,
  std::vector<FloatType>
> internalAngles(const std::vector<FloatType>& edgeLengths) {
  assert(
    exists(edgeLengths)
    && "The passed sequence of lengths cannot be used to construct a polygon. "
      "An edge length surpassed the sum of the lengths of all others. This is "
      "the necessary condition for the existence of a cyclic polygon."
  );
  assert(
    edgeLengths.size() >= 3
    && "It is unreasonable to call this for less than three edges."
  );
  assert(
    Temple::all_of(
      edgeLengths,
      [&](const FloatType length) -> bool {
        return length != 0;
      }
    ) && "At least one edge length in the sequence is zero "
      "Perhaps consider removing it from the set and approximating it as the "
      "next smaller polygon."
  );

  if(edgeLengths.size() == 3) {
    return Detail::triangleShortcut<FloatType>(edgeLengths);
  }

  if(edgeLengths.size() == 4) {
    return Detail::quadrilateralShortcut<FloatType>(edgeLengths);
  }

  auto circumradiusResult = Detail::convexCircumradius(edgeLengths);

  // General solving scheme
  return Detail::generalizedInternalAngles(
    edgeLengths,
    circumradiusResult.first,
    circumradiusResult.second
  );
}

} // namespace CyclicPolygons
} // namespace Molassembler
} // namespace Scine

#endif
