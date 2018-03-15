#ifndef INCLUDE_CYCLIC_POLYGONS_LIB_H
#define INCLUDE_CYCLIC_POLYGONS_LIB_H

#include "boost/math/tools/roots.hpp"

#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"

#include <cassert>
#include <vector>

/*! @file
 *
 * Provides functionality to calculate the circumradii and internal angles of
 * cyclic polygons of any size.  Contains shortcut calculations of the internal angles
 * for triangles and quadrilaterals, which are easier to treat.
 */

namespace CyclicPolygons {

namespace detail {

template<typename FloatType>
inline FloatType inverseLawOfCosines(
  const FloatType& opposingSideLength,
  const FloatType& adjacentSideLengthA,
  const FloatType& adjacentSideLengthB
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

template<typename FloatType>
std::vector<FloatType> triangleShortcut(const std::vector<FloatType>& edgeLengths) {
  return {
    inverseLawOfCosines(edgeLengths[2], edgeLengths[0], edgeLengths[1]),
    inverseLawOfCosines(edgeLengths[0], edgeLengths[1], edgeLengths[2]),
    inverseLawOfCosines(edgeLengths[1], edgeLengths[0], edgeLengths[2])
  };
}

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

template<typename FloatType>
std::vector<FloatType> quadrilateralShortcut(const std::vector<FloatType>& edgeLengths) {
  return {
    quadrilateralInternalAngle(edgeLengths, {0, 1}, {2, 3}),
    quadrilateralInternalAngle(edgeLengths, {1, 2}, {3, 0}),
    quadrilateralInternalAngle(edgeLengths, {2, 3}, {0, 1}),
    quadrilateralInternalAngle(edgeLengths, {3, 0}, {1, 2}),
  };
}

template<typename FloatType>
FloatType centralAnglesHelper(const FloatType& edgeLength, const FloatType& circumradius) {
  return std::acos(
    1 - (edgeLength * edgeLength) / (
      2 * circumradius * circumradius
    )
  );
}

namespace circumcenterInside {

template<typename FloatType>
std::vector<FloatType> centralAngles(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  return temple::map(
    edgeLengths,
    [&](const FloatType& edgeLength) -> FloatType {
      return centralAnglesHelper(edgeLength, circumradius);
    }
  );
}

template<typename FloatType>
FloatType centralAnglesDeviation(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType& longestEdge
) {
  assert(circumradius > longestEdge / 2);

  return temple::sum(
    centralAngles(circumradius, edgeLengths)
  ) - 2 * M_PI;
}

template<typename FloatType>
FloatType centralAnglesDeviationDerivative(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  return temple::sum(
    temple::map(
      edgeLengths,
      [&](const FloatType& a) -> FloatType {
        return -2 * a / (
          circumradius * std::sqrt(
            4 * circumradius * circumradius - a * a
          )
        );
      }
    )
  );
}

template<typename FloatType>
FloatType centralAnglesDeviationSecondDerivative(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths
) {
  const FloatType squareCircumradius = circumradius * circumradius;

  return temple::sum(
    temple::map(
      edgeLengths,
      [&](const FloatType& a) -> FloatType {
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

template<typename FloatType>
std::vector<FloatType> centralAngles(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType& longestEdge
) {
  return temple::map(
    edgeLengths,
    [&](const FloatType& edgeLength) -> FloatType {
      if(edgeLength == longestEdge) {
        return 2 * M_PI - centralAnglesHelper(edgeLength, circumradius);
      }

      return centralAnglesHelper(edgeLength, circumradius);
    }
  );
}

template<typename FloatType>
FloatType centralAnglesDeviation(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType& longestEdge
) {
  assert(circumradius > longestEdge / 2);

  return temple::sum(
    centralAngles(circumradius, edgeLengths, longestEdge)
  ) - 2 * M_PI;
}

template<typename FloatType>
FloatType centralAnglesDeviationDerivative(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType& longestEdge
) {
  return temple::sum(
    temple::map(
      edgeLengths,
      [&](const FloatType& a) -> FloatType {
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

template<typename FloatType>
FloatType centralAnglesDeviationSecondDerivative(
  const FloatType& circumradius,
  const std::vector<FloatType>& edgeLengths,
  const FloatType& longestEdge
) {
  const FloatType squareCircumradius = circumradius * circumradius;

  return temple::sum(
    temple::map(
      edgeLengths,
      [&](const FloatType& a) -> FloatType {
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

// TODO try upping the precision requirements using numeric_limits
template<typename FloatType>
FloatType regularCircumradius(
  const unsigned& nSides,
  const FloatType& a
) {
  return 0.5 * a / std::sin(M_PI / nSides);
}

//! Precise only to 6 * sizeof(FloatType) bits
template<typename FloatType>
std::pair<FloatType, bool> convexCircumradius(const std::vector<FloatType>& edgeLengths) {
  const unsigned bitPrecision = 6 * sizeof(FloatType);

  const FloatType longestEdge = temple::max(edgeLengths);

  const FloatType lowerBound = longestEdge / 2 + 1e-10;

  const FloatType upperBound = std::numeric_limits<FloatType>::max();

  FloatType rootGuess = regularCircumradius(
    edgeLengths.size(),
    std::max(
      temple::average(edgeLengths),
      lowerBound
    )
  );

  if(rootGuess < lowerBound) {
    rootGuess = lowerBound;
  } 

  assert(lowerBound <= rootGuess && rootGuess <= upperBound);

  auto rootSearchLambda = [&](const FloatType& circumradius) -> std::tuple<FloatType, FloatType, FloatType> {
    return std::make_tuple<FloatType, FloatType, FloatType>(
      circumcenterInside::centralAnglesDeviation(circumradius, edgeLengths, longestEdge),
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

  if(std::fabs(circumcenterInside::centralAnglesDeviation(root, edgeLengths, longestEdge)) >= 1e-6) {
    // Perhaps the circumcenter is outside the polygon?

    auto alternateSearchLambda = [&](const FloatType& circumradius) -> std::tuple<FloatType, FloatType, FloatType> {
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

template<typename FloatType>
std::vector<FloatType> generalizedInternalAngles(
  const std::vector<FloatType>& edgeLengths,
  const FloatType& circumradius,
  bool circumcenterInside
) {
  // Immediately multiply with 2 to avoid doing so in every calculation
  const FloatType doubleR = 2 * circumradius;

  /* Add the first edge length onto the list so we can use pairwiseMap to get
   * all subsequent pairs
   */
  auto lengthsCopy = edgeLengths;
  lengthsCopy.emplace_back(lengthsCopy.front());

  double longestEdge = temple::max(edgeLengths);

  if(circumcenterInside) {
    return temple::mapSequentialPairs(
      lengthsCopy,
      [&](const FloatType& a, const FloatType& b) -> FloatType {
        return std::acos(a / doubleR) + std::acos(b / doubleR);
      }
    );
  }

  auto delta = [&](const FloatType& a) -> FloatType {
    FloatType value = std::acos(a / doubleR);

    if(a == longestEdge) {
      return -value;
    }

    return value;
  };

  return temple::mapSequentialPairs(
    lengthsCopy,
    [&](const FloatType& a, const FloatType& b) -> FloatType {
      return delta(a) + delta(b);
    }
  );
}

} // namespace detail


/*!
 * Returns whether a cyclic polygon exists for the specified sequence of edge
 * lengths.
 */
template<typename FloatType>
bool exists(const std::vector<FloatType>& edgeLengths) {
  /* If a1, a2, ..., aN satisfy: Each edge length smaller than sum of others
   * -> There exists a convex cyclic polygon (Iosif Pinelis, 2005)
   *
   * Equivalent to saying largest value in set of edge lengths smaller than 
   * sum of remainder, no need to check all of them.
   */

  const FloatType maxValue = temple::max(edgeLengths);

  return maxValue < temple::sum(edgeLengths) - maxValue - 0.01 * maxValue;
}

/*!
 * Returns internal angles of the convex cyclic polygon specified by the passed
 * edge lengths. Angles are returned in the following sequence:
 *
 *   edge lengths: a1, a2, ..., aN
 *   angles: a1 ∡ a2, a2 ∡ a3, ..., a(N-1) ∡ aN, aN ∡ a1
 *
 * Requires that the passed vector of edge lengths contains at minimum 3 edges.
 * The cyclic polygon must exist (use exists to check beforehand)
 * and no edge length may be zero (assumes logical error in calling code).
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
    temple::all_of(
      edgeLengths,
      [&](const FloatType& length) -> bool {
        return length != 0;
      }
    ) && "At least one edge length in the sequence is zero "
      "Perhaps consider removing it from the set and approximating it as the "
      "next smaller polygon."
  );

  if(edgeLengths.size() == 3) {
    return detail::triangleShortcut<FloatType>(edgeLengths);
  }

  if(edgeLengths.size() == 4) {
    return detail::quadrilateralShortcut<FloatType>(edgeLengths);
  }

  auto circumradiusResult = detail::convexCircumradius(edgeLengths);

  // General solving scheme
  return detail::generalizedInternalAngles(
    edgeLengths,
    circumradiusResult.first,
    circumradiusResult.second
  );
}

} // namespace CyclicPolygons

#endif
