#include "CyclicPolygons.h"

#include "template_magic/VectorView.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/tools/roots.hpp>

#include <fstream>
#include <cassert>

namespace CyclicPolygons {

namespace detail {

template<class CallFunction>
double recursiveLoopSummationImpl(
  CallFunction&& callFunction,
  std::vector<unsigned> indices,
  const unsigned& nLoops,
  const unsigned& upperLimit
) {
  if(indices.size() == nLoops) {
    return callFunction(indices);
  }
  
  double sumValues = 0;

  for(
    unsigned i = indices.back() + 1;
    i <= upperLimit - nLoops + indices.size();
    ++i
  ) {
    auto indicesCopy = indices;
    indicesCopy.emplace_back(i);

    sumValues += recursiveLoopSummationImpl(
      std::forward<CallFunction>(callFunction),
      indicesCopy,
      nLoops,
      upperLimit
    );
  }

  return sumValues;
}

template<class CallFunction>
double recursiveLoopSummation(
  CallFunction&& callFunction,
  const unsigned& nLoops,
  const unsigned upperLimit
) {
  assert(
    nLoops > 0 
    && "Number of loops ought to be > 0, else why are you calling this function?"
  );

  double sumValues = 0;
  for(unsigned i = 0; i <= upperLimit - nLoops; ++i) {
    sumValues += recursiveLoopSummationImpl(
      std::forward<CallFunction>(callFunction),
      {i},
      nLoops,
      upperLimit
    );
  }

  return sumValues;
}

} // namespace detail

namespace analysis {

void writeAnalysisFiles(
  const std::vector<double>& edgeLengths,
  const std::string& baseName
) {
  using namespace std::string_literals;

  const auto squaredEdgeLengths = TemplateMagic::map(
    edgeLengths,
    CyclicPolygons::math::square<double>
  );

  const auto epsilon = TemplateMagic::map(
    CyclicPolygons::math::intSeq(0, 5),
    [&](const unsigned& k) -> double {
      return CyclicPolygons::math::elementarySymmetricPolynomial(
        k,
        squaredEdgeLengths
      );
    }
  );

  const auto unaryPolynomial = [&](const double& rho) -> double {
    return CyclicPolygons::Pentagon::svrtanPolynomial(rho, epsilon);
  };

  std::string scanFilename = baseName + ".csv"s;
  std::string metaFilename = baseName + "-meta.csv"s;

  std::ofstream outFile(scanFilename);
  outFile << std::scientific << std::setprecision(8);

  const double trialFactor = 1.2;

  const double average = TemplateMagic::numeric::average(edgeLengths);
  const double stddev = TemplateMagic::numeric::stddev(edgeLengths, average);

  const double rhoGuess = 1 / CyclicPolygons::math::square(
    0.5 * average / std::sin(M_PI / 5)
  );
  const double rhoStddev = 2 * stddev / (
    std::pow(rhoGuess, 3) * std::sin(M_PI / edgeLengths.size())
  );
  const bool triedFirstStrategy = (stddev < average * 0.1);

  boost::optional<double> rhoRoot;

  if(triedFirstStrategy) {
    const double lowerValue = unaryPolynomial(rhoGuess / trialFactor);
    const double upperValue = unaryPolynomial(rhoGuess * trialFactor);
    bool increasing = lowerValue < upperValue;

    const unsigned maxIter = 1000;
    boost::uintmax_t iterCount = maxIter; // maxIter
    auto returnBounds = boost::math::tools::bracket_and_solve_root(
      unaryPolynomial,
      rhoGuess,
      trialFactor,
      increasing,
      boost::math::tools::eps_tolerance<double>(32),
      iterCount
    );

    if(iterCount < static_cast<boost::uintmax_t>(maxIter)) {
      const double rhoGuess = (returnBounds.first + returnBounds.second) / 2;
      if(CyclicPolygons::Pentagon::validateRhoGuess(edgeLengths, rhoGuess)) {
        // Success! 
        rhoRoot = rhoGuess;
      }
    }
  }

  const bool usedFallbackStrategy = (!rhoRoot);

  const double scanLower = (
    CyclicPolygons::math::radiusRhoConversionConstant / CyclicPolygons::math::square(
      TemplateMagic::numeric::max(edgeLengths)
    )
  );
  const double scanUpper = (
    CyclicPolygons::math::radiusRhoConversionConstant / CyclicPolygons::math::square(
      TemplateMagic::numeric::min(edgeLengths)
    )
  );
  const unsigned nSteps = 1000;
  const double stepLength = (scanUpper - scanLower) / nSteps;
  const double initialValue = unaryPolynomial(scanLower);

  bool previousSignBit = std::signbit(initialValue);

  outFile << scanLower << ", " << initialValue << std::endl;

  for(unsigned i = 1; i < nSteps; i++) {
    const double rho = scanLower + i * stepLength;
    const double value = unaryPolynomial(rho);

    outFile << rho << ", " << value << std::endl;

    if(
      std::signbit(unaryPolynomial(rho)) != previousSignBit
      && !rhoRoot
    ) {
      const double lowerBound = rho - stepLength;

      const double lowerValue = unaryPolynomial(lowerBound);
      const double upperValue = unaryPolynomial(rho);
      
      const unsigned maxIter = 1000;
      boost::uintmax_t iterCount = maxIter;

      auto returnBounds = boost::math::tools::toms748_solve(
        unaryPolynomial,
        lowerBound,
        rho,
        lowerValue,
        upperValue,
        boost::math::tools::eps_tolerance<double>(32),
        iterCount
      );

      if(iterCount < static_cast<boost::uintmax_t>(maxIter)) {
        const double rhoGuess = (returnBounds.first + returnBounds.second) / 2;
        if(CyclicPolygons::Pentagon::validateRhoGuess(edgeLengths, rhoGuess)) {
          // Success! 
          rhoRoot = rhoGuess;
        }
      } 

      previousSignBit = std::signbit(upperValue);
    }
  }

  outFile.close();

  boost::optional<bool> validRoot;

  // Analyze returned root
  if(rhoRoot) {
    const double doubleR = 2 * sqrt(1 / rhoRoot.value());

    auto lengthsCopy = edgeLengths;
    lengthsCopy.emplace_back(lengthsCopy.front());

    auto internalAngles = TemplateMagic::pairwiseMap(
      lengthsCopy,
      [&doubleR](const double& a, const double& b) -> double {
        return (
          std::acos(
            a / doubleR
          ) + std::acos(
            b / doubleR
          )
        );
      }
    );

    const double angleSum = TemplateMagic::numeric::sum(internalAngles);
    validRoot = (std::fabs(angleSum - 3 * M_PI) < 1e-6);
  }

  std::ofstream metaFile(metaFilename);
  for(unsigned i = 0; i < edgeLengths.size(); i++) {
    metaFile << edgeLengths[i];
    if(i != edgeLengths.size() - 1) {
      metaFile << ", ";
    }
  }

  metaFile << std::endl;
  metaFile << triedFirstStrategy << ", " << usedFallbackStrategy << std::endl;
  metaFile << rhoGuess << ", " << rhoStddev << std::endl;
  metaFile << rhoRoot.value_or(0) << ", " << validRoot.value_or(0) << std::endl;

  metaFile.close();
}

} // namespace analysis

namespace math {

double elementarySymmetricPolynomial(
  const unsigned& k,
  const std::vector<double>& values
) {
  if(k == 0) {
    return 1;
  }

  return detail::recursiveLoopSummation(
    [&values](const std::vector<unsigned>& indices) -> double {
      return TemplateMagic::reduce(
        TemplateMagic::subset(
          values,
          indices
        ),
        1.0,
        std::multiplies<double>()
      );
    },
    k,
    values.size()
  );
}

/* Generates a vector containing the unsigned integer sequence 
 * [start, start+ 1, ..., end] (end inclusive). If start and end are the same,
 * the sequence [start] is returned.
 *
 * This could be replaced with a struct that offers iterators that merely
 * generate the sequence, but do not store it. The implementation is okay for
 * short sequences, but potentially incredibly wasteful for large ones. 
 */
std::vector<unsigned> intSeq(
  const unsigned& start,
  const unsigned& end
) {
  assert(start <= end);

  std::vector<unsigned> sequence (end - start + 1);

  std::iota(
    sequence.begin(),
    sequence.end(),
    start
  );

  return sequence;
}

double inverseLawOfCosines(
  const double& opposingSideLength,
  const double& adjacentSideLengthA,
  const double& adjacentSideLengthB
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

} // namespace math

namespace Pentagon {

double lambda(
  const unsigned& k,
  const double& rho, 
  const std::vector<double>& epsilon
) {
  assert(epsilon.size() == 6);
  assert(k <= 4);

  return TemplateMagic::reduce(
    TemplateMagic::map(
      math::intSeq(k, 5),
      [&](const auto& i) -> double {
        return (
          boost::math::binomial_coefficient<double>(2 * i, i - k)
          * std::pow(-1, 5 - i)
          * epsilon.at(5 - i)
          * std::pow(rho, 5 - i)
        );
      }
    ),
    0.0,
    std::plus<double>()
  );
}

double svrtanPolynomial(
  const double& rho,
  const std::vector<double>& epsilon
) {
  const auto lambdas = TemplateMagic::map(
    math::intSeq(0, 4),
    [&](const unsigned& k) -> double {
      return lambda(
        k,
        rho,
        epsilon
      );
    }
  );

  return (
    math::square(
      A5(lambdas)
    ) 
    - math::square(
      B5(lambdas)
    ) * Delta5(lambdas)
  );
}

boost::optional<double> convexCircumradiusSvrtan(const std::vector<double>& edgeLengths) {
  assert(edgeLengths.size() == 5);

  assert(
    cyclicPolygonConstructible(edgeLengths) 
    && "No pentagon constructible with given edge lengths! You should've "
      "checked for this edge case before calling convexCircumradiusSvrtan!"
  );

  const auto squaredEdgeLengths = TemplateMagic::map(
    edgeLengths,
    math::square<double>
  );

  const auto epsilon = TemplateMagic::map(
    math::intSeq(0, 5),
    [&](const unsigned& k) -> double {
      return math::elementarySymmetricPolynomial(
        k,
        squaredEdgeLengths
      );
    }
  );

  auto unaryPolynomial = [&](const double& rho) -> double {
    return Pentagon::svrtanPolynomial(rho, epsilon);
  };

  const double average = TemplateMagic::numeric::average(edgeLengths);
  const double stddev = TemplateMagic::numeric::stddev(edgeLengths, average);

  if(stddev < average * 0.1) {
    /* In the special case where the standard deviation is small compared to the
     * average, we can approximate the circumradius with a regular polygon and
     * then search in the local vicinity.
     *
     * We use the formula for regular polygons with side length a and number of
     * vertices n:
     * 
     *   r = 0.5 * a * csc(π / n)
     *
     * We approximate a ≈ ā (average of all a_n) and propagate the standard
     * deviation of the edge length average through to the calculation of rho,
     * which is rho = 1 / r²
     */
    const double rhoGuess = 1 / math::square(
      0.5 * average / std::sin(M_PI / 5)
    );
    const double trialFactor = 1.1;

    const double lowerValue = unaryPolynomial(rhoGuess / trialFactor);
    const double upperValue = unaryPolynomial(rhoGuess * trialFactor);

    // Ensure there is a crossing at all, if not calling bracket_solve would fail
    if(std::signbit(lowerValue) != std::signbit(upperValue)) { 
      bool increasing = lowerValue < upperValue;

      const unsigned maxIter = 1000;
      boost::uintmax_t iterCount = maxIter; // maxIter
      auto returnBounds = boost::math::tools::bracket_and_solve_root(
        unaryPolynomial,
        rhoGuess,
        trialFactor,
        increasing,
        boost::math::tools::eps_tolerance<double>(32),
        iterCount
      );

      if(iterCount < static_cast<boost::uintmax_t>(maxIter)) {
        // Validate the rho guess
        const double rhoRoot = (returnBounds.first + returnBounds.second) / 2;
        if(validateRhoGuess(edgeLengths, rhoRoot)) {
          // Success! We can return the circumradius
          return std::sqrt(1 / rhoRoot);
        }
      }
    }
  }

  // Fallback strategy -> Scan the function until the first crossing
  const double scanLower = (
    math::radiusRhoConversionConstant / CyclicPolygons::math::square(
      TemplateMagic::numeric::max(edgeLengths)
    )
  );
  const double scanUpper = (
    math::radiusRhoConversionConstant / CyclicPolygons::math::square(
      TemplateMagic::numeric::min(edgeLengths)
    )
  );
  const unsigned Nsteps = 1000;
  const double stepLength = (scanUpper - scanLower) / Nsteps;

  bool previousSignBit = std::signbit(unaryPolynomial(scanLower));
  double rho;
  for(unsigned i = 1; i < Nsteps; i++) {
    rho = scanLower + i * stepLength;
    if(std::signbit(unaryPolynomial(rho)) != previousSignBit) {
      const double lowerBound = rho - stepLength;

      const double lowerValue = unaryPolynomial(lowerBound);
      const double upperValue = unaryPolynomial(rho);
      
      const unsigned maxIter = 1000;
      boost::uintmax_t iterCount = maxIter;

      auto returnBounds = boost::math::tools::toms748_solve(
        unaryPolynomial,
        lowerBound,
        rho,
        lowerValue,
        upperValue,
        boost::math::tools::eps_tolerance<double>(32),
        iterCount
      );

      if(iterCount < static_cast<boost::uintmax_t>(maxIter)) {
        // Validate the root
        const double rhoRoot = (returnBounds.first + returnBounds.second) / 2;
        if(validateRhoGuess(edgeLengths, rhoRoot)) {
          // Success! We can return the circumradius
          return std::sqrt(1 / rhoRoot);
        }
      }

      /* Overwrite the previous signbit in case crossing was achieved but the 
       * found root is invalid!
       */
      previousSignBit = std::signbit(upperValue);
    }
  }

  return boost::none;
}

double circumradius(const std::vector<double>& edgeLengths) {
  const double minR = TemplateMagic::numeric::max(edgeLengths) / 2 + 1e-10;
  const double lowerBound = std::max(
    regularCircumradius(
      TemplateMagic::numeric::min(edgeLengths)
    ),
    minR
  );
  const double upperBound = std::max(
    regularCircumradius(
      TemplateMagic::numeric::max(edgeLengths)
    ), 
    minR
  );
  const double rootGuess = regularCircumradius(
    std::max(
      TemplateMagic::numeric::average(edgeLengths),
      minR
    )
  );

  auto rootSearchLambda = [&](const double& circumradius) -> std::tuple<double, double, double> {
    return std::make_tuple<double, double, double>(
      centralAnglesDeviation(circumradius, edgeLengths),
      centralAnglesDeviationDerivative(circumradius, edgeLengths),
      centralAnglesDeviationSecondDerivative(circumradius, edgeLengths)
    );
  };

  const unsigned maxIter = 1000;
  boost::uintmax_t iterCount = maxIter;

  auto root = boost::math::tools::schroder_iterate(
    rootSearchLambda,
    rootGuess,
    lowerBound,
    upperBound,
    32, // bits precision
    iterCount
  );

  if(iterCount == static_cast<boost::uintmax_t>(maxIter)) {
    throw std::logic_error("Could not find pentagon circumradius!");
  }

  return root;
}

/* In a cyclic pentagon, if the circumradius is known, then isosceles triangles
 * can be spanned from every edge to the center of the circle and the base
 * angles calculated via inverting the law of cosines. The angle between two
 * edges of the pentagon is then merely the sum of the neighboring internal
 * triangle angles.
 */
std::vector<double> internalAngles(
  const std::vector<double>& edgeLengths,
  const double& circumradius
) {
  // Immediately multiply with 2 to avoid doing so in every calculation
  const double doubleR = 2 * circumradius;

  /* Add the first edge length onto the list so we can use pairwiseMap to get
   * all subsequent pairs
   */
  auto lengthsCopy = edgeLengths;
  lengthsCopy.emplace_back(lengthsCopy.front());

  return TemplateMagic::pairwiseMap(
    lengthsCopy,
    [&doubleR](const double& a, const double& b) -> double {
      return (
        std::acos(
          a / doubleR
        ) + std::acos(
          b / doubleR
        )
      );
    }
  );
}

std::vector<double> centralAngles(
  const double& circumradius,
  const std::vector<double>& edgeLengths
) {
  return TemplateMagic::map(
    edgeLengths,
    [&](const double& edgeLength) -> double {
      return std::acos(
        1 - (edgeLength * edgeLength) / (
          2 * circumradius * circumradius
        )
      );
    }
  );
}

double centralAnglesDeviation(
  const double& circumradius,
  const std::vector<double>& edgeLengths
) {
  assert(edgeLengths.size() == 5);
  assert(circumradius > TemplateMagic::numeric::max(edgeLengths) / 2);

  return TemplateMagic::numeric::sum(
    centralAngles(circumradius, edgeLengths)
  ) - 2 * M_PI;
}

double centralAnglesDeviationDerivative(
  const double& circumradius,
  const std::vector<double>& edgeLengths
) {
  return TemplateMagic::numeric::sum(
    TemplateMagic::map(
      edgeLengths,
      [&](const double& a) -> double {
        return -2 * a / (
          circumradius * std::sqrt(
            4 * circumradius * circumradius - a * a
          )
        );
      }
    )
  );
}

double centralAnglesDeviationSecondDerivative(
  const double& circumradius,
  const std::vector<double>& edgeLengths
) {
  const double squareCircumradius = circumradius * circumradius;

  return TemplateMagic::numeric::sum(
    TemplateMagic::map(
      edgeLengths,
      [&](const double& a) -> double {
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

bool validateRhoGuess(
  const std::vector<double>& edgeLengths,
  const double& rhoGuess
) {
  const double circumradius = std::sqrt(1 / rhoGuess);
  auto internalAngles = Pentagon::internalAngles(edgeLengths, circumradius);
  return (
    TemplateMagic::numeric::sum(internalAngles) - 3 * M_PI
    < 1e-6
  );
}


} // namespace Pentagon

namespace Quadrilateral {

double convexCircumradius(const std::vector<double>& edgeLengths) {
  assert(edgeLengths.size() == 5);

  assert(
    cyclicPolygonConstructible(edgeLengths) 
    && "No quadrilateral constructible with given edge lengths! You should've "
      "checked for this edge case before calling maximumQuadrilateralCircumradius!"
  );

  double s = TemplateMagic::numeric::sum(edgeLengths) / 2;
  /* (
   *   (a * c + b * d)
   *   * (a * d + b * c)
   *   * (a * b + c * d)
   * ) / (
   *   (s - a)
   *   * (s - b)
   *   * (s - c)
   *   * (s - d)
   * )
   */
  return std::sqrt(
    (
      (edgeLengths[0] * edgeLengths[2] + edgeLengths[1] * edgeLengths[3])
      * (edgeLengths[0] * edgeLengths[3] + edgeLengths[1] * edgeLengths[2])
      * (edgeLengths[0] * edgeLengths[1] + edgeLengths[2] * edgeLengths[3])
    ) / TemplateMagic::reduce(
      TemplateMagic::map(
        edgeLengths,
        [&](const double& edgeLength) -> double {
          return s - edgeLength;
        }
      ),
      1.0,
      std::multiplies<double>()
    )
  );
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
double quadrilateralInternalAngle(
  const std::pair<double, double>& adjacentEdgeLengths,
  const std::pair<double, double>& nonAdjacentEdgeLengths
) {
  return std::acos(
    (
      math::square(adjacentEdgeLengths.first)
      + math::square(adjacentEdgeLengths.second)
      - math::square(nonAdjacentEdgeLengths.first)
      - math::square(nonAdjacentEdgeLengths.second)
    ) / (
      2 * (
        adjacentEdgeLengths.first * adjacentEdgeLengths.second 
        + nonAdjacentEdgeLengths.first * nonAdjacentEdgeLengths.second 
      )
    )
  );
}

std::vector<double> internalAngles(const std::vector<double>& edgeLengths) {
  return {
    quadrilateralInternalAngle(
      {edgeLengths[0], edgeLengths[1]},
      {edgeLengths[2], edgeLengths[3]}
    ),
    quadrilateralInternalAngle(
      {edgeLengths[1], edgeLengths[2]},
      {edgeLengths[3], edgeLengths[0]}
    ),
    quadrilateralInternalAngle(
      {edgeLengths[2], edgeLengths[3]},
      {edgeLengths[0], edgeLengths[1]}
    ),
    quadrilateralInternalAngle(
      {edgeLengths[3], edgeLengths[0]},
      {edgeLengths[1], edgeLengths[2]}
    )
  };
}

} // namespace Quadrilateral

namespace Triangle {

std::vector<double> internalAngles(const std::vector<double>& edgeLengths) {
  return {
    math::inverseLawOfCosines(edgeLengths[2], edgeLengths[0], edgeLengths[1]),
    math::inverseLawOfCosines(edgeLengths[0], edgeLengths[1], edgeLengths[2]),
    math::inverseLawOfCosines(edgeLengths[1], edgeLengths[0], edgeLengths[2])
  };
}

} // namespace Triangle

bool cyclicPolygonConstructible(const std::vector<double>& edgeLengths) {
  /* If a1, a2, ..., aN satisfy: Each edge length smaller than sum of others
   * -> There exists a convex cyclic polygon (Iosif Pinelis, 2005)
   *
   * Equivalent to saying largest value in set of edge lengths smaller than 
   * remainder, no need to check all of them.
   */

  const double maxValue = TemplateMagic::numeric::max(edgeLengths);

  return (
    maxValue <
    TemplateMagic::numeric::sum(edgeLengths) - maxValue 
  );
}

/*! Returns internal angles in the following sequence
 * edge lengths: a1, a2, ..., aN
 * angles: a1 ∡ a2, a2 ∡ a3, ..., a(N-1) ∡ aN, aN ∡ a1
 */
std::vector<double> internalAngles(const std::vector<double>& edgeLengths) {
  assert(
    cyclicPolygonConstructible(edgeLengths)
    && "The passed sequence of lengths cannot be used to construct a polygon. "
      "An edge length surpassed the sum of the lengths of all others. This is "
      "the necessary condition for the existence of a cyclic polygon."
  );
  assert(
    3 <= edgeLengths.size() && edgeLengths.size() <= 5
    && "This function can only handle triangles, quadrilaterals and pentagons."
  );
  assert( 
    TemplateMagic::all_of(
      TemplateMagic::map(
        edgeLengths,
        [&](const double& length) -> bool {
          return length != 0;
        }
      )
    ) && "At least one edge length in the sequence is zero "
      "Perhaps consider removing it from the set and approximating it as the "
      "next smaller polygon."
  );

  if(edgeLengths.size() == 3) {
    return Triangle::internalAngles(edgeLengths);
  } 
  
  if(edgeLengths.size() == 4) {
    return Quadrilateral::internalAngles(edgeLengths);
  } 

  return Pentagon::internalAngles(
    edgeLengths, 
    Pentagon::convexCircumradius(edgeLengths)
  );
}

} // namespace CyclicPolygons
