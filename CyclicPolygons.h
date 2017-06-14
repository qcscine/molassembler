#ifndef INCLUDE_CYCLIC_POLYGONS_LIB_H
#define INCLUDE_CYCLIC_POLYGONS_LIB_H

#include "constexpr_magic/Math.h"
#include "template_magic/TemplateMagic.h"

#include <boost/optional.hpp>

#include <vector>

namespace CyclicPolygons {

namespace detail {

template<class CallFunction>
double recursiveLoopSummationImpl(
  CallFunction&& callFunction,
  std::vector<unsigned> indices,
  const unsigned& nLoops,
  const unsigned& upperLimit
);

template<class CallFunction>
double recursiveLoopSummation(
  CallFunction&& callFunction,
  const unsigned& nLoops,
  const unsigned upperLimit
);

} // namespace detail

namespace analysis {

void writeAnalysisFiles(
  const std::vector<double> edgeLengths,
  const std::string& baseName
);

} // namespace analysis

namespace math {

static constexpr double radiusRhoConversionConstant = (
  10 - 2 * ConstexprMagic::Math::sqrt(5.0)
) / 4.0;

double elementarySymmetricPolynomial(
  const unsigned& k,
  const std::vector<double>& values
);

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
);

template<typename T>
inline T square(const T& value) {
  return value * value;
}

} // namespace math

namespace svrtan {

double lambda(
  const unsigned& k,
  const double& rho, 
  const std::vector<double>& epsilon
);

inline double A5(const std::vector<double>& lambdas) {
  return TemplateMagic::numeric::kahanSum(
    std::vector<double>{
      std::pow(lambdas[4], 4),
      (
        - 3 * lambdas[3]
        + 2 * lambdas[2]
        + lambdas[1]
        - 3
      ) * math::square(lambdas[4]),
      (
        - 2 * lambdas[3]
        - 4 * lambdas[1]
        + 2
      ) * lambdas[4],
      2 * math::square(lambdas[3]),
      (
        - 2 * lambdas[2]
        - 2 * lambdas[1]
        + 4 
      ) * lambdas[3],
      math::square(lambdas[2]),
      2 * lambdas[2],
      - 2 * lambdas[1],
      (
        lambdas[3] + 3
      ) * lambdas[0],
      2
    }
  );
}

inline double B5(const std::vector<double>& lambdas) {
   return TemplateMagic::numeric::kahanSum(
    std::vector<double>{
      - std::pow(lambdas[4], 3),
      2 * math::square(lambdas[4]),
      (
         2 * lambdas[3] - lambdas[2]
      ) * lambdas[4],
      - 2 * lambdas[3],
      2 * lambdas[1],
      - lambdas[0],
      - 2
    }
  );
}

inline double Delta5(const std::vector<double>& lambdas) {
  return (
    lambdas[0]
    + 2 * TemplateMagic::numeric::kahanSum(
      std::vector<double>{
        lambdas[1],
        lambdas[2],
        lambdas[3],
        lambdas[4],
        1
      }
    )
  );
}

double pentagon(
  const double& rho,
  const std::vector<double>& epsilon
);

} // namespace svrtan

bool cyclicPolygonConstructible(const std::vector<double>& edgeLengths);

boost::optional<double> maximumPentagonCircumradius(const std::vector<double>& edgeLengths);

double maximumQuadrilateralCircumradius(const std::vector<double> edgeLengths);

double inverseLawOfCosines(
  const double& opposingSideLength,
  const double& adjacentSideLengthA,
  const double& adjacentSideLengthB
);

std::vector<double> internalAnglesTriangle(const std::vector<double>& edgeLengths);

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
);

std::vector<double> internalAnglesQuadrilateral(const std::vector<double>& edgeLengths);

/* In a cyclic pentagon, if the circumradius is known, then isosceles triangles
 * can be spanned from every edge to the center of the circle and the base
 * angles calculated via inverting the law of cosines. The angle between two
 * edges of the pentagon is then merely the sum of the neighboring internal
 * triangle angles.
 */
std::vector<double> internalAnglesPentagon(
  const std::vector<double>& edgeLengths,
  const double& circumradius
);

/*! 
 * Returns whether a guess for a root is valid or not by calculating all
 * internal angles and checking whether the sum is very close to 3π.
 */
bool validateRhoGuess(
  const std::vector<double>& edgeLengths,
  const double& rhoGuess
);

/*! Returns internal angles in the following sequence
 * edge lengths: a1, a2, ..., aN
 * angles: a1 ∡ a2, a2 ∡ a3, ..., a(N-1) ∡ aN, aN ∡ a1
 */
std::vector<double> internalAngles(const std::vector<double> edgeLengths);

} // namespace CyclicPolygons

#endif
