#ifndef INCLUDE_CYCLIC_POLYGONS_LIB_H
#define INCLUDE_CYCLIC_POLYGONS_LIB_H

#include <boost/math/special_functions/binomial.hpp>

#include <vector>
#include <array>
#include <cassert>

#include "template_magic/VectorView.h"

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

double A5(const std::vector<double>& lambdas) {
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

double B5(const std::vector<double>& lambdas) {
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

double Delta5(const std::vector<double>& lambdas) {
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

} // namespace svrtan

double maximumPentagonCircumradius(const std::vector<double>& edgeLengths) {
  assert(edgeLengths.size() == 5);

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
}

} // namespace CyclicPolygons

#endif
