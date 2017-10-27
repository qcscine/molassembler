#ifndef INCLUDE_DG_DLIB_REFINEMENT_PROBLEM_H
#define INCLUDE_DG_DLIB_REFINEMENT_PROBLEM_H

#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "common_typedefs.h"

#include <dlib/optimization.h>

/*! @file
 *
 * Specifies the refinement minimization problem for the dlib library.
 */

namespace MoleculeManip {

namespace DistanceGeometry {

namespace errfDetail {
  using Vector = dlib::matrix<double, 0, 1>;

  inline dlib::vector<double, 3> getPos3D(const Vector& positions, const unsigned& i) {
    assert(4 * i + 3 < positions.size());

    return dlib::rowm(
      positions,
      dlib::range(4 * i, 4 * i + 2)
    );
  }

  inline Vector getPos(const Vector& positions, const unsigned& i) {
    assert(4 * i + 3 < positions.size());

    return dlib::rowm(
      positions, 
      dlib::range(4 * i, 4 * i + 3)
    );
  }

  inline double volume(const Vector& positions, const std::array<AtomIndexType, 4>& indices) {
    return (
      getPos3D(positions, indices[0]) - getPos3D(positions, indices[3])
    ).dot(
      (
        getPos3D(positions, indices[1]) - getPos3D(positions, indices[3])
      ).cross(
        getPos3D(positions, indices[2]) - getPos3D(positions, indices[3])
      )
    );
  }

  double proportionChiralityConstraintsCorrectSign(
    const std::vector<ChiralityConstraint>& chiralityConstraints,
    const Vector& positions
  );

  bool finalStructureAcceptable(
    const DistanceBoundsMatrix& bounds,
    const std::vector<ChiralityConstraint> chiralityConstraints,
    const Vector& positions
  );
} // namespace errfDetail

template<bool compress>
struct errfValue {
  // Typedefs
  using Vector = dlib::matrix<double, 0, 1>;

  // State
  const unsigned N;
  const dlib::matrix<double, 0, 0> squaredBounds;
  const std::vector<ChiralityConstraint>& chiralityConstraints;

  // Constructor
  explicit errfValue(
    const DistanceBoundsMatrix& distanceBounds,
    const std::vector<ChiralityConstraint>& chiralityConstraints
  ) : N(distanceBounds.N),
      squaredBounds(
        dlib::mat( // Converts Eigen::MatrixXd to dlib::matrix
          distanceBounds.makeSquaredBoundsMatrix()
        )
      ),
      chiralityConstraints(chiralityConstraints)
  {}

  // Information
  double distanceError(const Vector& positions) const {
    double error = 0, lowerTerm, upperTerm;

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = i + 1; j < N; j++) {
        const double& upperBoundSquared = squaredBounds(i, j);
        const double& lowerBoundSquared = squaredBounds(j, i);
        // Since i < j, upper Bound is (i, j), lower Bound is (j, i)
        assert(lowerBoundSquared < upperBoundSquared);

        const double diffLength = dlib::length_squared(
          errfDetail::getPos(positions, j) - errfDetail::getPos(positions, i)
        );

        // Upper bound term
        upperTerm = diffLength / upperBoundSquared - 1;

        assert(!( // Ensure early continue logic is correct
          upperTerm > 0 
          && (
            2 * lowerBoundSquared / (
              lowerBoundSquared + diffLength
            ) - 1
          ) > 0
        ));

        if(upperTerm > 0) {
          error += upperTerm * upperTerm;

          // If the upper term contributes, the lower certainly doesn't
          continue;
        }

        // Lower bound term
        lowerTerm = 2 * lowerBoundSquared / (
          lowerBoundSquared + diffLength
        ) - 1;

        if(lowerTerm > 0) {
          error += lowerTerm * lowerTerm;
        }
      }
    }

    return error;
  }

  double chiralError(const Vector& positions) const {
    double error = 0, volume, upperTerm, lowerTerm;

    for(const auto& constraint : chiralityConstraints) {
      volume = errfDetail::volume(positions, constraint.indices);

      // Upper bound term
      upperTerm = volume - constraint.upper;

      assert(!( // Ensure early continue logic is correct
        upperTerm > 0
        && constraint.lower - volume > 0
      ));

      if(upperTerm > 0) {
        error += upperTerm * upperTerm;
        
        // If the upper term contributes, the lower certainly doesn't
        continue;
      }

      // Lower bound term
      lowerTerm = constraint.lower - volume;
      if(lowerTerm > 0) {
        error += lowerTerm * lowerTerm;
      }
    }

    return error;
  }

  double extraDimensionError(const Vector& positions) const {
    if(!compress) {
      return 0;
    }

    double error = 0;
    for(unsigned i = 0; i < N; i++) {
      const double& w = positions(4 * i + 3);
      error += w * w;
    }

    return error;
  }

  // Operators
  double operator() (const Vector& positions) const {
    assert(positions.size() % 4 == 0);
    assert(positions.size() / 4 == N);

    // Evaluation of addition proceeds left-to-right
    if(compress) {
      // After chiral inversion, generally: distance >> chiral > extraDim
      return (
        extraDimensionError(positions)
        + chiralError(positions)
        + distanceError(positions)
      );
    }

    // Before chiral inversion, typically: chiral >> distance > extraDim
    return (
      extraDimensionError(positions)
      + distanceError(positions) 
      + chiralError(positions) 
    );
  }
};

template<bool compress>
struct errfGradient {
public:
  // Typedefs
  using Vector = dlib::matrix<double, 0, 1>;

  // State
  const unsigned N;
  const dlib::matrix<double, 0, 0> squaredBounds;
  const std::vector<ChiralityConstraint>& chiralityConstraints;

  // Constructor
  explicit errfGradient(
    const DistanceBoundsMatrix& distanceBounds,
    const std::vector<ChiralityConstraint>& chiralityConstraints
  ) : N(distanceBounds.N),
      squaredBounds(
        dlib::mat( // Converts Eigen::MatrixXd to dlib::matrix
          distanceBounds.makeSquaredBoundsMatrix()
        )
      ),
      chiralityConstraints(chiralityConstraints)
  {}

  // Information
  Vector referenceA(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(unsigned i = 0; i < N; i++) {
        if(i == alpha) { // skip identical indices
          continue;
        }

        const Vector diff = errfDetail::getPos(positions, alpha) - errfDetail::getPos(positions, i);

        const double upperTerm = dlib::length_squared(diff) / squaredBounds(
          std::min(i, alpha),
          std::max(i, alpha)
        ) - 1;

        if(upperTerm > 0) {
          dlib::set_rowm(
            gradient, 
            dlib::range(4 * alpha, 4 * alpha + 3)
          ) += 4 * diff * upperTerm / squaredBounds(
            std::min(i, alpha),
            std::max(i, alpha)
          );
        }
      }
    }

    return gradient;
  }

  Vector referenceB(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(unsigned i = 0; i < N; i++) {
        if(i == alpha) { // skip identical indices
          continue;
        }

        const Vector diff = errfDetail::getPos(positions, i) - errfDetail::getPos(positions, alpha);

        const double& lowerBoundSquared = squaredBounds(
          std::max(i, alpha),
          std::min(i, alpha)
        );

        const double quotient = lowerBoundSquared + dlib::length_squared(diff);

        const double lowerTerm = 2 * lowerBoundSquared / quotient - 1;

        if(lowerTerm > 0) {
          dlib::set_rowm(
            gradient,
            dlib::range(4 * alpha, 4 * alpha + 3)
          ) += 8 * lowerBoundSquared * diff * lowerTerm / (
            quotient * quotient
          );
        }
      }
    }

    return gradient;
  }

  Vector referenceC(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(const auto& constraint : chiralityConstraints) {
      const double volume = errfDetail::volume(positions, constraint.indices);

      const double upperTerm = volume - constraint.upper;
      if(upperTerm > 0) {
        // Set it for i -> constraint.indices[0]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[0], 4 * constraint.indices[0] + 2)
        ) += 2 * upperTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[2]) // k
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for j -> constraint.indices[1]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[1], 4 * constraint.indices[1] + 2)
        ) += 2 * upperTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[2]) // k
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for k -> constraint.indices[2]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[2], 4 * constraint.indices[2] + 2)
        ) += 2 * upperTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for l -> constraint.indices[3]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[3], 4 * constraint.indices[3] + 2)
        ) += 2 * upperTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[2]) // k
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[2]) // k
          )
        );
      }
    }

    return gradient;
  }

  Vector referenceD(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(const auto& constraint : chiralityConstraints) {
      const double volume = errfDetail::volume(positions, constraint.indices);

      const double lowerTerm = constraint.lower - volume;
      if(lowerTerm > 0) {
        // Set it for i -> constraint.indices[0]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[0], 4 * constraint.indices[0] + 2)
        ) -= 2 * lowerTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[2]) // k
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for j -> constraint.indices[1]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[1], 4 * constraint.indices[1] + 2)
        ) -= 2 * lowerTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[2]) // k
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for k -> constraint.indices[2]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[2], 4 * constraint.indices[2] + 2)
        ) -= 2 * lowerTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for l -> constraint.indices[3]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[3], 4 * constraint.indices[3] + 2)
        ) -= 2 * lowerTerm * (
          (
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[2]) // k
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[2]) // k
          )
        );

      }
    }

    return gradient;
  }

  Vector referenceE(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    if(compress) {
      for(unsigned i = 0; i < N; i++) {
        gradient(4 * i + 3) += 2 * positions(4 * i + 3);
      }
    }

    return gradient;
  }

  Vector reference(const Vector& positions) const {
    return (
      referenceA(positions)
      + referenceB(positions)
      + referenceC(positions)
      + referenceD(positions)
      + referenceE(positions)
    );
  }

  inline void gradientDistanceContribution(
    const Vector& positions,
    Vector& gradient,
    const double& lowerBoundSquared,
    const double& upperBoundSquared,
    const unsigned& alpha,
    const unsigned& i
  ) const {
    assert(lowerBoundSquared < upperBoundSquared);

    // For both
    const Vector alphaMinusI = (
      errfDetail::getPos(positions, alpha) 
      - errfDetail::getPos(positions, i)
    );

    const double diffLength = dlib::length_squared(alphaMinusI);

    // Upper term
    const double upperTerm = diffLength / upperBoundSquared - 1;

    assert(!( // Ensure early return logic is correct
      upperTerm > 0
      && ( // lowerTerm
        2 * lowerBoundSquared / (
          lowerBoundSquared + diffLength
        ) - 1
      ) > 0
    ));

    if(upperTerm > 0) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alpha, 4 * alpha + 3)
      ) += 4 * alphaMinusI * upperTerm / upperBoundSquared;

      /* if the upper term is set, there is no way the lower term needs to be
       * applied -> early return
       */
      return;
    }

    // Lower term
    const double quotient = lowerBoundSquared + diffLength;
    const double lowerTerm = 2 * lowerBoundSquared / quotient - 1;

    if(lowerTerm > 0) {
      /* We use -= because the lower term needs the position vector
       * difference (i - alpha), so we reuse alphaMinusI and just subtract 
       * from the gradient instead of adding to it
       */
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alpha, 4 * alpha + 3)
      ) -= 8 * lowerBoundSquared * alphaMinusI * lowerTerm / (
        quotient * quotient
      );
    }
  }

  // Operators
  Vector operator() (const Vector& positions) const {
    assert(positions.size() % 4 == 0);
    assert(positions.size() / 4 == N);

    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    // Distance gradient contributions (A and B)
    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(unsigned i = 0; i < alpha; i++) {
        // Here, i < alpha, so lower bound is (alpha, i), upper bound is (i, alpha)
        gradientDistanceContribution(
          positions,
          gradient,
          squaredBounds(alpha, i),
          squaredBounds(i, alpha),
          alpha,
          i
        );
      }

      for(unsigned i = alpha + 1; i < N; i++) {
        // Here, i > alpha, so lower bound is (i, alpha), upper bound is (alpha, i)
        gradientDistanceContribution(
          positions,
          gradient,
          squaredBounds(i, alpha),
          squaredBounds(alpha, i),
          alpha,
          i
        );
      }
    }

    // Chirality gradient contributions (C and D)
    for(const auto& constraint : chiralityConstraints) {
      const double volume = errfDetail::volume(positions, constraint.indices);

      // It may not occur that both terms contribute!
      assert(!(volume - constraint.upper > 0 && constraint.lower - volume > 0));

      const double factor = 2 * (
        std::max(0.0, volume - constraint.upper) 
        - std::max(0.0, constraint.lower - volume)
      );

      /* Make sure computing all cross products is worth it by checking that
       * one of both terms is actually greater than 0
       */
      if(factor != 0.0) {
        // Set for i -> constraint.indices[0]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[0], 4 * constraint.indices[0] + 2)
        ) += factor * (
          (
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[2]) // k
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for j -> constraint.indices[1]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[1], 4 * constraint.indices[1] + 2)
        ) += factor * (
          (
            errfDetail::getPos3D(positions, constraint.indices[2]) // k
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for k -> constraint.indices[2]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[2], 4 * constraint.indices[2] + 2)
        ) += factor * (
          (
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[3]) // l
          )
        );

        // Set for l -> constraint.indices[3]
        dlib::set_rowm(
          gradient,
          dlib::range(4 * constraint.indices[3], 4 * constraint.indices[3] + 2)
        ) += factor * (
          (
            errfDetail::getPos3D(positions, constraint.indices[1]) // j
            - errfDetail::getPos3D(positions, constraint.indices[2]) // k
          ).cross(
            errfDetail::getPos3D(positions, constraint.indices[0]) // i
            - errfDetail::getPos3D(positions, constraint.indices[2]) // k
          )
        );
      }
    }

    // Fourth dimension contribution (E)
    if(compress) {
      for(unsigned i = 0; i < N; i++) {
        gradient(4 * i + 3) += 2 * positions(4 * i + 3);
      }
    }

    return gradient;
  }
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
