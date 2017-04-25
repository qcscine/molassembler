#ifndef INCLUDE_DG_GENERATE_CONFORMATION_H
#define INCLUDE_DG_GENERATE_CONFORMATION_H

#include "Molecule.h"
#include "Delib/PositionCollection.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/BFSConstraintCollector.h"
#include "Log.h"

#include <vector>
#include <Eigen/Core>

/* TODO
 */

namespace MoleculeManip {

namespace DistanceGeometry {

namespace detail {

Delib::PositionCollection convertToPositionCollection(
  const Eigen::VectorXd& vectorizedPositions,
  const EmbeddingOption& embedding
);

template<int segmentSize>
auto getEigen(
  const Eigen::VectorXd& v,
  const unsigned& index, 
  const unsigned& dimensionality
) {
  /* Return a fixed-size const reference to a part of the vector
   *
   * Interesting tidbit to the syntax below:
   * If you write: return v.segment<3>(3 * index);
   *             member function -^^^- integer 3
   *                               |
   *                          operator <
   * 
   * so to disambiguate, you write template before the name of the member
   * function.
   */
  return v.template segment<segmentSize>(dimensionality * index);
}

/* Functor to mimic a HOF that takes a distance Matrix as partial application
 * Using auto makePrototypePropagator(Matrix) {
 *   return [&matrix]()(Prototype) -> ChiralityConstraint {
 *     ...
 *   };
 * }
 * made it impossible to separate header and implementation, as the header is 
 * unusable since auto cannot be deduced.
 */
template<class DistanceGetter>
class Propagator {
private:
  DistanceGetter distanceGetter;

public:
  Propagator(DistanceGetter&& distanceGetter) 
    : distanceGetter(std::forward<DistanceGetter>(distanceGetter)) 
  {}

  ChiralityConstraint operator () (
    const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
  ) {
    using namespace Stereocenters;

    if(prototype.second == Stereocenter::ChiralityConstraintTarget::Flat) {
      return ChiralityConstraint {
        prototype.first,
        0.0
      };
    } else {
      /*
       * Calculate the volume. The target volume of the chirality constraint
       * created by the tetrahedron is calculated using internal coordinates (the
       * Cayley-Menger determinant), always leading to V > 0, so depending on the
       * current assignment, the sign of the result is switched.  The formula
       * used later in chirality constraint calculation for explicit coordinates
       * is adjusted by V' = 6 V to avoid an unnecessary factor, so we do that
       * here too:
       *               
       *    288 V²  = |...|               | substitute V' = 6 V
       * -> 8 (V')² = |...|               
       * ->      V' = sqrt(|...| / 8)
       *
       * where the Cayley-Menger determinant |...| is square symmetric:
       *   
       *          |   0    1    1    1    1  |
       *          |        0  d12² d13² d14² |
       *  |...| = |             0  d23² d24² |
       *          |                  0  d34² |
       *          |  ...                  0  |
       *
       */

      Eigen::Matrix<double, 5, 5> cayleyMenger;
      cayleyMenger.setZero();

      for(unsigned i = 0; i < 4; i++) {
        for(unsigned j = i + 1; j < 4; j++) {

          cayleyMenger(i + 1, j + 1) = pow(
            distanceGetter(
              prototype.first[i],
              prototype.first[j]
            ),
            2
          );
        }
      }

      // top row of cayleyMenger matrix
      for(unsigned i = 1; i < 5; i++) {
        cayleyMenger(0, i) = 1;
      }

#ifndef NDEBUG
      // Log in Debug builds, provided particular is set
      Log::log(Log::Particulars::PrototypePropagatorDebugInfo) 
        << "Cayley-menger matrix: " << std::endl << cayleyMenger << std::endl;
#endif

      auto determinant = static_cast<
        Eigen::Matrix<double, 5, 5>
      >(
        cayleyMenger.selfadjointView<Eigen::Upper>()
      ).determinant();


      /* If this fails, then we have not done distance selection properly. There
       * must be triangle inequality violations in the distances matrix.
       */
      assert(
        std::fabs(determinant) < 1e-8 // this is saying equal to zero
        || determinant > 0
      );

      if(std::fabs(determinant) < 1e-8) {
        // determinant is zero -> atoms involved are flat, this is alright
        return ChiralityConstraint {
          prototype.first,
          0
        };
      } else { // positive, greater than zero determinant
        auto chiralityTarget = sqrt(
          determinant / 8.0
        );

        if(
          prototype.second 
          == Stereocenters::Stereocenter::ChiralityConstraintTarget::Negative
        ) {
          chiralityTarget *= -1;
        }

        return ChiralityConstraint {
          prototype.first,
          chiralityTarget
        };
      }
    }

  }
};

struct DGStatistics {
  unsigned failures = 0;
  std::vector<double> energies;
};

struct DGResult {
  std::list<Delib::PositionCollection> ensemble;
  DGStatistics statistics;
};

// Generator to help with functional-style instantiation
template<class DistanceGetter> 
auto makePropagator(
  DistanceGetter&& distanceGetter
) {
  return Propagator<DistanceGetter>(
    std::forward<DistanceGetter>(distanceGetter)
  );
}

DGResult runDistanceGeometry(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization,
  const EmbeddingOption& embedding,
  const bool& useYInversionTrick = true,
  const BFSConstraintCollector::DistanceMethod& distanceMethod = BFSConstraintCollector::DistanceMethod::UFFLike
);

} // namespace detail

struct MoleculeDGInformation {
  DistanceBoundsMatrix distanceBounds;
  std::vector<
    Stereocenters::Stereocenter::ChiralityConstraintPrototype
  > chiralityConstraintPrototypes;

  MoleculeDGInformation(const unsigned& N);
};

MoleculeDGInformation gatherDGInformation(
  const Molecule& molecule,
  const BFSConstraintCollector::DistanceMethod& distanceMethod = BFSConstraintCollector::DistanceMethod::UFFLike
);

// Public functions
std::list<Delib::PositionCollection> generateEnsemble(
  const Molecule& molecule,
  const unsigned& numStructures,
  const MetrizationOption& metrization = MetrizationOption::off,
  const EmbeddingOption& embedding = EmbeddingOption::threeDimensional
);

Delib::PositionCollection generateConformation(
  const Molecule& molecule,
  const MetrizationOption& metrization = MetrizationOption::off,
  const EmbeddingOption& embedding = EmbeddingOption::threeDimensional
);

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
