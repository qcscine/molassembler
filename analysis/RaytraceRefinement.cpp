#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"

#include "BoundsFromSymmetry.h"
#include "IO.h"
#include "Log.h"

#include "AnalysisHelpers.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

int main() {
  Log::level = Log::Level::None;
  Log::particulars = {Log::Particulars::DGRefinementChiralityNumericalDebugInfo};

  /* Re-implementation of DG procedure so that we can use step() in the 
   * conjugated descent gradient solver
   */

  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGConjugatedGradientDescentSolver;

  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  // TODO this will need adjustment when some experience exists
  stopCriteria.iterations = 1000; 
  stopCriteria.fDelta = 1e-5;

  DGConjugatedGradientDescentSolver.setStopCriteria(stopCriteria);

  const unsigned nStructures = 10;

  for(const auto& symmetryName : Symmetry::allNames) {
    // Make a space-free string from the name
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Generate distance bounds
    auto simpleMol = DGDBM::symmetricMolecule(symmetryName);

    // Begin
    const auto DGData = gatherDGInformation(simpleMol);
    unsigned optimizationFailures = 0;
    const double failureRatio = 0.5;

    for(
      unsigned currentStructureNumber = 0;
      // Failed optimizations do not count towards successful completion
      currentStructureNumber - optimizationFailures < nStructures;
      currentStructureNumber += 1
    ) {

      const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix(
        MetrizationOption::off
      );
      
      /* Get the chirality constraints by converting the prototypes found by the
       * collector into full chiralityConstraints using the distances matrix
       */
      auto chiralityConstraints = TemplateMagic::map(
        DGData.chiralityConstraintPrototypes,
        /* Partial application with distances matrix so we have a unary function
         * to perform the mapping with
         */
        detail::makePropagator(
          [&distancesMatrix](const AtomIndexType& i, const AtomIndexType& j) {
            return distancesMatrix(
              std::min(i, j),
              std::max(i, j)
            );
          }
        )
      );

      /* Instantiantiate the refinement problem and its solver, set the stop 
       * criteria
       */
      DGRefinementProblem<double> problem(
        chiralityConstraints,
        DGData.distanceBounds
      );

      // Make a metric matrix from the distances matrix
      MetricMatrix metric(distancesMatrix);

      // Get a position matrix by embedding the metric matrix
      auto embeddedPositions = metric.embed(EmbeddingOption::threeDimensional);

      /* If a count of chirality constraints reveals that more than half are
       * incorrect, we can invert the structure (by multiplying e.g. all y 
       * coordinates with -1) and then have more than half of chirality 
       * constraints correct! In the count, chirality constraints with a target
       * value of zero are not considered (this would skew the count as those
       * chirality constraints should not have to pass an energetic maximum to 
       * converge properly as opposed to tetrahedra with volume).
       */
      if(detail::moreThanHalfChiralityConstraintsIncorrect(
        embeddedPositions,
        chiralityConstraints
      )) {
        embeddedPositions.row(2) *= -1;
      }

      // Vectorize the positions for use with cppoptlib
      Eigen::VectorXd vectorizedPositions {
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      };

      // Run the minimization, but step-wise!
      auto stepResult = DGConjugatedGradientDescentSolver.step(problem, vectorizedPositions);
      unsigned iterations = 1;

      writePOVRayFile(
        spaceFreeName,
        currentStructureNumber,
        iterations,
        problem,
        simpleMol,
        vectorizedPositions,
        stepResult
      );

      while(stepResult.status == cppoptlib::Status::Continue) {
        stepResult = DGConjugatedGradientDescentSolver.step(problem, vectorizedPositions);
        iterations += 1;

        writePOVRayFile(
          spaceFreeName,
          currentStructureNumber,
          iterations,
          problem,
          simpleMol,
          vectorizedPositions,
          stepResult
        );
      }

      std::cout << spaceFreeName << "-" << currentStructureNumber << ": "
        << iterations << " iterations, exited with code " << static_cast<int>(stepResult.status) << std::endl;

      // What to do if the optimization fails
      if(DGConjugatedGradientDescentSolver.status() == cppoptlib::Status::Continue) {
        optimizationFailures += 1;

        if(
          static_cast<double>(optimizationFailures) / currentStructureNumber 
          >= failureRatio
        ) {
          throw std::runtime_error("Refinement failures exceeded threshold!");
        }
      } 
    }

    std::cout << optimizationFailures << " failures." << std::endl;
  }
}
