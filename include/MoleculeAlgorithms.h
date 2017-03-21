#ifndef INCLUDE_MOLECULE_ALGORITHMS_H
#define INCLUDE_MOLECULE_ALGORITHMS_H

// Delib
#include "Types/PositionCollection.h"

#include "Molecule.h"

// Implementation
#include "CNStereocenter.h"
#include "template_magic/templateMagic.h"
#include "CommonTrig.h"

namespace MoleculeManip {

namespace detail {
  double getAngle(
    const Delib::PositionCollection& positions,
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) {
    auto a = positions[i].asEigenVector() - positions[j].asEigenVector(),
         b = positions[k].asEigenVector() - positions[j].asEigenVector();

    return std::acos(
      a.dot(b) / (
        a.norm() * b.norm()
      )
    );
  }

  double getDistance(
    const Delib::PositionCollection& positions,
    const AtomIndexType& i,
    const AtomIndexType& j
  ) {
    return (
      positions[i].asEigenVector()
      - positions[j].asEigenVector()
    ).norm();
  }

  double toRadians(const double& inDegrees) {
    return M_PI * inDegrees / 180;
  }
}

Molecule fitToPositions(
  Molecule&& molecule,
  const Delib::PositionCollection& positions
) {

  const auto& adjacencies = molecule.getAdjacencyList();

  for(auto& stereocenterPtr : molecule.stereocenters) {

    // How to find out which Assignment best reflects the read positions?
    if(stereocenterPtr -> type() == Stereocenters::Type::CNStereocenter) {
      /* STEPS
       * - Reduce the local geometry to some combination of 
       *   internal coordinates, 1-3 distances and signed tetrahedron values
       *   (chirality constraints) so that it can be compared easily to gathered
       *   distance geometry constraints
       * - Cycle through all Symmetries of appropriate size, checking the 
       *   gathered 3D data against supposed DG constraints of the candidate
       *   Symmetry and assignment (if present)
       * - Cycle through the Stereocenter's assignments and calculate total
       *   deviation for all angles
       * - Pick the assignment that has lowest deviation, if it has multiplicity 1
       */
      // Downcast the Stereocenter to a CNStereocenter
      auto CNStereocenterPtr = std::dynamic_pointer_cast<
        Stereocenters::CNStereocenter
      >(stereocenterPtr);

      std::map<
        Symmetry::Name, 
        std::vector<double> // Deviation for every assignment in symmetry
      > symmetryDeviations;

      for(const auto& symmetryName : Symmetry::allNames) {
        if( // Skip any Symmetries of different size
          Symmetry::size(symmetryName) != Symmetry::size(
            CNStereocenterPtr -> symmetry
          )
        ) continue;

        // Change the symmetry of the CNStereocenter
        CNStereocenterPtr -> changeSymmetry(symmetryName);

        // Get the adjacent indices to the central atom
        std::vector<AtomIndexType> adjacentAtoms = adjacencies.getAdjacencies(
          CNStereocenterPtr -> centerAtom
        );

        std::vector<double> deviations;

        for(
          unsigned assignment = 0;
          assignment < (CNStereocenterPtr -> assignments());
          assignment++
        ) {
          // Assign the stereocenter
          CNStereocenterPtr -> assign(assignment);


          // Calculate the deviation from the positions
          /* Three contributions to deviations
           * - i-j-k angles
           * - 1-3 distances (via 1-2 bond distances from positions and 
           *   symmetry-ideal angles)
           * - chirality constraints (if applicable)
           */
          // i-j-k angles
          double angleDeviation = TemplateMagic::sum(
            TemplateMagic::allPairsMap(
              adjacentAtoms,
              [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
                return std::fabs(
                  detail::getAngle( // The angle from the positions
                    positions,
                    i,
                    CNStereocenterPtr -> centerAtom,
                    k
                  ) - detail::toRadians( // The ideal angle from the Stereocenter
                    CNStereocenterPtr -> angle(
                      i,
                      CNStereocenterPtr -> centerAtom,
                      k
                    )
                  )
                );
              }
            )
          );
          
          // 1-3 distances
          double oneThreeDeviation = TemplateMagic::sum(
            TemplateMagic::allPairsMap(
              adjacentAtoms,
              [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
                return std::fabs(
                  detail::getDistance( // i-k 1-3 distance from positions
                    positions,
                    i,
                    k
                  ) - CommonTrig::lawOfCosines( // idealized 1-3 distance from
                    detail::getDistance( // i-j 1-2 distance from positions
                      positions,
                      i,
                      CNStereocenterPtr -> centerAtom
                    ),
                    detail::getDistance( // j-k 1-2 distance from positions
                      positions,
                      CNStereocenterPtr -> centerAtom,
                      k
                    ),
                    detail::toRadians(
                      CNStereocenterPtr -> angle( // idealized Stereocenter angle
                        i,
                        CNStereocenterPtr -> centerAtom,
                        k
                      )
                    )
                  )
                );
              }
            )
          );

          double chiralityConstraintDeviation = 0; 
          /* TODO as soon as we know how chirality constraints arise from
           * stereocenters
           */

          // Emit structured data
          std::cout << CNStereocenterPtr -> centerAtom
            << ", " << (
              std::find(
                Symmetry::allNames.begin(),
                Symmetry::allNames.end(),
                CNStereocenterPtr -> symmetry
              ) - Symmetry::allNames.begin()
            ) << ", " << assignment
            << ", " << CNStereocenterPtr -> assignments() << ", "
            << std::setprecision(4) << std::fixed
            << angleDeviation << ", "
            << oneThreeDeviation << ", "
            << chiralityConstraintDeviation 
            << std::endl;

          // Add to deviations vector
          deviations.emplace_back(
            angleDeviation 
            + oneThreeDeviation 
            + chiralityConstraintDeviation
          );
        }
        
        // Add it to the deviations
        symmetryDeviations[symmetryName] = deviations;
      }

      // Group analyze the assignments, assuming equal deviations do not happen
      // for different symmetries -> dangerous!
      std::map<
        std::pair<
          double, // deviation
          Symmetry::Name // Symmetry
        >, 
        std::vector<unsigned> // which assignments
      > groups;

      for(const auto& iterPair : symmetryDeviations) {
        for(unsigned i = 0; i < iterPair.second.size(); i++) {
          auto& deviation = iterPair.second[i];

          // If key doesn't exist yet
          if(groups.count(
            std::make_pair(deviation, iterPair.first)
          ) == 0) {
            groups.insert(
              std::make_pair(
                std::make_pair(deviation, iterPair.first),
                std::vector<unsigned>({i})
              )
            );
          } else {
            groups.at(
              std::make_pair(deviation, iterPair.first)
            ).push_back(i);
          }
        }
      }

      /*for(const auto& iterPair : groups) {
        std::cout << "(" << iterPair.first.first << ", " << Symmetry::name(iterPair.first.second) << ") => vec{";
        for(const auto& value: iterPair.second) {
          std::cout << value;
          if(value != iterPair.second.back()) std::cout << ", ";
        }
        std::cout << "}" << std::endl;
      }*/

      // Find Symmetry with lowest deviation -> initialize to SOMETHING
      auto lowestDeviation = (groups.begin() -> first).first;
      auto bestSymmetry = (groups.begin() -> first).second;

      for(const auto& iterPair : groups) {
        if(iterPair.first.first < lowestDeviation) {
          lowestDeviation = iterPair.first.first;
          bestSymmetry = iterPair.first.second;
        }
      }


      // Set it to the best symmetry
      CNStereocenterPtr -> changeSymmetry(bestSymmetry);

      /* If that pair with lowest deviation has only a single assignment,
       * assign it, else leave it unassigned
       */
      if(groups.at({lowestDeviation, bestSymmetry}).size() == 1) {
        CNStereocenterPtr -> assign(
          groups.at({lowestDeviation, bestSymmetry}).front()
        );
      } /*else { // Debug message
        unsigned countGroupsInBestSymmetry = 0;
        for(const auto& iterPair : groups) {
          if(iterPair.first.second == bestSymmetry) {
            countGroupsInBestSymmetry += 1;
          }
        }
        std::cout << "Best symmetry has " << countGroupsInBestSymmetry 
          << " deviation group(s)." << std::endl;
      }*/
      
    } else if(stereocenterPtr -> type() == Stereocenters::Type::EZStereocenter) {
      /* STEPS
       * - Calculate dihedral angle of high-priority pair from 3D
       *   -> Select E/Z within tolerance of 0° / 180° endpoints
       *   -> Throw outside of those tolerances
       */
    }
  } // end for every stereocenter

  return molecule;
}

} // eo namespace MoleculeManip

#endif
