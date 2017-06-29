#include "VSEPR.h"
#include "Delib/ElementInfo.h"
#include <numeric>

/* TODO
 * - Changes needed to ensure correctness in complex environments:
 *   
 *   - concept of dative bonds is required
 *     -> How to detect if a bond is dative?
 * 
 * - Some of these exceptions would classify as "boneheaded" and may be better 
 *   off as asserts with a final cout message.
 */

namespace LocalGeometry {

Symmetry::Name VSEPR::determineGeometry(
  const Delib::ElementType& centerAtomType,
  const unsigned& nSites,
  const std::vector<LigandType>& ligands,
  const int& formalCharge
) {
  if(nSites <= 1) {
    throw std::logic_error(
      "Don't use a model on terminal atoms! Single bonds don't "
      "have stereochemistry!"
    );
  }

  if(
    std::any_of(
      ligands.begin(),
      ligands.end(),
      [](const auto& ligand) {
        return std::get<2>(ligand).size() > 1;
      }
    )
  ) {
    throw std::logic_error(
      "The ligand set includes ligands with multiple atoms on a site!"
      " This is just VSEPR, there should be no eta bond situations!"
    );
  }
  
  // get uncharged VE count
  auto VEOption = MoleculeManip::AtomInfo::mainGroupVE(centerAtomType);

  // return 
  if(!VEOption) {
    std::stringstream ss;
    ss << "You used VSEPR on a non-main group atom type: "
      << Delib::ElementInfo::symbol(centerAtomType)
      << "! That's not allowed.";
    throw std::logic_error(
      ss.str().c_str()
    );
  }

  // calculate X, E (VSEPR parameters)
  const unsigned& X = nSites;
  const int E = std::ceil(
    (
      static_cast<double>(VEOption.value()) 
      - formalCharge 
      - std::accumulate(
        ligands.begin(),
        ligands.end(),
        0.0,
        [](const double& carry, const auto& ligand) {
          /* can abort multiple ways:
           * vector front() is end() -> no ligands => API misuse
           * bondWeights has no entry for bty => error in bondWeights
           */
          return carry + bondWeights.at(
            std::get<2>(ligand).front().second
            // vec (pairs) ---^ pair -^ bty -^
          );
        }
      )
    ) / 2.0
  );

  if(E < 0) {
    throw std::logic_error(
      "For some reason, E is < 0 in VSEPR. That shouldn't happen."
    );
  }

  using Symmetry::Name;

  const auto XESum = X + E;

  if(XESum == 2) {
    return Name::Linear;
  }
  
  if(XESum == 3) {
    if(X == 3) {
      return Name::TrigonalPlanar;
    }

    return Name::Bent;
  }

  if(XESum == 4) {
    if(X == 4) {
      return Name::Tetrahedral;
    }

    if(X == 3) {
      return Name::TrigonalPyramidal;
    }

    return Name::Bent;
  }

  if(XESum == 5) {
    if(X == 5) {
      return Name::TrigonalBiPyramidal;
    }

    if(X == 4) {
      return Name::Seesaw;
    }

    if(X == 3) {
      return Name::TShaped;
    }

    return Name::Linear;
  }

  if(XESum == 6) {
    if(X == 6) {
      return Name::Octahedral;
    }

    if(X == 5) {
      return Name::SquarePyramidal;
    }

    return Name::SquarePlanar;
  }

  if(XESum == 7) {
    if(X == 7) {
      return Name::PentagonalBiPyramidal;
    }

    if(X == 6) {
      return Name::PentagonalPyramidal;
    }

    return Name::PentagonalPlanar;
  }

  if(XESum == 8) {
    return Name::SquareAntiPrismatic;
  }

  std::stringstream ss;
  ss << "Could not find a fitting symmetry for your X + E case: "
    << "X = " << X << ", E = " << E << ". Maybe your molecular graph is "
    << " too weird for VSEPR. Have another look at it.";

  throw std::logic_error(
    ss.str().c_str()
  );

  // Never runs, for static analysis
  return {};
}

} // namespace LocalGeometry
