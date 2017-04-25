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
  if(nSites <= 1) throw std::logic_error(
    "Don't use a model on terminal atoms! Single bonds don't "
    "have stereochemistry!"
  );

  if(
    std::any_of(
      ligands.begin(),
      ligands.end(),
      [](const auto& ligand) {
        return std::get<2>(ligand).size() > 1;
      }
    )
  ) throw std::logic_error(
    "The ligand set includes ligands with multiple atoms on a site!"
    " This is just VSEPR, there should be no eta bond situations!"
  );
  
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

  if(E < 0) throw std::logic_error(
    "For some reason, E is < 0 in VSEPR. That shouldn't happen."
  );

  using Symmetry::Name;

  Name result;

  switch(X + E) {

    case 2: {
      result = Name::Linear;
    }; break;
    case 3: {
      if(X == 3) result = Name::TrigonalPlanar;
      else result = Name::Bent;
    }; break;
    case 4: {
      if(X == 4) result = Name::Tetrahedral;
      else if(X == 3) result = Name::TrigonalPyramidal;
      else result = Name::Bent;
    }; break;
    case 5: {
      if(X == 5) result = Name::TrigonalBiPyramidal;
      else if(X == 4) result = Name::Seesaw;
      else if(X == 3) result = Name::TShaped;
      else result = Name::Linear;
    }; break;
    case 6: {
      if(X == 6) result = Name::Octahedral;
      else if(X == 5) result = Name::SquarePyramidal;
      else result = Name::SquarePlanar;
    }; break;
    case 7: {
      if(X == 7) result = Name::PentagonalBiPyramidal;
      else if(X == 6) result = Name::PentagonalPyramidal;
      else result = Name::PentagonalPlanar;
    }; break;
    case 8: {
      result = Name::SquareAntiPrismatic;
    }; break;
    default: {
      std::stringstream ss;
      ss << "Could not find a fitting symmetry for your X + E case: "
        << "X = " << X << ", E = " << E << ". Maybe your molecular graph is "
        << " too weird for VSEPR. Have another look at it.";

      throw std::logic_error(
        ss.str().c_str()
      );
    }; break;
  }

  return result;
}

} // namespace LocalGeometry
