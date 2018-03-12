#include "LocalGeometryModel.h"
#include "Delib/ElementInfo.h"
#include "AtomInfo.h"

namespace molassembler {

namespace LocalGeometry {

const std::map<BondType, double> bondWeights {
  {BondType::Single, 1.0},
  {BondType::Double, 2.0},
  {BondType::Triple, 3.0},
  {BondType::Quadruple, 4.0},
  {BondType::Quintuple, 5.0},
  {BondType::Sextuple, 6.0},
  {BondType::Aromatic, 1.5},
  {BondType::Eta, 0.0} // TODO is this wise? (duplicate in BondDistance!)
};

boost::optional<Symmetry::Name> vsepr(
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

  if(!AtomInfo::isMainGroupElement(centerAtomType)) {
    return boost::none;
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
    return boost::none;
    /*"The ligand set includes ligands with multiple atoms on a site!"
    " This is just VSEPR, there should be no eta bond situations!"*/
  }
  
  // get uncharged VE count
  auto VEOption = molassembler::AtomInfo::mainGroupVE(centerAtomType);

  // return 
  if(!VEOption) {
    return boost::none;
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
    // "For some reason, E is < 0 in VSEPR. That shouldn't happen."
    return boost::none;
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

  /* "Could not find a fitting symmetry for your X + E case: "
    << "X = " << X << ", E = " << E << ". Maybe your molecular graph is "
    << " too weird for VSEPR. Have another look at it.";*/

  // Never runs, for static analysis
  return boost::none;
}

boost::optional<Symmetry::Name> firstOfSize(const unsigned& size) {
  // Pick the first Symmetry of fitting size
  auto findIter = std::find_if(
    Symmetry::allNames.begin(),
    Symmetry::allNames.end(),
    [&size](const auto& symmetryName) -> bool {
      return Symmetry::size(symmetryName) == size;
    }
  );

  if(findIter == Symmetry::allNames.end()) {
    return boost::none;
  }

  return *findIter;
}

} // namespace LocalGeometry

} // namespace molassembler
