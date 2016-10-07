#ifndef UNUSED
#define UNUSED(x) (void)(x)
#endif 

/* TODO
 * - First finish Molecule interface. Accessors to private properties must 
 *   exist, it is impractical as hell to pass only the necessary information as 
 *   copies into every GraphFeature
 */

#include <math.h>

#include "ElementTypes.h" // Delib

#include "GraphFeature.h"
#include "CommonTrig.h"

namespace MoleculeManip {

class AromaticRing : public GraphFeature {
private:
  std::vector<Delib::ElementType> _elementTypes;
  std::vector<AtomIndexType> _members;
public:
  /* Modification */
  virtual void assign(const Assignment& assignment) override final {
    UNUSED(assignment);
  }

  /* Information */
  virtual std::string type() const override final {
    return "AromaticRing";
  }

  virtual std::set<AtomIndexType> involvedAtoms() const override final {
    return std::set<AtomIndexType>(
      _members.begin(),
      _members.end()
    );
  }

  virtual std::vector<DistanceConstraint> distanceConstraints() const override final {
    const double aromaticBondLengthVariation = 0.01; // in °A
    // inside angle sum of polygon divided by number of vertices
    const double insideAngle = 180.0 - 360.0 / (double) _members.size();

    std::vector<DistanceConstraint> constraints;

    // for every conceivable pair of indices
    for(unsigned i = 0; i < _members.size(); i++) {
      for(unsigned j = i + 1; j < _members.size(); j++) {
        // what's their distance? "DFS"
        unsigned m = i;
        unsigned p = i;
        std::vector<unsigned> mChain = {i};
        std::vector<unsigned> pChain = {i};

        unsigned distance = 0;
        while(
          m != j
          && p != j 
        ) {
          m = (m == 0) 
            ? _members.size() 
            : m - 1;
          mChain.push_back(m);
          p = (p == _members.size()) 
            ? 0
            : p + 1;
          pChain.push_back(p);

          distance++;
        }

        std::vector<unsigned> chain = (m == j) 
          ? std::move(mChain)
          : std::move(pChain);

        if(distance == 1) {
          double bondDistance = Bond::calculateBondDistance(
            _elementTypes[i],
            _elementTypes[j],
            BondType::Aromatic
          );
          constraints.push_back(
            std::make_tuple(
              _members[i],
              _members[j],
              bondDistance - aromaticBondLengthVariation,
              bondDistance + aromaticBondLengthVariation
            )
          );
        } else if(distance == 2) {
          // bond - angle - bond
          // use law of cosines
          double a = Bond::calculateBondDistance(
            _elementTypes[chain[0]],
            _elementTypes[chain[1]],
            BondType::Aromatic
          );
          double b = Bond::calculateBondDistance(
            _elementTypes[chain[1]],
            _elementTypes[chain[2]],
            BondType::Aromatic
          );

          constraints.push_back(
            std::make_tuple(
              _members[i],
              _members[j],
              CommonTrig::lawOfCosines(
                a - aromaticBondLengthVariation,
                b - aromaticBondLengthVariation,
                insideAngle
              ),
              CommonTrig::lawOfCosines(
                a + aromaticBondLengthVariation,
                b + aromaticBondLengthVariation,
                insideAngle
              )
            )
          );
        } else if(distance == 3) {
          // dihedral angles are strictly 0°
          // bond-angle-bond-angle-bond
        }
      }
    }

    

    return constraints;
  }

};

} // eo namespace
