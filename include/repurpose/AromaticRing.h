#ifndef INCLUDE_GRAPH_FEATURES_AROMATIC_RING_H
#define INCLUDE_GRAPH_FEATURES_AROMATIC_RING_H

#ifndef UNUSED
#define UNUSED(x) (void)(x)
#endif 

  /* TODO
   * - Test
   */

#include <math.h>

#include "Delib/ElementTypes.h" // Delib

#include "GraphFeature.h"
#include "CommonTrig.h"
#include "Molecule.h"

namespace MoleculeManip {

namespace GraphFeatures {

/* AromaticRing is deprecated as it should NOT be a GraphFeature */
class [[deprecated]] AromaticRing : public GraphFeature {
private:
  //std::vector<Delib::ElementType> _elementTypes;
  const std::vector<AtomIndexType> _members;
  const Molecule* _molPtr;
public:
  /* Constructor */
  AromaticRing() = delete;
  AromaticRing(
    const std::vector<AtomIndexType> members,
    Molecule* molPtr
  ) :
    _members(members),
    _molPtr(molPtr)
  {};

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
    /* NOTES
     * - Assumes _members is ordered! bonded atoms are next to one another
     */
    const double aromaticBondLengthVariation = 0.01; // in °A
    // inside angle sum of polygon divided by number of vertices
    const double insideAngleRadians = (
      180.0 - 360.0 / (double) _members.size()
    ) * ( M_PI / 180.0 );

    std::vector<DistanceConstraint> constraints;

    // for every conceivable pair of indices
    for(unsigned i = 0; i < _members.size(); i++) {
      for(unsigned j = i + 1; j < _members.size(); j++) {
        /* what's their distance? 
         *
         * For example:
         *
         * _members:
         * [] 0 1 2 3 4  size = 5
         *    4 8 2 9 1
         *
         * represents:
         *
         *     .9.   
         *   .°   °.
         *  1       2
         *   \     / 
         *    4 - 8
         *
         * so to find the distance between e.g. i = 0 and j = 4,
         * we go both ways around the ring from one of them and take the
         * chain that reached it first: m is the chain going backwards
         * (minus), p the chain going forwards (plus)
         *
         * [] 0 1 2 3 4  size = 5
         *    4 8 2 9 1
         *    ^
         *    (m,p)
         *
         * step: m -= 1 (wraps around), p += 1
         *
         * [] 0 1 2 3 4  size = 5
         *    4 8 2 9 1
         *      ^     ^
         *      p     m
         *
         * -> distance = 1, chain = {0, 4}
         *
         */
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
            ? _members.size() - 1
            : m - 1;
          mChain.push_back(_members[m]);
          p = (p == _members.size() - 1) 
            ? 0
            : p + 1;
          pChain.push_back(_members[p]);

          distance++;
        }

        std::vector<unsigned> chain = (m == j) 
          ? std::move(mChain)
          : std::move(pChain);

        if(distance == 1) {
          double bondDistance = Bond::calculateBondDistance(
            _molPtr -> getElementType(_members[i]),
            _molPtr -> getElementType(_members[j]),
            BondType::Aromatic
          );
          constraints.emplace_back(
            _members[i],
            _members[j],
            bondDistance - aromaticBondLengthVariation,
            bondDistance + aromaticBondLengthVariation
          );
        } else if(distance == 2) {
          // bond - angle - bond
          // use law of cosines
          double a = Bond::calculateBondDistance(
            _molPtr -> getElementType(chain[0]),
            _molPtr -> getElementType(chain[1]),
            BondType::Aromatic
          );
          double b = Bond::calculateBondDistance(
            _molPtr -> getElementType(chain[1]),
            _molPtr -> getElementType(chain[2]),
            BondType::Aromatic
          );

          constraints.emplace_back(
            _members[i],
            _members[j],
            CommonTrig::lawOfCosines(
              a - aromaticBondLengthVariation,
              b - aromaticBondLengthVariation,
              insideAngleRadians
            ),
            CommonTrig::lawOfCosines(
              a + aromaticBondLengthVariation,
              b + aromaticBondLengthVariation,
              insideAngleRadians
            )
          );
        } else if(distance == 3) {
          /* dihedral angles are strictly 0°
           * bond-angle-bond-angle-bond
           *
           *            d
           *   1 - - - - - - - - 4
           *    \        _ .·°  /
           *   a \   _.·°      / c
           *      \ α       α /
           *       2 ––––––– 3
           *            b
           *
           * 1-2: a
           * 2-3: b
           * 3-4: c
           * 2-4: x <- law of cosines (b, c, α)
           * angle β (3-2-4) <- law of sines (α, x, c)
           * 1-4: d <- law of cosines (a, x, α-β)
           *
           * Although this figure suggests it, b is not parallel to d!
           */

          const double a = Bond::calculateBondDistance(
            _molPtr -> getElementType(chain[0]),
            _molPtr -> getElementType(chain[1]),
            BondType::Aromatic
          );
          const double b = Bond::calculateBondDistance(
            _molPtr -> getElementType(chain[1]),
            _molPtr -> getElementType(chain[2]),
            BondType::Aromatic
          );
          const double c = Bond::calculateBondDistance(
            _molPtr -> getElementType(chain[2]),
            _molPtr -> getElementType(chain[3]),
            BondType::Aromatic
          );

          constraints.emplace_back(
            _members[i],
            _members[j],
            CommonTrig::getRingOneFourDistance(
              a - aromaticBondLengthVariation,
              b - aromaticBondLengthVariation,
              c - aromaticBondLengthVariation,
              insideAngleRadians
            ),
            CommonTrig::getRingOneFourDistance(
              a + aromaticBondLengthVariation,
              b + aromaticBondLengthVariation,
              c + aromaticBondLengthVariation,
              insideAngleRadians
            )
          );
        } // all larger distances are skipped
      } // end j loop
    } // end i loop

    return constraints;
  }

  virtual std::vector<ChiralityConstraint> chiralityConstraints() const override final {
    std::vector<ChiralityConstraint> constraints;
    // TODO continue here
    /* There are two types of chirality constraints arising here:
     * 1. Atoms in the ring -> These get constraints in sets of 4 atoms (ring 
     *  size must be > 3
     * 2. Ring substituents
     *
     *
     * Example:
     *
     *   1       2
     *    \     /
     *     5 = 6
     *    /     \
     *   3       4
     *   \\     //
     *     7 – 8
     *
     * Type 1 constraints: 3-5-6-4, 5-6-4-8, 6-4-8-7
     * Type 2 constraints: 1-5-6-3, 2-6-4-5
     */

    /* Type 1
     * Walk through the ring in groups of four, no overlap with initial atom 
     * required for coplanarity.
     */
    if(_members.size() > 3) {
      for(unsigned i = 0; i < _members.size() - 3; i++) {
        constraints.emplace_back(
          _members[i],
          _members[i + 1],
          _members[i + 2],
          _members[i + 3],
          0.0
        );
      }
    }

    /* Type 2
     * For every ring atom, check for non-ring substituents. If there is one 
     * and only one, create a chirality constraint for it.
     */
    for(unsigned i = 0; i < _members.size(); i++) {
      // get the adjacent atoms
      auto adjacent = _molPtr -> getBondedAtomIndices(_members[i]);
      // remove any ring atoms from adjacent
      adjacent.erase(
        std::remove_if(
          adjacent.begin(),
          adjacent.end(),
          [this](const auto& atomIndex) {
            return std::find(
              _members.begin(),
              _members.end(),
              atomIndex
            ) != _members.end();
          }
        ),
        adjacent.end()
      );
      if(adjacent.size() == 1) {
        constraints.emplace_back(
          adjacent[0], // the out-of-ring adjacent atom
          _members[i], // the bonded ring atom
          _members[ // the expression below selects the next ring atom
            (i + 1 == _members.size()) 
              ? 0
              : i + 1
          ],
          _members[ // the expression below selects the previous ring atom
            (i == 0)
              ? _members.size() - 1
              : i - 1
          ],
          0.0
        );
      }
    }

    return constraints;
  }

  virtual std::vector<Assignment> assignments() const override final{
    return {};
  }

  virtual bool hasAssignments() const override final{
    return false;
  }

  virtual bool assigned() const override final{
    return true;
  }

};

} // eo namespace GraphFeatures

} // eo namespace

#endif
