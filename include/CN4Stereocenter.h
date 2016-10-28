#ifndef INCLUDE_CN4_STEREOCENTER_H
#define INCLUDE_CN4_STEREOCENTER_H

#include <math.h>
#include <experimental/optional>
#include "ElementTypes.h"

#include "Stereocenter.h"
#include "CommonTrig.h"
#include "Molecule.h"
#include "Cache.h"
#include "StdlibTypeAlgorithms.h"
#include "UniqueAssignments/GenerateUniques.h"
#include "UniqueAssignments/SymmetryInformation.h"

namespace MoleculeManip {

namespace Stereocenters {

class CN4Stereocenter : public Stereocenter {
private:
/* Typedefs */
  using AssignmentType = UniqueAssignments::Assignment<
    PermSymmetry::Tetrahedral
  >;

/* Private data */
  const Molecule* _molPtr;
  const AtomIndexType _centerAtom;

  // State of whether it is assigned, and if so, in which (of _cache["assignmentlist"])
  std::experimental::optional<unsigned> _assignment;

  /* mutable so as to allow const qualified member functions to read and edit
   * the cache
   */
  mutable Cache _cache {
    std::make_pair(
      "assignmentsList",
      [this]() {
        return UniqueAssignments::uniqueAssignments(
          AssignmentType(
            _reduceSubstituents()
          )
        );
      }
    )
  };

/* Private member functions */
  std::vector<char> _reduceSubstituents() const {
    // OLD ––––––––––––––––––
    std::vector<AtomIndexType> rankedSubstituentNextAtomIndices;
    std::set<
      std::pair<
        AtomIndexType,
        AtomIndexType
      >
    > equalSubstituentPairsSet;

    std::tie(
      rankedSubstituentNextAtomIndices,
      equalSubstituentPairsSet
    ) = _molPtr -> rankCIPPriority(
      _centerAtom,
      {} // nothing to exclude
    );

    /* C++17 structured binding improvement assuming tuples are still
     * convertible to pairs
    auto [
      rankedSubstituentNextAtomIndices,
      equalSubstituentPairsSet
    ] = _molPtr -> rankCIPPriority(
      _centerAtom,
      {}
    );
    */

    /* restructure the set of pairs to a vector of non-overlapping sets
     * e.g. set{pair{1, 3}, pair{1, 4}, pair{2, 5}} 
     *  -> vector{set{1, 3, 4}, set{2, 5}}
     */
    auto setsVector = StdlibTypeAlgorithms::makeIndividualSets(
      equalSubstituentPairsSet
    );

    // Add lone substituents to setsVector
    for(const auto& index: rankedSubstituentNextAtomIndices) {
      // if the current substituent index is not in any of the sets
      if(!std::accumulate(
        setsVector.begin(),
        setsVector.end(),
        false,
        [&index](const bool& carry, const std::set<AtomIndexType>& set) {
          return (
            carry
            || set.count(index) == 1
          );
        }
      )) {
        // add a single-atom set
        setsVector.push_back(
          std::move(
            std::set<AtomIndexType>{index}
          )
        );
      }
    }

    /* so now we have e.g.
     * setsVector = vector{set{1, 4}, set{2}, set{3}};
     * rankedSubstituentNextAtomIndices = vector{ 2, 1, 4, 3};
     * -> reduce to {A, B, B, C}
     */

    // create a mapping between indices and ligand symbols
    std::map<
      AtomIndexType,
      char
    > indexSymbolMap;

    const char initialChar = 'A';
    for(const auto& index : rankedSubstituentNextAtomIndices) {
      // find position in setsVector
      unsigned posInSetsVector = 0;
      while(
        setsVector[posInSetsVector].count(index) == 0 
        && posInSetsVector < setsVector.size()
      ) {
        posInSetsVector++;
      }
      indexSymbolMap[index] = initialChar + posInSetsVector;
    }

    std::vector<char> ligandSymbols (rankedSubstituentNextAtomIndices.size());
    std::transform(
      rankedSubstituentNextAtomIndices.begin(),
      rankedSubstituentNextAtomIndices.end(),
      ligandSymbols.begin(),
      [&indexSymbolMap](const auto& index) {
        return indexSymbolMap.at(index);
      }
    );

    return ligandSymbols;

    /* TODO no use of connectivity information as of yet to determine whether 
     * ligands are bridged!
     */
  }

public:
/* Public member functions */
  /* Construction */
  CN4Stereocenter(
    const Molecule* molPtr,
    const AtomIndexType& center
  ) : 
    _molPtr(molPtr),
    _centerAtom (center)
  {};

  /* Modification */
  /*!
   * Assign this feature
   */
  virtual void assign(const unsigned& assignment) override final {
    assert(
      assignment 
      < _cache.getOrGenerate<
        std::vector<AssignmentType>
      >("assignmentsList").size()
    );
    _assignment = assignment;
  }

  /* Information */
  /*!
   * Return a string specifying the type of stereocenter
   */
  virtual std::string type() const override final {
    return "CN4Stereocenter";
  }

  /*!
   * Return a set of involved atom indices
   */
  virtual std::set<AtomIndexType> involvedAtoms() const override final {
    return std::set<AtomIndexType>({_centerAtom});
  }

  /*!
   * Return a list of distance constraints
   */
  virtual std::vector<DistanceConstraint> distanceConstraints() const override final {
    std::vector<DistanceConstraint> constraints;
    const auto neighbors = _molPtr -> getBondedAtomIndices(_centerAtom);
    for(unsigned i = 0; i < neighbors.size(); i++) {
      for(unsigned j = i + 1; j < neighbors.size(); j++) {

      }
    }
  }

  /*!
   * Return a list of chirality constraints
   */
  virtual std::vector<ChiralityConstraint> chiralityConstraints() const = 0;

  /*!
   * Return the number of possible assignments at this feature
   */
  virtual unsigned assignments() const override final {
    return _cache.getOrGenerate<
      std::vector<AssignmentType>
    >("assignmentsList").size();
  }

  /*!
   * Return whether this feature has been assigned or not
   */
  virtual bool assigned() const override final {
    // Fundamentals TS
    return (bool) _assignment;
    // C++17
    // return _assignment.has_value();
  }
};

}

}

#endif
