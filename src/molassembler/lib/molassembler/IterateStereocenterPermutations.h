#ifndef INCLUDE_ITERATE_STEREOCENTER_PERMUTATIONS_H
#define INCLUDE_ITERATE_STEREOCENTER_PERMUTATIONS_H

#include "Molecule.h"
#include "StdlibTypeAlgorithms.h"

#include "temple/Containers.h"

/*! @file
 *
 * Contains functionality in order to systematically vary stereocenter
 * assignments in order to exhaustively enumerate all possible conformers.
 */

/* TODO
 * - Maybe make a set of logical requirements that the iterators must fulfill to
 *   test against, similar for the comparison operators? e.g. self-consistency
 *   a == a, some of the concept requirements specified in the STL, etc.
 * - As soon as reranking is implemented, this code can generate extraneous
 *   combinations of stereocenter assignment. Since a stereocenter's assignment
 *   can effect ranking, the number of assignments of ANY stereocenter could
 *   change. This can fail spectacularly by trying to set an assignment on some
 *   stereocenter that no longer exists. See meso forms.
 *   Rewriting this so that it tracks which combinations it has examined by
 *   keeping a full mapping indices of to assignments might be possible, or
 *   maybe a full tree exploration starting from each stereocenter once.
 * - Move this to be a part of StereocenterList
 */

namespace molassembler {

/*!
 * Range-for temporary object to allow iterator-like behavior in permutation of
 * stereocenters
 */
/*class StereocenterPermutationTemporary {
private:
  const Molecule _molecule;

public:
  StereocenterPermutationTemporary(
    const Molecule& molecule
  ) : _molecule(molecule) 
  {}

  using BaseIteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    Molecule,                  // value_type
    unsigned,                  // difference_type
    const Molecule*,           // pointer
    const Molecule             // reference
  >;

  template<typename PointerType> 
  class iterator : public BaseIteratorType {
  private:
    mutable Molecule _molecule;
    std::vector<unsigned> _currentCombination;
    std::vector<unsigned> _limits;
    std::vector<
      std::shared_ptr<Stereocenters::Stereocenter>
    > _stereocenterList;
    bool _done;

  public:
    explicit iterator(
      const Molecule& molecule,
      const bool& endPosition
    ) : _molecule(molecule),
        _done(endPosition)
    {
      if(!endPosition) {
        for(const auto& stereocenterPtr : _molecule.getStereocenterList()) {

          if(stereocenterPtr -> numStereopermutations() > 1) {
            _stereocenterList.push_back(stereocenterPtr);
            _limits.push_back(
              stereocenterPtr -> numStereopermutations() - 1
              // limits are 0-based, numStereopermutations is 1-based
            );
          }
        }

        // appropriately initialize combination
        _currentCombination = std::vector<unsigned>(_limits.size(), 0);
      }
    }

    iterator& operator ++ () { 
      if(!_done) {
        _done = !StdlibTypeAlgorithms::nextCombinationPermutation(
          _currentCombination,
          _limits
        );
      }
      return *this; 
    }

    iterator operator++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval; 
    }

    bool operator == (iterator other) const {
      if(
        _done == other._done
        && _done
      ) {
        return true;
      } else {
        return _currentCombination == other._currentCombination;
      }
    }

    bool operator != (iterator other) const {
      return !(
        *this == other
      );
    }

    std::string toString() const {
      return (
        "combination: vec{"s
        + temple::condenseIterable(_currentCombination)
        + "}, limits: vec{"s
        + temple::condenseIterable(_limits)
        + "}, done: "s + std::to_string(_done)
      );
    }

    typename BaseIteratorType::reference operator * () const { 
      for(unsigned i = 0; i < _currentCombination.size(); ++i) {
        if( // only reassign if unset or at different assignment
          !_stereocenterList.at(i) -> assigned()
          || (
            _stereocenterList.at(i) -> assigned().value() 
            != _currentCombination.at(i)
          )
        ) {
          _stereocenterList.at(i) -> assign(
            _currentCombination.at(i)
          );
        }
      }

      return _molecule;
    }
  };

  iterator<const Molecule*> begin() {
    return iterator<const Molecule*>(
      _molecule,
      false
    );
  }

  iterator<const Molecule*> end() {
    return iterator<const Molecule*>(
      _molecule,
      true
    );
  }
};

StereocenterPermutationTemporary iterateStereocenterPermutations(
  const Molecule& molecule
) {
  return StereocenterPermutationTemporary(molecule);
}*/

} // namespace molassembler

#endif
