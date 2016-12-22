#ifndef LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H
#define LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H

#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <set>
#include <map>
#include <cassert>

#include "UniqueAssignments/AssignmentColumn.h"

/* TODO
 */

/* NOTES
 * - The implementation of Assignment's virtual members cannot be put outside
 *   the class definition since templates may not be virtual
 */

namespace UniqueAssignments {

/*!
 * This class represents a simplified model of a sterically unique assignment
 * of a set of ligands to a stereocenter. It exists to uniquely identify the 
 * steric configuration at this stereocenter, and provides methods to assist
 * a systematic generation of all possible configurations. It is generalized 
 * over a number of symmetries which are encoded elsewhere and that serve as 
 * template parameter to this class.
 * \tparam Symmetry The template class specifying which symmetry is to be 
 *  enforced.
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
>
struct Assignment {
private:
/* Private member functions */
  std::vector<
    std::vector<unsigned>
  > _reduceGroups() const;

  bool _reducedGroupsAreEqual(
    const std::vector<
      std::vector<unsigned>
    >& a,
    const std::vector<
      std::vector<unsigned>
    >& b
  ) const; 

public:
/* Public members */
  std::vector<AssignmentColumn> positionOccupations;

  /* Constructors */
  Assignment() = delete;
  Assignment(
    const std::vector<char>& characters
  );
  Assignment(
    const std::vector<char>& characters,
    const std::vector<
      std::pair<unsigned, unsigned>
    >& pairedIndices
  );

  /* Modification */
  void sortOccupations() {
    std::sort(
      positionOccupations.begin(),
      positionOccupations.end()
    );
  }

  /*!
   * Generates the next permutation of positionOccupations in which the 
   * ligand connections are ordered.
   */
  bool nextPermutation() {
    while(std::next_permutation(
      positionOccupations.begin(),
      positionOccupations.end()
    )) {
      if(ligandConnectionsAreOrdered()) {
        return true;
      }
    }
    return false;
  };

  /*!
   * Applies a rotation function to the positionOccupations
   */
  void applyRotation(
    std::function<
      std::vector<AssignmentColumn>(
        const std::vector<AssignmentColumn>&
      )
    > rotationFunction
  ) {
    positionOccupations = rotationFunction(positionOccupations);
  }


  /* Information */
  /*!
   * Gets a map of ligand symbol character to position in the permutational 
   * symmetry.
   */
  std::map<
    char,
    std::vector<
      uint8_t
    >
  > getCharMap() const {
    std::map<
      char,
      std::vector<
        uint8_t
      >
    > returnMap;

    for(uint8_t i = 0; i < positionOccupations.size(); i++) {
      const char& columnChar = positionOccupations[i].character;
      if(returnMap.count(columnChar) == 0) {
        returnMap[columnChar] = {i};
      } else {
        returnMap[columnChar].push_back(i);
      }
    }

    return returnMap;
  }
  /*!
   * Checks whether the list of AssignmentColumns (positionOccupations) is
   * ordered.
   */
  bool occupationsAreOrdered() const {
    return std::is_sorted(
      positionOccupations.begin(),
      positionOccupations.end()
    );
  }

  /*!
   * Checks whether the ligand connections specified in the AssignmentColumns
   * are ordered. In a row wise view of the data in the columns:
   * 1010
   * 0101
   * is ordered, but
   * 0101
   * 1010
   * is not, although both sets of bit vectors have the same meaning.
   */
  bool ligandConnectionsAreOrdered() const {
    // shortcut if rowView will be empty
    if(positionOccupations[0].groups.size() == 0) return true;

    unsigned groupsRowSize = positionOccupations[0].groups.size();

    for(unsigned i = 0; i < groupsRowSize - 1; i++) {
      // compare row i with row i+1
      for(unsigned j = 0; j < Symmetry<>::size; j++) { 
        if(
          positionOccupations[j].groups[i] 
          < positionOccupations[j].groups[i + 1]
        ) {
          break;
        } else if(
          positionOccupations[j].groups[i] 
          > positionOccupations[j].groups[i + 1]
        ) {
          return false;
        }
      }
    }

    return true;
  }

  /*!
   * Checks whether this Assignment is rotationally superimposable with another.
   */
  bool isRotationallySuperimposable(
    const Assignment<Symmetry>& other
  ) const;

  /*!
   * Generates a set of all rotational equivalents of this Assignment as 
   * defined by its symmetry template parameter.
   */
  std::set<
    Assignment<Symmetry>
  > generateAllRotations() const;

  /*!
   * Implementation of the generation of a set of all rotational equivalents of
   * this Assignment as defined by its symmetry template parameter. Takes an 
   * interrupt callback as an argument to which it passes *this and a new 
   * rotational structure every time one is found. If the callback returns
   * true, the generation of assignments is terminated and a pair containing
   * the set of generated assignments and a boolean with the value true is 
   * returned. If the generation is allowed to finish, the full set and the 
   * boolean false are returned.
   */
  std::pair<
    std::set<
      Assignment<Symmetry>
    >,
    bool
  > _generateAllRotations(
    std::function<
      bool(const Assignment<Symmetry>&, const Assignment<Symmetry>&)
    > interruptCallbackOnNewAssignment
  ) const;

  //! Converts positionOccupations into a string for display
  std::string toString() const;

  bool reducedIsEqual(
    const Assignment<Symmetry>& other
  ) const;

  /* Operators */
  bool operator < (
    const Assignment<Symmetry>& other
  ) const;
  bool operator > (
    const Assignment<Symmetry>& other
  ) const;
  bool operator == (
    const Assignment<Symmetry>& other
  ) const;
  bool operator != (
    const Assignment<Symmetry>& other
  ) const;
};

#include "UniqueAssignments/Assignment.hxx"

} // eo namespace

#endif
