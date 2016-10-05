#ifndef LIB_ASSIGNMENT_HPP
#define LIB_ASSIGNMENT_HPP

#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <set>

#include "AssignmentColumn.h"

/* TODO
 * - reduceGroups needs serious testing
 * - make Assignment a template of the underlying SymmetryInformation derived 
 *   class, use its static functions to generalize the methods involved
 */

/* NOTES
 * - The implementation of Assignment's virtual members cannot be put outside
 *   the class definition since templates may not be virtual
 */

struct AbstractAssignment {
  /* public members */
  virtual void sortOccupations() = 0;
  virtual bool nextPermutation() = 0;
  virtual void applyRotation(
    std::function<
      std::vector<AssignmentColumn>(
        const std::vector<AssignmentColumn>&
      )
    > rotationFunction
  ) = 0;
  virtual bool occupationsAreOrdered() const = 0;
  virtual bool ligandConnectionsAreOrdered() const = 0;
  /* NOTE: the function below is needed, but the derived class should have 
   * two functions: one for const shared_ptr<Abstract>&, * one for const Derived&
   * solution to casting problem is in testing/cpp/inheritance_function_types
   */
  /*virtual bool isRotationallySuperimposable(
    const std::shared_ptr<AbstractAssignment>& other
  ) const = 0; */
};

template<
  template<typename T = AssignmentColumn>
  class Symmetry
>
struct Assignment : public AbstractAssignment {
private:
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
  virtual void sortOccupations() override final {
    std::sort(
      positionOccupations.begin(),
      positionOccupations.end()
    );
  }

  virtual bool nextPermutation() override final {
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
  virtual void applyRotation(
    std::function<
      std::vector<AssignmentColumn>(
        const std::vector<AssignmentColumn>&
      )
    > rotationFunction
  ) override final {
    positionOccupations = rotationFunction(positionOccupations);
  }


  /* Information */
  virtual bool occupationsAreOrdered() const override final {
    return std::is_sorted(
      positionOccupations.begin(),
      positionOccupations.end()
    );
  }

  virtual bool ligandConnectionsAreOrdered() const override final {
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

  bool isRotationallySuperimposable(
    const Assignment<Symmetry>& other
  ) const;

  /* Operators */
  bool operator < (
    const Assignment<Symmetry>& other
  ) const;
  bool operator == (
    const Assignment<Symmetry>& other
  ) const;
};

#include "Assignment.hxx"

#endif
