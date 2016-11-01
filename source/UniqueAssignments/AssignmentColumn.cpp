#include "UniqueAssignments/AssignmentColumn.h"

namespace UniqueAssignments {

bool AssignmentColumn::operator < (const AssignmentColumn& other) const {
  assert(this -> groups.size() == other.groups.size());

  if(this -> character < other.character) return true;
  else if(this -> character > other.character) return false;
  else {
    for(unsigned i = 0; i < this -> groups.size(); i++) {
      if(this -> groups[i] < other.groups[i]) return true; 
      else if(this -> groups[i] > other.groups[i]) return false; 
      else continue;
    }

    return false;
  }
}

bool AssignmentColumn::operator > (const AssignmentColumn& other) const {
  assert(this -> groups.size() == other.groups.size());

  // TODO too lazy! inefficient
  return(!(
    *this < other
    || *this == other
  ));
}

bool AssignmentColumn::operator == (const AssignmentColumn& other) const {
  if(this -> character != other.character) return false;

  for(unsigned i = 0; i < this -> groups.size(); i++) {
    if(this -> groups[i] != other.groups[i]) return false;
  }

  return true;
}

bool AssignmentColumn::operator != (const AssignmentColumn& other) const {
  return !(
    (*this) == other
  );
}

} // eo namespace UniqueAssignments
