#ifndef INCLUDE_LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_COLUMN_H
#define INCLUDE_LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_COLUMN_H

#include <vector>
#include <cassert>

namespace UniqueAssignments {

struct AssignmentColumn {
  char character;
  std::vector<bool> groups;

  AssignmentColumn() = delete;
  AssignmentColumn(
    const char& passCharacter,
    const std::vector<bool> passGroups
  ) : 
    character(passCharacter),
    groups(passGroups) 
  {};

  bool operator < (const AssignmentColumn& other) const;
  bool operator > (const AssignmentColumn& other) const;
  bool operator == (const AssignmentColumn& other) const;
  bool operator != (const AssignmentColumn& other) const;
};

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

} // eo namespace

#endif
