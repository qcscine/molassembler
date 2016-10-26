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

#include "AssignmentColumn.hxx"

} // eo namespace

#endif
