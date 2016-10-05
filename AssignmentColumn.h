#ifndef INCLUDE_LIB_ASSIGNMENT_COLUMN_HPP
#define INCLUDE_LIB_ASSIGNMENT_COLUMN_HPP

#include <vector>
#include <cassert>

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

#endif
