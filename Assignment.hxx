/* Public members */
/*  Constructors */
/*!
 * Constructs an Assignment from a list of ligand characters.
 * \tparam Symmetry A SymmetryInformation derived class template.
 * \param characters A vector of chars signifying abstract ligands.
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
>
Assignment<Symmetry>::Assignment(
  const std::vector<char>& characters
) {
  assert(characters.size() == Symmetry<>::size);

  for(unsigned i = 0; i < Symmetry<>::size; i++) {
    positionOccupations.emplace_back(
      characters[i],
      std::vector<bool>()
    );
  }

  sortOccupations();
}


/*!
 * Construct an Assignment from a list of ligand characters and a list of 
 * bonded indices referencing the ligand characters.
 * \tparam Symmetry A SymmetryInformation derived class template.
 * \param characters A vector of chars signifying abstract ligands.
 * \param pairedIndices A vector of pairs. Describes which ligand characters 
 *  are bonded to one another.
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
>
Assignment<Symmetry>::Assignment(
  const std::vector<char>& characters,
  const std::vector<
    std::pair<unsigned, unsigned>
  >& pairedIndices
) {
  // make sure the number of characters matches the current symmetry
  assert(characters.size() == Symmetry<>::size);

  // create a row representation of the groups
  std::vector<
    std::vector<bool>
  > groups;

  for(const auto& indexPair : pairedIndices) {
    // make sure the indices are valid
    assert(
      indexPair.first < Symmetry<>::size
      && indexPair.second < Symmetry<>::size
      && indexPair.first != indexPair.second
    );

    // create a group from the pair
    std::vector<bool> group (Symmetry<>::size, false);
    group.at(indexPair.first) = true;
    group.at(indexPair.second) = true;

    // add the group 
    groups.push_back(group);
  }

  // create Columns
  for(unsigned i = 0; i < Symmetry<>::size; i++) {
    // create the column's group bit vector
    std::vector<bool> columnGroupBitVector;
    for(const auto& groupRow : groups) {
      columnGroupBitVector.push_back(groupRow.at(i));
    }
    positionOccupations.emplace_back(
      characters[i],
      columnGroupBitVector
    );
  }

  // sort the positionOccupations
  sortOccupations();
}

/* Public members */
/*!
 * Determines if two Assignment instances are rotationally superimposable. This
 * function exploits the rotations defined in the template parameter Symmetry
 * to check equivalence.
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
> bool Assignment<Symmetry>::isRotationallySuperimposable(
  const Assignment<Symmetry>& other
) const {

  // add the initial structure to a set of Assignments
  std::set<Assignment> enumeratedAssignments = {other};   

  // Systematically explore all rotations
  // maximum element is the size of the rotation vector
  unsigned linkLimit = Symmetry<>::rotations.size();

  // initialize 
  std::vector<unsigned> chain = {0};
  std::vector<
    Assignment<Symmetry>
  > chainStructures = {other};
  unsigned depth = 0;

  // begin loop
  while(chain.at(0) < linkLimit) {
    // TEMP
    /*std::cout << "chain: ";
    for(const auto& link : chain) {
      std::cout << link << " ";
    }
    std::cout << std::endl;*/

    // perform rotation
    // copy the last element in chainStructures
    Assignment<Symmetry> generated = chainStructures.at(
      chainStructures.size() - 1
    );
    // apply the rotation referenced by the last link in chain
    generated.applyRotation(
      Symmetry<>::rotations[
        chain.at(
          chain.size() -1
        )
      ].first
    );

    // is it something new?
    if(enumeratedAssignments.count(generated) == 0) {
      // is it the same as this?
      if(*this == generated) {
        return true;
      }
      // add it to the set
      enumeratedAssignments.insert(generated);
      // add it to chainStructures
      chainStructures.push_back(generated);
      // increase depth, add a link
      depth++;
      chain.emplace_back(0);
    } else {
      // if we are not at the maximum instruction
      if(chain.at(depth) < linkLimit - 1) {
        chain.at(depth)++;
      } else {
        // collapse the chain until we are at an incrementable position
        while(
          depth > 0
          && chain.at(depth) == linkLimit - 1
        ) {
          chain.pop_back();
          chainStructures.pop_back();
          depth--;
        }

        chain.at(depth)++;
      }
    }
  }

  /* we've enumerated all possible rotations of the Assignment other, none of 
   * those have matched *this. So...
   */
  return false;
}

/*template<
  template<typename T = AssignmentColumn>
  class Symmetry
> std::set<
  Assignment<Symmetry>
> Assignment<Symmetry>::generateAllRotations() const {
  std::set<
    Assignment<Symmetr>
  > rotations = {*this};

  // Systematically explore all rotations
  // maximum element is the size of the rotation vector
  unsigned linkLimit = Symmetry<>::rotations.size();
  // initialize 
  std::vector<unsigned> chain = {0};
  unsigned depth = 0;
  while(chain.at(0) < linkLimit) {
    // perform instruction
    Assignment<Symmetry> generated = other;
    for(const auto& link: chain) {
      generated.applyRotation(
        Symmetry<>::rotations[link].first
      );
    }

    // is it something new?
    if(enumeratedAssignments.count(generated) == 0) {
      // is it the same as this?
      if(*this == generated) {
        return true;
      }
      // add it to the set
      enumeratedAssignments.insert(generated);
      // increase depth, add a link
      depth++;
      chain.emplace_back(0);
    } else {
      // if we are not at the maximum instruction
      if(chain.at(depth) < linkLimit - 1) {
        chain.at(depth)++;
      } else {
        // collapse the chain until we are at an incrementable position
        while(
          depth > 0
          && chain.at(depth) == linkLimit - 1
        ) {
          chain.pop_back();
          depth--;
        }

        chain.at(depth)++;
      }
    }
  }

} */

/* Operators */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
> bool Assignment<Symmetry>::operator < (
  const Assignment<Symmetry>& other
) const {
  if(
    this -> positionOccupations.size()
    < other.positionOccupations.size()
  ) {
    return true;
  } else if(
    this -> positionOccupations.size()
    > other.positionOccupations.size()
  ) {
    return false;
  }
  // positionOccupations.size() is equal

  // compare AssignmentColumns
  for(unsigned i = 0; i < this -> positionOccupations.size(); i++) {
    if(
      this -> positionOccupations[i]
      < other.positionOccupations[i]
    ) {
      return true;
    } else if(
      this -> positionOccupations[i]
      > other.positionOccupations[i]
    ) {
      return false;
    }
  }

  // all AssignmentColumns are equal
  return false;
}

template<
  template<typename T = AssignmentColumn>
  class Symmetry
> bool Assignment<Symmetry>::operator == (
  const Assignment<Symmetry>& other
) const {
  // compare characters
  for(unsigned i = 0; i < Symmetry<>::size; i++) {
    if(
      this -> positionOccupations[i].character != 
        other.positionOccupations[i].character
    ) {
      return false;
    }
  }

  // compare reduced groups
  if(!_reducedGroupsAreEqual(
    _reduceGroups(),
    other._reduceGroups()
  )) {
    return false;
  }
  
  return true;
}

/* Private members */
/*!
 * Compares two reduced group representations
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
> bool Assignment<Symmetry>::_reducedGroupsAreEqual(
  const std::vector<
    std::vector<unsigned>
  >& a,
  const std::vector<
    std::vector<unsigned>
  >& b
) const {
  // every element in a must be in b
  if(a.size() != b.size()) return false;
  
  for(const auto& aVectorElement : a) {
    if(!std::accumulate(
      b.begin(),
      b.end(),
      false,
      [&aVectorElement](
        const bool& carry,
        const std::vector<unsigned>& bVectorElement
      ) {
        if(carry) return carry;
        else return (
          carry
          || std::equal(
            aVectorElement.begin(),
            aVectorElement.end(),
            bVectorElement.begin()
          )
        );
      }
    )) {
      return false;
    }
  }

  return true;
}

/*!
 * Reduces the bool vector groups to pairs of connected indices.
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
> std::vector<
  std::vector<unsigned>
> Assignment<Symmetry>::_reduceGroups() const {
  std::vector<
    std::vector<unsigned>
  > groupReduction;

  for(unsigned i = 0; i < Symmetry<>::size; i++) {
    if(groupReduction.size() < positionOccupations[i].groups.size()) {
      groupReduction.resize(
        positionOccupations[i].groups.size()
      );
    }

    for(unsigned j = 0; j < positionOccupations[i].groups.size(); j++) {
      if(positionOccupations[i].groups[j]) {
        groupReduction[j].push_back(i);
      }
    }
  }

  // sort! this is required from _reducedGroupsAreEqual
  for(auto& reducedGroup : groupReduction) {
    std::sort(
      reducedGroup.begin(),
      reducedGroup.end()
    );
  }

  return groupReduction;
}

/*!
 * ostream operator for easier debugging
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
> std::ostream& operator << (
  std::ostream& os,
  const Assignment<Symmetry>& a
) {
  // make group.size + 1 stringstreams
  std::vector<
    std::stringstream
  > streams (
    1 + a.positionOccupations[0].groups.size()
  );

  for(const auto& column : a.positionOccupations ) {
    streams[0] << column.character << " ";
    for(unsigned i = 0; i < column.groups.size(); i++) {
      streams[i + 1] << (column.groups[i] ? "T" : "F") << " ";
    }
  }

  for(const auto& stream : streams) {
    os << stream.str() << std::endl;
  }

  return os;
}
