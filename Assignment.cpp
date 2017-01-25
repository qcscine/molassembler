#include "Assignment.h"
#include "Util.h"

#include <cassert>
#include <algorithm>

namespace UniqueAssignments {

/* Public members */
/*  Constructors */
Assignment::Assignment(
  const Symmetry::Name& passSymmetryName,
  const std::vector<char>& passCharacters
) : symmetryName(passSymmetryName),
    characters(passCharacters)
{
  assert(characters.size() == Symmetry::size(passSymmetryName));
}

Assignment::Assignment(
  const Symmetry::Name& passSymmetryName,
  const std::vector<char>& passCharacters,
  const LinksSetType& passLinks
) : symmetryName(passSymmetryName),
    characters(passCharacters),
    links(passLinks) 
{
  // make sure the number of characters matches the current symmetry
  assert(characters.size() == Symmetry::size(passSymmetryName));
}

/* Public members */
/* Modifiers ––––––––––––––––––––––– */
void Assignment::columnSwap(
  const unsigned& a,
  const unsigned& b
) {
  // exchange characters, adjust all LinksSetType accordingly
  std::swap(characters[a], characters[b]);

  auto determineNewIndex = [&a, &b](const unsigned& index) -> unsigned {
    if(index == a) return b;
    else if(index == b) return a;
    else return index;
  };

  LinksSetType newLinks;
  unsigned i, j;
  for(const auto& pair: links) {
    i = determineNewIndex(pair.first);
    j = determineNewIndex(pair.second);

    newLinks.emplace(
      std::min(i, j),
      std::max(i, j)
    );
  }

  // overwrite links
  links = std::move(newLinks);
}

bool Assignment::nextPermutation() {
  /* This is where it gets interesting. The elements we are permuting are the
   * characters ALONG with their indices in passLinks, and a permutation must
   * be a different operation than a rotation.
   *
   * The algorithm here is a copy of std::next_permutation, modified from
   * iterator-based to index-based and using the custom comparison and
   * swapping functions defined in-class.
   */
  unsigned i = characters.size() - 1, j, k;

  while(true) {
    j = i;

    if(
      i != 0 
      && columnSmaller(--i, j)
    ) {
      k = characters.size();

      while(
        k != 0 
        && !columnSmaller(i, --k)
      ) continue;

      columnSwap(i, k);
      reverseColumns(j, characters.size());
      return true;
    }

    if(i == 0) {
      reverseColumns(0, characters.size());
      return false;
    }
  }
}

bool Assignment::previousPermutation() {

  unsigned i = characters.size() - 1, j, k;

  while(true) {
    j = i;

    if(
      i != 0 
      && columnSmaller(j, --i)
    ) {
      k = characters.size();

      while(
        k != 0 
        && !columnSmaller(--k, i)
      ) continue;

      columnSwap(i, k);
      reverseColumns(j, characters.size());
      return true;
    }

    if(i == 0) {
      reverseColumns(0, characters.size());
      return false;
    }
  }

}

typename Assignment::LinksSetType Assignment::rotateLinks(
  const LinksSetType& links,
  const unsigned& rotationFunctionIndex
) {

  auto rotateIndex = [&rotationFunctionIndex, this](
    const unsigned& from
  ) -> unsigned {
    const auto& symVec = Symmetry::rotations(symmetryName).at(
      rotationFunctionIndex
    );

    return std::find(
      symVec.begin(),
      symVec.end(),
      from
    ) - symVec.begin();
  };

  LinksSetType retSet;

  for(const auto& pair : links) {
    retSet.emplace(
      Util::minMaxAdapt(
        std::function<
          std::pair<unsigned, unsigned>(
            unsigned,
            unsigned
          )
        >(std::make_pair<unsigned, unsigned>),
        rotateIndex(pair.first),
        rotateIndex(pair.second)
      )
    );
  }

  return retSet;
}

std::set<Assignment> Assignment::generateAllRotations() const {
  return _generateAllRotations(
    [](
      const Assignment& a __attribute__ ((unused)),
      const Assignment& b __attribute__ ((unused))
    ) -> bool {
      return false;
    }
  ).first;
}

bool Assignment::isRotationallySuperimposable(
  const Assignment& other
) const {
  return (
    *this == other
    || _generateAllRotations(
      [&other](
        const Assignment& a __attribute__ ((unused)),
        const Assignment& b
      ) -> bool {
        return b == other;
      }
    ).second
  );
}

//! Converts positionOccupations into a string
std::string Assignment::toString() const {

  std::stringstream out;

  out << "chars {";
  for(unsigned i = 0; i < characters.size(); i++) {
    out << characters[i];
    if(i != characters.size() - 1) {
      out << ", ";
    }
  }
  out << "}, links {";

  unsigned pairs = links.size();
  for(const auto& pair: links) {
    out << "[" << pair.first << ", " << pair.second << "]";
    if(--pairs != 0) out << ", ";
  }

  out << "}";

  return out.str();
}

/* Operators */
bool Assignment::operator < (
  const Assignment& other
) const {
  return Util::compareSmaller(
    this -> characters,
    other.characters
  ).value_or(
    this -> links < other.links // size and then lexicographical comparison
  );
}

bool Assignment::operator > (
  const Assignment& other
) const {
  return other < *this;
}

bool Assignment::operator == (
  const Assignment& other
) const {
  return (
    this -> characters == other.characters
    && this -> links == other.links
  );
}

bool Assignment::operator != (
  const Assignment& other
) const {
  return !(*this == other);
}

/* Private members */
std::pair<
  std::set<Assignment>,
  bool
> Assignment::_generateAllRotations(
  std::function<
    bool(const Assignment&, const Assignment&)
  > interruptCallbackOnNewAssignment
) const {

  // add the initial structure to a set of Assignments
  std::set<Assignment> enumeratedAssignments = {*this};   

  /* Systematically explore all rotations:
   * The chain exists to keep track of which rotations starting from the base
   * structure we have already explored, and which to perform next on which 
   * intermediate structure. Say we have 4 rotations defined in the current
   * symmetry, then linkLimit is 4.
   *
   * Chain initially is {0} and chainStructures is {copy of *this}. 
   * So in the loop we copy out the last element of chainStructures and perform
   * the rotation 0 on it. 
   *
   * If it's something new, we add it to a set of generated structures, push it
   * onto the chain of structures, and add an instruction to perform the
   * rotation 0 on the new structure.
   *
   * If it isn't something new, we increment the instruction at the end of 
   * the chain. If that happens to be linkLimit, then we pop off the last 
   * elements of both chains and increment the instruction at the new end 
   * position.
   *
   * This leads to a tree traversal that prunes the tree whenever an
   * instruction produces a structure that has already been seen.
   */
  // maximum element is the size of the rotation vector
  unsigned linkLimit = Symmetry::rotations(symmetryName).size();

  // initialize 
  std::vector<unsigned> chain = {0};
  std::vector<Assignment> chainStructures = {*this};
  unsigned depth = 0;

  // begin loop
  while(chain.front() < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    Assignment generated = chainStructures.back();

    // apply the rotation referenced by the last link in chain
    generated.applyRotation(
      chain.back()
    );

    // is it something new?
    if(enumeratedAssignments.count(generated) == 0) {
      /* give a chance to interrupt if a condition for *this and the newly
       *  generated structures is fulfilled
       */
      if(interruptCallbackOnNewAssignment(*this, generated)) {
        return make_pair(enumeratedAssignments, true);
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

  return make_pair(enumeratedAssignments, false);
}

void Assignment::lowestPermutation() {
  /* laziest way to implement is to call nextPermutation until it returns 
   * false, at which point the data structure is reset to its lowest
   * permutation.
   */
  if(!isSortedAsc()) {
    while(nextPermutation()) continue;
  }
}

void Assignment::reverseColumns(
  const unsigned& from,
  const unsigned& to
) {
  unsigned a = from, b = to;
  while(a != b && a != --b) columnSwap(a++, b);
}

std::vector<char> Assignment::rotateCharacters(
  const std::vector<char>& characters,
  const unsigned& rotationFunctionIndex
) {
  std::vector<char> retv;

  for(
    const auto& index : 
    Symmetry::rotations(symmetryName).at(rotationFunctionIndex)
  ) {
    retv.push_back(
      characters.at(index)
    );
  }
  
  return retv;
}

void Assignment::applyRotation(const unsigned& index) {
  characters = rotateCharacters(characters, index);
  links = rotateLinks(links, index);
}

bool Assignment::columnSmaller(
  const unsigned& a,
  const unsigned& b
) const {
  if(links.size() == 0) return characters[a] < characters[b];

  return Util::compareSmaller(
    characters[a],
    characters[b]
  ).value_or(
    Util::compareSmaller(
      // see if this does the trick
      Util::removeElementFromSet(makeConnectedIndicesSet(a), b),
      Util::removeElementFromSet(makeConnectedIndicesSet(b), a)
    ).value_or(
      false
    )
  );
}

std::map<
  char,
  std::vector<unsigned>
> Assignment::getCharMap() const {
  std::map<
    char,
    std::vector<
      unsigned
    >
  > returnMap;

  for(unsigned i = 0; i < Symmetry::size(symmetryName); i++) {
    const char& columnChar = characters[i];
    if(returnMap.count(columnChar) == 0) {
      returnMap[columnChar] = {i};
    } else {
      returnMap[columnChar].push_back(i);
    }
  }

  return returnMap;
}

bool Assignment::isSortedAsc() const {
  bool isSorted = true;
  for(unsigned i = 0; i < characters.size() - 1; i++) {
    /* It is a mistake with to test !columnSmaller(i, i + 1) as equal columns
     * can arise, e.g. chars {A, A}, links {[0, 1]}, in which 
     * columnSmaller(0, 1) = false. Using abovementioned test also places the
     * restriction of monotonic increase upon the columns. It is better to
     * test columnSmaller(i + 1, i).
     */
    if(columnSmaller(i + 1, i)) { 
      isSorted = false;
      break;
    }
  }
  
  return isSorted;
}

std::set<unsigned> Assignment::makeConnectedIndicesSet(
  const unsigned& index
) const {
    std::set<unsigned> connectedIndices;

    for(const auto& pair: links) {
      if(pair.first == index) {
        connectedIndices.insert(pair.second);
      }
      if(pair.second == index) {
        connectedIndices.insert(pair.first);
      }
    }

    return connectedIndices;
}

/*!
 * ostream operator for easier debugging
 */
std::ostream& operator << (
  std::ostream& os,
  const Assignment& a
) {
  os << a.toString();

  return os;
}

} // eo namespace
