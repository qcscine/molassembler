#include "Stereopermutation.h"

#include "boost/functional/hash.hpp"
#include "chemical_symmetries/Symmetries.h"

#include "Util.h"

#include <algorithm>
#include <cassert>

namespace stereopermutation {

/* Public members */
/*  Constructors */
Stereopermutation::Stereopermutation(
  const Symmetry::Name passSymmetryName,
  std::vector<char> passCharacters
) : characters(std::move(passCharacters))
{
  assert(characters.size() == Symmetry::size(passSymmetryName));
}

Stereopermutation::Stereopermutation(
  const Symmetry::Name passSymmetryName,
  std::vector<char> passCharacters,
  LinksSetType passLinks
) : characters(std::move(passCharacters)),
    links(std::move(passLinks))
{
  // make sure the number of characters matches the current symmetry
  assert(characters.size() == Symmetry::size(passSymmetryName));

  /* make sure all links are properly self-referential, i.e. only contain valid
   * indices to the characters
   */
#ifndef NDEBUG
  for(const auto& linkPair : links) {
    assert(
      linkPair.first < characters.size()
      && linkPair.second < characters.size()
    );
  }
#endif
}

/* Public members */
/* Modifiers ––––––––––––––––––––––– */
void Stereopermutation::columnSwap(
  const unsigned a,
  const unsigned b
) {
  assert(a < characters.size() && b < characters.size());

  // exchange characters, adjust all LinksSetType accordingly
  std::swap(characters[a], characters[b]);

  auto determineNewIndex = [&a, &b](const unsigned index) -> unsigned {
    if(index == a) {
      return b;
    }

    if(index == b) {
      return a;
    }

    return index;
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

bool Stereopermutation::nextPermutation() {
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
      ) { }

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

bool Stereopermutation::previousPermutation() {

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
      ) { }

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

typename Stereopermutation::LinksSetType Stereopermutation::rotateLinks(
  const LinksSetType& links,
  const std::vector<unsigned>& rotationIndices
) const {
  auto rotateIndex = [&rotationIndices](const unsigned from) -> unsigned {
    return std::find(
      rotationIndices.begin(),
      rotationIndices.end(),
      from
    ) - rotationIndices.begin();
  };

  LinksSetType rotatedSet;

  for(const auto& pair : links) {
    rotatedSet.emplace(
      Util::sortBinaryArgs(
        std::make_pair<unsigned, unsigned>,
        rotateIndex(pair.first),
        rotateIndex(pair.second)
      )
    );
  }

  return rotatedSet;
}

std::set<Stereopermutation> Stereopermutation::generateAllRotations(
  const Symmetry::Name& symmetryName
) const {
  return _generateAllRotations(
    [](
      const Stereopermutation& /* a */,
      const Stereopermutation& /* b */
    ) -> bool {
      // Never interrupt
      return false;
    },
    symmetryName
  ).first;
}

bool Stereopermutation::isRotationallySuperimposable(
  const Stereopermutation& other,
  const Symmetry::Name& symmetryName
) const {
  return (
    *this == other
    || _generateAllRotations(
      [&other](
        const Stereopermutation& /* a */,
        const Stereopermutation& b
      ) -> bool {
        return b == other;
      },
      symmetryName
    ).second
  );
}

//! Converts positionOccupations into a string
std::string Stereopermutation::toString() const {

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
    if(--pairs != 0) {
      out << ", ";
    }
  }

  out << "}";

  return out.str();
}

/* Operators */
bool Stereopermutation::operator < (
  const Stereopermutation& other
) const {
  return Util::compareSmaller(
    this -> characters,
    other.characters
  ).value_or(
    this -> links < other.links // size and then lexicographical comparison
  );
}

bool Stereopermutation::operator > (
  const Stereopermutation& other
) const {
  return other < *this;
}

bool Stereopermutation::operator == (
  const Stereopermutation& other
) const {
  return (
    this -> characters == other.characters
    && this -> links == other.links
  );
}

bool Stereopermutation::operator != (
  const Stereopermutation& other
) const {
  return !(*this == other);
}

/* Private members */
std::pair<
  std::set<Stereopermutation>,
  bool
> Stereopermutation::_generateAllRotations(
  const std::function<
    bool(const Stereopermutation&, const Stereopermutation&)
  >& interruptCallbackOnNewStereopermutation,
  const Symmetry::Name& symmetryName
) const {

  // add the initial structure to a set of Stereopermutations
  std::set<Stereopermutation> enumeratedStereopermutations = {*this};

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
  std::vector<Stereopermutation> chainStructures = {*this};
  unsigned depth = 0;

  // begin loop
  while(chain.front() < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    Stereopermutation generated = chainStructures.back();

    // apply the rotation referenced by the last link in chain
    generated.applyRotation(
      Symmetry::rotations(symmetryName).at(
        chain.back()
      )
    );

    // is it something new?
    if(enumeratedStereopermutations.count(generated) == 0) {
      /* give a chance to interrupt if a condition for *this and the newly
       *  generated structures is fulfilled
       */
      if(interruptCallbackOnNewStereopermutation(*this, generated)) {
        return make_pair(enumeratedStereopermutations, true);
      }

      // add it to the set
      enumeratedStereopermutations.insert(generated);

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

  return make_pair(enumeratedStereopermutations, false);
}

void Stereopermutation::lowestPermutation() {
  /* laziest way to implement is to call nextPermutation until it returns
   * false, at which point the data structure is reset to its lowest
   * permutation.
   */
  if(!isSortedAsc()) {
    while(nextPermutation()) { }
  }
}

void Stereopermutation::reverseColumns(
  const unsigned from,
  const unsigned to
) {
  unsigned a = from, b = to;
  while(a != b && a != --b) {
    columnSwap(a++, b);
  }
}

std::vector<char> Stereopermutation::rotateCharacters(
  const std::vector<char>& characters,
  const std::vector<unsigned>& rotationIndices
) const {
  std::vector<char> retv;
  retv.reserve(characters.size());

  for(const auto& index : rotationIndices) {
    retv.push_back(
      characters.at(index)
    );
  }

  return retv;
}

void Stereopermutation::applyRotation(const std::vector<unsigned>& rotationIndices) {
  characters = rotateCharacters(characters, rotationIndices);
  links = rotateLinks(links, rotationIndices);
}

bool Stereopermutation::columnSmaller(
  const unsigned a,
  const unsigned b
) const {
  if(links.empty()) {
    return characters[a] < characters[b];
  }

  return Util::compareSmaller(
    characters[a],
    characters[b]
  ).value_or(
    Util::compareSmaller(
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
> Stereopermutation::getCharMap() const {
  std::map<
    char,
    std::vector<unsigned>
  > returnMap;

  for(unsigned i = 0; i < characters.size(); i++) {
    const char& columnChar = characters[i];
    if(returnMap.count(columnChar) == 0) {
      returnMap[columnChar] = {i};
    } else {
      returnMap[columnChar].push_back(i);
    }
  }

  return returnMap;
}

bool Stereopermutation::isSortedAsc() const {
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

std::set<unsigned> Stereopermutation::makeConnectedIndicesSet(
  const unsigned index
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

std::size_t hash_value(const Stereopermutation& assignment) {
  std::size_t seed = 0;

  boost::hash_combine(
    seed,
    boost::hash_range(assignment.characters.begin(), assignment.characters.end())
  );

  boost::hash_combine(
    seed,
    boost::hash_range(assignment.links.begin(), assignment.links.end())
  );

  return seed;
}

/*!
 * ostream operator for easier debugging
 */
std::ostream& operator << (
  std::ostream& os,
  const Stereopermutation& a
) {
  os << a.toString();

  return os;
}

} // namespace stereopermutation
