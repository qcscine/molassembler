#ifndef LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H
#define LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H

#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <set>
#include <map>
#include <cassert>

#include "SymmetryRotations.h"

#include "boost/optional.hpp"

template<typename Comparable>
boost::optional<bool> compareSmaller(
  const Comparable& a,
  const Comparable& b
) {
  if(a < b) return true;
  else if(b < a) return false;
  else return {};
}

/* TODO
 * - update documentation
 * - experiment with an alternative implementation where links are two-element
 *   sets to see performance impact
 * - reimplement getCharMap
 */

/* NOTES
 * - One path to implementation success of implementing permutations may be to
 *   implement iterators on Assignment that will then allow 
 *   std::next_permutation to run on Assignment! Positions, length, and when 
 *   characters are "different" must be clearly defined so as to mimic a 
 *   container 
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
template<class Symmetry>
struct Assignment {
private:
/* Private member functions */

public:
/* Typedefs */
  using LinksSetType = std::set<
    std::pair<unsigned, unsigned>
  >;

/* Public members */
  std::vector<char> characters;
  LinksSetType links;

  /* Constructors */
  Assignment() = delete;
  Assignment(
    const std::vector<char>& passCharacters
  );
  Assignment(
    const std::vector<char>& passCharacters,
    const LinksSetType& passLinks
  );

  void columnSwap(
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

  std::set<unsigned> makeConnectedIndicesSet(const unsigned& index) {
      std::set<unsigned> connectedIndices;

      for(const auto& pair: links) {
        if(pair.first == index) {
          connectedIndices.insert(pair.second);
        }
        if(pair.second == index) {
          connectedIndices.insert(pair.first);
        }
      }

      return std::move(connectedIndices);
  }

  bool columnSmaller(
    const unsigned& a,
    const unsigned& b
  ) {
    if(links.size() == 0) return characters[a] < characters[b];

    return compareSmaller(
      characters[a],
      characters[b]
    ).value_or(
      compareSmaller(
        makeConnectedIndicesSet(a),
        makeConnectedIndicesSet(b)
      ).value_or(
        false
      )
    );
  }

  void reverseColumns(
    const unsigned& from,
    const unsigned& to
  ) {
    unsigned a = from, b = to;
    while(a != b && a != --b) columnSwap(a++, b);
  }

  /* Modification */
  void lowestPermutation() {
    /* laziest way to implement is to call nextPermutation until it returns 
     * false, at which point the data structure is reset to lowest permutation.
     */
    while(nextPermutation()) continue;
  }

  /*!
   * Generates the next permutation of positionOccupations in which the 
   * ligand connections are ordered.
   */
  bool nextPermutation() {
    /* 
     * This is where it gets interesting. The elements we are permuting are the
     * characters ALONG with their indices in passLinks, and a permutation 
     * must be a different operation than a rotation.
     */
    unsigned i = characters.size() - 1, j, k;

    while(true) {
      j = i;

      if(
        i != 0 
        && columnSmaller(
          --i,
          j
        )
      ) {
        k = characters.size();

        while(
          k != 0 
          && !columnSmaller(
            i,
            --k
          )
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

  std::vector<char> rotateCharacters(
    const std::vector<char>& characters,
    const unsigned& rotationFunctionIndex
  ) {
    std::vector<char> retv;

    for(const auto& index: Symmetry::rotations.at(rotationFunctionIndex)) {
      retv.push_back(
        characters.at(index)
      );
    }
    
    return retv;
  }

  LinksSetType rotateLinks(
    const LinksSetType& links,
    const unsigned& rotationFunctionIndex
  ) {
    auto rotateIndex = [&rotationFunctionIndex](
      const unsigned& from
    ) -> unsigned {
      return Symmetry::rotations.at(
        rotationFunctionIndex
      ).at(from);
    };

    LinksSetType retSet;

    for(const auto& pair : links) {
      retSet.emplace(
        std::make_pair(
          rotateIndex(pair.first),
          rotateIndex(pair.second)
        )
      );
    }

    return retSet;
  }

  /*!
   * Applies a Symmetry rotation.
   */
  void applyRotation(const unsigned& index) {
    characters = rotateCharacters(characters, index);
    links = rotateLinks(links, index);
  }

  /* Information */
  /*!
   * Gets a map of ligand symbol character to position in the permutational 
   * symmetry. 
   *
   * ?? How to reimplement?
   */
  /*std::map<
    char,
    std::vector<
      unsigned
    >
  > getCharMap() const {
    std::map<
      char,
      std::vector<
        unsigned
      >
    > returnMap;

    for(unsigned i = 0; i < positionOccupations.size(); i++) {
      const char& columnChar = positionOccupations[i].character;
      if(returnMap.count(columnChar) == 0) {
        returnMap[columnChar] = {i};
      } else {
        returnMap[columnChar].push_back(i);
      }
    }

    return returnMap;
  }*/

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

#ifndef UNUSED
#define UNUSED(x) (void)(x)
#endif

/* Public members */
/*  Constructors */
/*!
 * Constructs an Assignment from a list of ligand characters.
 * \tparam Symmetry A SymmetryInformation derived class template.
 * \param characters A vector of chars signifying abstract ligands.
 */
template<class Symmetry>
Assignment<Symmetry>::Assignment(
  const std::vector<char>& passCharacters
) : characters(passCharacters) {
  assert(characters.size() == Symmetry::size);
}


/*!
 * Construct an Assignment from a list of ligand characters and a list of 
 * bonded indices referencing the ligand characters.
 * \tparam Symmetry A SymmetryInformation derived class template.
 * \param characters A vector of chars signifying abstract ligands.
 * \param pairedIndices A vector of pairs. Describes which ligand characters 
 *  are bonded to one another.
 */
template<class Symmetry>
Assignment<Symmetry>::Assignment(
  const std::vector<char>& passCharacters,
  const LinksSetType& passLinks
) : 
  characters(passCharacters),
  links(passLinks) {
  // make sure the number of characters matches the current symmetry
  assert(characters.size() == Symmetry::size);
}

/* Public members */
template<class Symmetry> 
std::set<
    Assignment<Symmetry>
> Assignment<Symmetry>::generateAllRotations() const {
  return _generateAllRotations(
    [](const Assignment<Symmetry>& a, const Assignment<Symmetry>& b) -> bool {
      UNUSED(a);
      UNUSED(b);
      return false;
    }
  ).first;
}

template<class Symmetry>
bool Assignment<Symmetry>::isRotationallySuperimposable(
  const Assignment<Symmetry>& other
) const {
  return (
    *this == other
    || _generateAllRotations(
      [&other](
        const Assignment<Symmetry>& a,
        const Assignment<Symmetry>& b
      ) -> bool {
        UNUSED(a);
        return b == other;
      }
    ).second
  );
}

template<class Symmetry> 
std::pair<
  std::set<
    Assignment<Symmetry>
  >,
  bool
> Assignment<Symmetry>::_generateAllRotations(
  std::function<
    bool(const Assignment<Symmetry>&, const Assignment<Symmetry>&)
  > interruptCallbackOnNewAssignment
) const {

  // add the initial structure to a set of Assignments
  std::set<Assignment> enumeratedAssignments = {*this};   

  // Systematically explore all rotations
  // maximum element is the size of the rotation vector
  unsigned linkLimit = Symmetry::rotations.size();

  // initialize 
  std::vector<unsigned> chain = {0};
  std::vector<
    Assignment<Symmetry>
  > chainStructures = {*this};
  unsigned depth = 0;

  // begin loop
  while(chain.at(0) < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    Assignment<Symmetry> generated = chainStructures.at(
      chainStructures.size() - 1
    );
    // apply the rotation referenced by the last link in chain
    generated.applyRotation(
      chain.at(
        chain.size() - 1
      )
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

/* Operators */
template<class Symmetry> 
bool Assignment<Symmetry>::operator < (
  const Assignment<Symmetry>& other
) const {
  return compareSmaller(
    this -> characters,
    other.characters
  ).value_or(
    this -> links < other.links // size and then lexicographical comparison
  );
}

template<class Symmetry>
bool Assignment<Symmetry>::operator > (
  const Assignment<Symmetry>& other
) const {
  return other < *this;
}

template<class Symmetry> 
bool Assignment<Symmetry>::operator == (
  const Assignment<Symmetry>& other
) const {
  return (
    this -> characters == other.characters
    && this -> links == other.links
  );
}

template<class Symmetry>
bool Assignment<Symmetry>::operator != (
  const Assignment<Symmetry>& other
) const {
  return !(*this == other);
}

/* Private members */
//! Converts positionOccupations into a string
template<class Symmetry>
std::string Assignment<Symmetry>::toString() const {

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

/*!
 * ostream operator for easier debugging
 */
template<class Symmetry> 
std::ostream& operator << (
  std::ostream& os,
  const Assignment<Symmetry>& a
) {
  os << a.toString();

  return os;
}

} // eo namespace

#endif
