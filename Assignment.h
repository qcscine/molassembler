#ifndef LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H
#define LIB_UNIQUE_ASSIGNMENTS_ASSIGNMENT_H

#include <vector>
#include <algorithm>
#include <map>
#include <cassert>

#include "Util.h"

/* TODO
 * - update documentation
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

  /* Modifiers ––––––––––––––––––––––––––––*/
  //! Swap two "columns"
  void columnSwap(
    const unsigned& a,
    const unsigned& b
  );

  //! Transform this Assignment into its lowest permutation.
  void lowestPermutation() {
    /* laziest way to implement is to call nextPermutation until it returns 
     * false, at which point the data structure is reset to its lowest
     * permutation.
     */
    if(!isSortedAsc()) {
      while(nextPermutation()) continue;
    }
  }

  //! Modify the "columns" to the previous permutation
  bool nextPermutation();

  //! Modify the "columns" to the previous permutation
  bool previousPermutation();

  //! Reverse a span of "columns"
  void reverseColumns(
    const unsigned& from,
    const unsigned& to
  ) {
    unsigned a = from, b = to;
    while(a != b && a != --b) columnSwap(a++, b);
  }

  //! Rotate charactes according to template symmetry
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

  //! Rotate links according to template symmetry
  LinksSetType rotateLinks(
    const LinksSetType& links,
    const unsigned& rotationFunctionIndex
  );

  /*!
   * Applies a Symmetry rotation.
   */
  void applyRotation(const unsigned& index) {
    characters = rotateCharacters(characters, index);
    links = rotateLinks(links, index);
  }

  /* Information –––––––––––––––––––––––– */
  /*!
   * Compares two "columns". An important note about how this works:
   * E.g. if we have: chars {A, A}, links {[0, 1]}, then columnSmaller(0, 1) is
   * false. These columns are considered equal size in order to avoid confusion
   * in the permutation code where the instruction columnSwap(0, 1) would have
   * no net effect.
   */
  bool columnSmaller(
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
        Util::removeIndexFromSet(makeConnectedIndicesSet(a), b),
        Util::removeIndexFromSet(makeConnectedIndicesSet(b), a)
      ).value_or(
        false
      )
    );
  }


  /*!
   * Generates a set of all rotational equivalents of this Assignment as 
   * defined by its symmetry template parameter.
   */
  std::set<
    Assignment<Symmetry>
  > generateAllRotations() const;

  /*!
   * Gets a map of ligand symbol character to position in the permutational 
   * symmetry. 
   */
  std::map<
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

    for(unsigned i = 0; i < Symmetry::size; i++) {
      const char& columnChar = characters[i];
      if(returnMap.count(columnChar) == 0) {
        returnMap[columnChar] = {i};
      } else {
        returnMap[columnChar].push_back(i);
      }
    }

    return returnMap;
  }

  //! Returns whether the "columns" are sorted in ascending order
  bool isSortedAsc() const {
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

  /*!
   * Checks whether this Assignment is rotationally superimposable with
   * another.
   */
  bool isRotationallySuperimposable(
    const Assignment<Symmetry>& other
  ) const;

  /*!
   * Makes a set of a "column"'s connected indices.
   * If e.g. chars {A, A, A, A, A, A}, links {[0, 1], [1, 2]}, then
   * makeConnectedIndicesSet(1) = set {0, 2}.
   */
  std::set<unsigned> makeConnectedIndicesSet(const unsigned& index) const {
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


  //! Converts positionOccupations into a string for display
  std::string toString() const;

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
/* Modifiers ––––––––––––––––––––––– */
template<class Symmetry>
void Assignment<Symmetry>::columnSwap(
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

template<class Symmetry>
bool Assignment<Symmetry>::nextPermutation() {
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

template<class Symmetry>
bool Assignment<Symmetry>::previousPermutation() {

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

template<class Symmetry>
typename Assignment<Symmetry>::LinksSetType Assignment<Symmetry>::rotateLinks(
  const LinksSetType& links,
  const unsigned& rotationFunctionIndex
) {

  auto rotateIndex = [&rotationFunctionIndex](
    const unsigned& from
  ) -> unsigned {
    const auto& symVec = Symmetry::rotations.at(rotationFunctionIndex);

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

template<class Symmetry> 
std::set<
    Assignment<Symmetry>
> Assignment<Symmetry>::generateAllRotations() const {
  return _generateAllRotations(
    [](
      const Assignment<Symmetry>& a __attribute__ ((unused)),
      const Assignment<Symmetry>& b __attribute__ ((unused))
    ) -> bool {
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
        const Assignment<Symmetry>& a __attribute__ ((unused)),
        const Assignment<Symmetry>& b
      ) -> bool {
        return b == other;
      }
    ).second
  );
}

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

/* Operators */
template<class Symmetry> 
bool Assignment<Symmetry>::operator < (
  const Assignment<Symmetry>& other
) const {
  return Util::compareSmaller(
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
  unsigned linkLimit = Symmetry::rotations.size();

  // initialize 
  std::vector<unsigned> chain = {0};
  std::vector<
    Assignment<Symmetry>
  > chainStructures = {*this};
  unsigned depth = 0;

  // begin loop
  while(chain.front() < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    Assignment<Symmetry> generated = chainStructures.back();

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
