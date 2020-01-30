/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "stereopermutation/Stereopermutation.h"

#include "boost/functional/hash.hpp"
#include "boost/optional.hpp"
#include "shapes/Data.h"

#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"

#include <algorithm>
#include <cassert>

namespace Scine {
namespace stereopermutation {

Stereopermutation::Stereopermutation(
  CharacterOccupation passCharacters,
  OrderedLinks passLinks
) : characters(std::move(passCharacters)),
    links(std::move(passLinks))
{
  /* make sure all links are properly self-referential, i.e. only contain valid
   * indices to the characters
   */
  const unsigned S = characters.size();
  if(
    temple::any_of(
      links,
      [S](const auto& link) { return link.first >= S || link.second >= S || link.first == link.second;}
    )
  ) {
    throw std::out_of_range("Links are invalid");
  }

  // Fix inverted links
  for(auto& link : links) {
    if(link.second < link.first) {
      std::swap(link.first, link.second);
    }
  }

  // Fix unsorted links
  if(!std::is_sorted(std::begin(links), std::end(links))) {
    temple::inplace::sort(links);
  }
}

typename Stereopermutation::OrderedLinks Stereopermutation::permuteLinks(
  const OrderedLinks& links,
  const shapes::Permutation& permutation
) {
  auto rotateIndex = [&permutation](const shapes::Vertex from) -> shapes::Vertex {
    const unsigned indexOf = std::find(
      std::begin(permutation),
      std::end(permutation),
      from
    ) - std::begin(permutation);
    return shapes::Vertex {indexOf};
  };

  return temple::sort(
    temple::map(
      links,
      [&](const auto& link) {
        auto mappedLink = std::make_pair(
          rotateIndex(link.first),
          rotateIndex(link.second)
        );

        if(mappedLink.second < mappedLink.first) {
          std::swap(mappedLink.first, mappedLink.second);
        }

        return mappedLink;
      }
    )
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

Stereopermutation::CharacterOccupation Stereopermutation::permuteCharacters(
  const CharacterOccupation& characters,
  const shapes::Permutation& permutation
) {
  return temple::map(permutation, temple::functor::at(characters));
}

Stereopermutation Stereopermutation::applyPermutation(const shapes::Permutation& permutation) const {
  return {
    permuteCharacters(characters, permutation),
    permuteLinks(links, permutation)
  };
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
    const char columnChar = characters[i];
    // C++17 insert_or_update
    if(returnMap.count(columnChar) == 0) {
      returnMap[columnChar] = {i};
    } else {
      returnMap[columnChar].push_back(i);
    }
  }

  return returnMap;
}

std::size_t hash_value(const Stereopermutation& assignment) {
  std::size_t seed = 0;

  boost::hash_combine(
    seed,
    boost::hash_range(
      std::begin(assignment.characters),
      std::end(assignment.characters)
    )
  );

  boost::hash_combine(
    seed,
    boost::hash_range(
      std::begin(assignment.links),
      std::end(assignment.links)
    )
  );

  return seed;
}

} // namespace stereopermutation
} // namespace Scine
