/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutation/Stereopermutation.h"

#include "boost/functional/hash.hpp"
#include "boost/optional.hpp"
#include "Molassembler/Shapes/Data.h"

#include "Molassembler/Temple/Adaptors/Iota.h"
#include "Molassembler/Temple/Functional.h"

#include <algorithm>
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {

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
    Temple::any_of(
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
    Temple::sort(links);
  }
}

typename Stereopermutation::OrderedLinks Stereopermutation::permuteLinks(
  const OrderedLinks& links,
  const Shapes::Permutation& permutation
) {
  const auto rotateIndex = [&permutation](const Shapes::Vertex from) -> Shapes::Vertex {
    return Shapes::Vertex(
      Temple::find(permutation, from) - std::begin(permutation)
    );
  };

  return Temple::sorted(
    Temple::map(
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

  out << "chars '";
  for(const char character : characters) {
    out << character;
  }
  out << "', links [";

  unsigned pairs = links.size();
  for(const auto& pair: links) {
    out << "(" << pair.first << ", " << pair.second << ")";
    if(--pairs != 0) {
      out << ", ";
    }
  }

  out << "]";

  return out.str();
}

Stereopermutation::CharacterOccupation Stereopermutation::permuteCharacters(
  const CharacterOccupation& characters,
  const Shapes::Permutation& permutation
) {
  return Temple::map(permutation, Temple::Functor::at(characters));
}

Stereopermutation Stereopermutation::applyPermutation(const Shapes::Permutation& permutation) const {
  return {
    permuteCharacters(characters, permutation),
    permuteLinks(links, permutation)
  };
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

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine
