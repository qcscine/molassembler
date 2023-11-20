/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutation/Stereopermutation.h"

#include "boost/functional/hash.hpp"
#include "boost/optional.hpp"
#include "Molassembler/Shapes/Data.h"

#include "Molassembler/Temple/Adaptors/Iota.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Permutations.h"

#include <algorithm>
#include <cassert>

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {

Stereopermutation::Stereopermutation(
  Occupation passOccupation,
  OrderedLinks passLinks
) : occupation(std::move(passOccupation)),
    links(std::move(passLinks))
{
  /* make sure all links are properly self-referential, i.e. only contain valid
   * indices to the characters
   */
  const unsigned S = occupation.size();
  if(
    Temple::any_of(
      links,
      [S](const auto& link) { return link.first >= S || link.second >= S || link.first == link.second;}
    )
  ) {
    throw std::out_of_range("Links are invalid");
  }

  // Fix inverted links
  for(Link& link : links) {
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
  using ShapePermutation = Temple::StrongIndexPermutation<Shapes::Vertex, Shapes::Vertex>;
  auto p = ShapePermutation::from(permutation).inverse();

  return Temple::sorted(
    Temple::map(
      links,
      [&](const Link& link) {
        Link mappedLink {
          p(link.first),
          p(link.second)
        };

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

  out << "occupation '";
  for(const auto& pair : occupation) {
    out << static_cast<char>('A' + pair.second);
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

Stereopermutation::Occupation Stereopermutation::permuteOccupation(
  const Occupation& occupation,
  const Shapes::Permutation& permutation
) {
  using ShapePermutation = Temple::StrongIndexPermutation<Shapes::Vertex, Shapes::Vertex>;
  return ShapePermutation::from(permutation).compose(occupation);
}

Stereopermutation::Occupation Stereopermutation::occupationFromChars(const std::string& chars) {
  return Occupation {
    Temple::map(chars, [](const char c) -> unsigned {
      const int result = c - 'A';
      if(result < 0) {
        throw std::runtime_error("Cannot generate occupation from chars less than 'A'");
      }
      return static_cast<unsigned>(result);
    })
  };
}

Stereopermutation Stereopermutation::applyPermutation(const Shapes::Permutation& permutation) const {
  return {
    permuteOccupation(occupation, permutation),
    permuteLinks(links, permutation)
  };
}

std::size_t hash_value(const Stereopermutation& assignment) {
  std::size_t seed = 0;

  boost::hash_combine(
    seed,
    boost::hash_range(
      std::begin(assignment.occupation.permutation.sigma),
      std::end(assignment.occupation.permutation.sigma)
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
