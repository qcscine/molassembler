#ifndef LIB_SYMMETRY_ROTATIONS_H
#define LIB_SYMMETRY_ROTATIONS_H

#include <cassert>
#include <functional>
#include <vector>
#include <set>

namespace SymmetryRotations {

  using LinksSetType = std::set<
    std::set<unsigned>
  >;

template<typename T, typename Symmetry>
std::function<
  std::vector<T>(const std::vector<T>&)
> makeVectorRotationFunction(const unsigned& rotationFunctionIndex) {
  return [&rotationFunctionIndex](
    const std::vector<T>& source
  ) -> std::vector<T> {
    assert(source.size() == Symmetry::size());

    std::vector<T> retv;

    for(const auto& index: Symmetry::rotations.at(rotationFunctionIndex)) {
      retv.push_back(
        source.at(index)
      );
    }
    
    return retv;
  };
}

template<typename Symmetry>
std::vector<char> rotateCharacters(
  const std::vector<char>& characters,
  const unsigned& rotationFunctionIndex
) {
  assert(characters.size() == Symmetry::size());

  std::vector<char> retv;

  for(const auto& index: Symmetry::rotations.at(rotationFunctionIndex)) {
    retv.push_back(
      characters.at(index)
    );
  }
  
  return retv;
}

template<typename Symmetry>
std::function<LinksSetType(const LinksSetType&)> makeLinksRotationFunction(
  const unsigned& rotationFunctionIndex
) {
  return [&rotationFunctionIndex](const LinksSetType& links) -> LinksSetType {
    LinksSetType retSets;

    for(const auto& set : links) {
      std::set<unsigned> transformedCopy;
      for(const auto& index: set) {
        transformedCopy.insert(
          Symmetry::rotations.at(
            rotationFunctionIndex
          ).at(index)
        );
      }

      retSets.insert(transformedCopy);
    }

    return retSets;
  };
}

template<typename Symmetry>
LinksSetType rotateLinks(
  const LinksSetType& links,
  const unsigned& rotationFunctionIndex
) {
  LinksSetType retSets;

  for(const auto& set : links) {
    std::set<unsigned> transformedCopy;
    for(const auto& index: set) {
      transformedCopy.insert(
        Symmetry::rotations.at(
          rotationFunctionIndex
        ).at(index)
      );
    }

    retSets.insert(transformedCopy);
  }

  return retSets;
}

} // eo namespace

#endif


