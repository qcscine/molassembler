#include "Symmetries.h"

#include "boost/hana.hpp"
#include "constexpr_magic/Containers.h"
#include "constexpr_magic/Set.h"

namespace Symmetry {

template<typename SymmetryClass, size_t... Indices>
constexpr std::array<unsigned, SymmetryClass::size> applyRotationImpl(
  const std::array<unsigned, SymmetryClass::size>& indices, 
  const unsigned& rotationFunctionIndex,
  std::index_sequence<Indices...>
) {
  return {{
    indices.at(
      SymmetryClass::rotations.at(rotationFunctionIndex).at(Indices)
    )...
  }};
}

template<typename SymmetryClass>
constexpr std::array<unsigned, SymmetryClass::size> applyRotation(
  const std::array<unsigned, SymmetryClass::size>& indices, 
  const unsigned& rotationFunctionIndex
) {
  return applyRotationImpl<SymmetryClass>(
    indices,
    rotationFunctionIndex,
    std::make_index_sequence<SymmetryClass::size>{}
  );
}

template<typename SymmetryClass>
constexpr ConstexprMagic::Vector getCoordinates(const unsigned& indexInSymmetry) {
  if(indexInSymmetry != replaceMe) {
    return SymmetryClass::coordinates.at(indexInSymmetry);
  }

  return {0, 0, 0};
}

constexpr double getTetrahedronVolume(
  const ConstexprMagic::Vector& i,
  const ConstexprMagic::Vector& j,
  const ConstexprMagic::Vector& k,
  const ConstexprMagic::Vector& l
) {
  return (i - l).dot(
    (j - l).cross(k - l)
  );
}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
constexpr double calculateAngleDistortion(
  const std::array<
    unsigned,
    ConstexprMagic::Math::max(
      SymmetryClassFrom::size,
      SymmetryClassTo::size
    )
  >& indexMapping
) {
  double angleDistortion = 0;

  for(unsigned i = 0; i < SymmetryClassFrom::size; ++i) {
    for(unsigned j = i + 1; j < SymmetryClassTo::size; ++j) {
      angleDistortion += ConstexprMagic::Math::abs(
        SymmetryClassFrom::angleFunction(i, j)
        - SymmetryClassTo::angleFunction(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angleDistortion;
}

template<size_t size>
unsigned propagateSymmetryPosition(
  const unsigned& symmetryPosition,
  const std::array<unsigned, size>& indexMapping
) {
  if(symmetryPosition != replaceMe) {
    return indexMapping.at(symmetryPosition);
  }

  return replaceMe;
}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
constexpr double calculateChiralDistortion(
  const std::array<
    unsigned,
    ConstexprMagic::Math::max(
      SymmetryClassFrom::size,
      SymmetryClassTo::size
    )
  >& indexMapping
) {
  double chiralDistortion = 0;
  
  for(const auto& tetrahedron : SymmetryClassFrom::tetrahedra) {
    chiralDistortion += ConstexprMagic::Math::abs(
      getTetrahedronVolume(
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(0)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(1)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(2)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(0), indexMapping)),
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(1), indexMapping)),
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(2), indexMapping)),
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(3), indexMapping))
      )
    );
  }

  return chiralDistortion;
}

template<size_t size>
constexpr std::array<unsigned, size> symPosMapping(
  const std::array<unsigned, size>& mapping
) {
  std::array<unsigned, size> symmetryPositions;

  for(unsigned i = 0; i < size; ++i) {
    symmetryPositions.at(
      mapping.at(i)
    ) = i;
  }

  return symmetryPositions;
}

template<size_t size>
struct DistortionInfo {
  std::array<unsigned, size> indexMapping;
  double totalDistortion;
  double chiralDistortion;

  constexpr DistortionInfo(
    const std::array<unsigned, size>& passIndexMapping,
    const double& passTotalDistortion,
    const double& passChiralDistortion
  ) : indexMapping(passIndexMapping),
      totalDistortion(passTotalDistortion),
      chiralDistortion(passChiralDistortion)
  {}
};

template<typename SymmetryClass>
using IndicesList = std::array<unsigned, SymmetryClass::size>;

enum class enabler_t {};

template<typename T>
using EnableIf = typename std::enable_if<T::value, enabler_t>::type;

template<typename SymmetryClass, size_t size>
constexpr std::enable_if_t<
  size == 0,
  std::pair<
    std::array<IndicesList<SymmetryClass>, 1>,
    std::array<unsigned, 1>
  >
> collapseChain(
  const std::array<IndicesList<SymmetryClass>, size>& chainStructures,
  const std::array<unsigned, size>& chain
) {
  return {};
}

template<typename SymmetryClass, size_t size>
constexpr auto collapseChain(
  const std::array<IndicesList<SymmetryClass>, size>& chainStructures,
  const std::array<unsigned, size>& chain
) {
  if(size == 1 || chain.back() == SymmetryClass::rotations.size()) {
    // Increment last item in chain via pop and push
    const auto lastItem = chain.back();
    const auto intermediateChain = ConstexprMagic::arrayPop(chain);

    return std::make_pair(
      chainStructures,
      ConstexprMagic::arrayPush(intermediateChain, lastItem + 1)
    );
  }

  return collapseChain(
    ConstexprMagic::arrayPop(chainStructures),
    ConstexprMagic::arrayPop(chain)
  );
}

template<
  typename SymmetryClass,
  class RotationsSetType,
  size_t chainLength
> constexpr auto generateAllRotationsImpl(
  const RotationsSetType& rotations,
  const std::array<IndicesList<SymmetryClass>, chainLength>& chainStructures,
  const std::array<unsigned, chainLength>& chain
) {
  if(chain.front() >= SymmetryClass::rotations.size()) {
    return rotations;
  }

  const auto newRotation = applyRotation<SymmetryClass>(
    chainStructures.back(),
    chain.back()
  );

  if(!rotations.contains(newRotation)) {
    return generateAllRotationsImpl<SymmetryClass>(
      rotations.insert(newRotation),
      ConstexprMagic::arrayPush(chainStructures, newRotation),
      ConstexprMagic::arrayPush(chain, 0u)
    );
  }

  if(chain.back() < SymmetryClass::rotations.size() - 1) {
    // increment last item via pop and push
    const auto lastItem = chain.back();
    auto intermediateChain = ConstexprMagic::arrayPop(chain);
    
    return generateAllRotationsImpl(
      rotations,
      chainStructures,
      ConstexprMagic::arrayPush(intermediateChain, lastItem + 1)
    );
  }

  const auto chainStructuresAndChainPair = collapseChain(
    chainStructures,
    chain
  );

  return generateAllRotationsImpl(
    rotations,
    chainStructuresAndChainPair.first,
    chainStructuresAndChainPair.second
  );
}

template<typename T, size_t size>
struct arrayComparator {
  static constexpr bool op(
    const std::array<T, size>& a,
    const std::array<T, size>& b
  ) {
    return ConstexprMagic::arraysEqual(a, b);
  }
};

template<typename SymmetryClass>
constexpr auto generateAllRotations(
  const std::array<unsigned, SymmetryClass::size>& indices
) {
  ConstexprMagic::Set<
    IndicesList<SymmetryClass>,
    arrayComparator<unsigned, SymmetryClass::size>
  > setInstance;

  return generateAllRotationsImpl<SymmetryClass>(
    setInstance,
    std::array<IndicesList<SymmetryClass>, 1> {{indices}},
    std::array<unsigned, 1> {{0}}
  );
}

constexpr auto maxAsymTetrRots = generateAllRotations<data::Tetrahedral>(
  std::array<unsigned, 4> {{0, 1, 2, 3}}
);

} // namespace Symmetry
