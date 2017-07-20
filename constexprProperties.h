#include "Symmetries.h"

#include "boost/hana.hpp"
#include "constexpr_magic/Containers.h"
#include "constexpr_magic/Set.h"
#include "constexpr_magic/Math.h"

namespace Symmetry {

// When C++17 rolls around, replace this with std::array!
template<typename T, size_t size> 
using ArrayType = ConstexprMagic::Array<T, size>;

template<size_t size, size_t... Indices>
constexpr ArrayType<unsigned, size> iotaImpl(
  std::index_sequence<Indices...>
) {
  return {
    Indices...
  };
}

template<size_t size>
constexpr ArrayType<unsigned, size> iota() {
  return iotaImpl<size>(
    std::make_index_sequence<size>{}
  );
}

template<typename SymmetryClass, size_t... Indices>
constexpr ArrayType<unsigned, SymmetryClass::size> applyRotationImpl(
  const ArrayType<unsigned, SymmetryClass::size>& indices, 
  const unsigned& rotationFunctionIndex,
  std::index_sequence<Indices...>
) {
  return {
    indices.at(
      SymmetryClass::rotations.at(rotationFunctionIndex).at(Indices)
    )...
  };
}

template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::size> applyRotation(
  const ArrayType<unsigned, SymmetryClass::size>& indices, 
  const unsigned& rotationFunctionIndex
) {
  return applyRotationImpl<SymmetryClass>(
    indices,
    rotationFunctionIndex,
    std::make_index_sequence<SymmetryClass::size>{}
  );
}

template<typename SymmetryClass>
constexpr unsigned rotationMultiplicityImpl(
  const unsigned& rotationFunctionIndex,
  const ArrayType<unsigned, SymmetryClass::size>& runningIndices,
  const unsigned& count
) {
  if(
    ConstexprMagic::arraysEqual(
      runningIndices,
      iota<SymmetryClass::size>()
    )
  ) {
    return count;
  } else {
    return rotationMultiplicityImpl<SymmetryClass>(
      rotationFunctionIndex,
      applyRotation<SymmetryClass>(
        runningIndices,
        rotationFunctionIndex
      ),
      count + 1
    );
  }
}

template<typename SymmetryClass>
constexpr unsigned rotationMultiplicity(
  const unsigned& rotationFunctionIndex
) {
  return rotationMultiplicityImpl<SymmetryClass>(
    rotationFunctionIndex,
    applyRotation<SymmetryClass>(
      iota<SymmetryClass::size>(),
      rotationFunctionIndex
    ),
    1
  );
}

template<typename SymmetryClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationMultiplicityImpl(
  const ArrayType<unsigned, SymmetryClass::size>& runningIndices,
  const unsigned& count
) {
  if(
    ConstexprMagic::arraysEqual(
      runningIndices,
      iota<SymmetryClass::size>()
    )
  ) {
    return count;
  } else {
    return rotationMultiplicityImpl<SymmetryClass, rotationFunctionIndex>(
      applyRotation<SymmetryClass>(
        runningIndices,
        rotationFunctionIndex
      ),
      count + 1
    );
  }
}

template<typename SymmetryClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationMultiplicity() {
  return rotationMultiplicityImpl<SymmetryClass, rotationFunctionIndex>(
    applyRotation<SymmetryClass>(
      iota<SymmetryClass::size>(),
      rotationFunctionIndex
    ),
    1
  );
}

template<typename SymmetryClass, size_t ...Inds>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> 
rotationMultiplicitiesImpl(std::index_sequence<Inds...>) {
  return {
    rotationMultiplicity<SymmetryClass, Inds>()...
  };
}

template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> rotationMultiplicities() {
  return rotationMultiplicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
}

template<typename SymmetryClass>
struct allRotationMultiplicities {
  static constexpr ArrayType<
    unsigned,
    SymmetryClass::rotations
  > value = rotationMultiplicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
};


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
  const ArrayType<
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
  const ArrayType<unsigned, size>& indexMapping
) {
  if(symmetryPosition != replaceMe) {
    return indexMapping.at(symmetryPosition);
  }

  return replaceMe;
}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
constexpr double calculateChiralDistortion(
  const ArrayType<
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
constexpr ArrayType<unsigned, size> symPosMapping(
  const ArrayType<unsigned, size>& mapping
) {
  ArrayType<unsigned, size> symmetryPositions;

  for(unsigned i = 0; i < size; ++i) {
    symmetryPositions.at(
      mapping.at(i)
    ) = i;
  }

  return symmetryPositions;
}

template<size_t size>
struct DistortionInfo {
  ArrayType<unsigned, size> indexMapping;
  double totalDistortion;
  double chiralDistortion;

  constexpr DistortionInfo(
    const ArrayType<unsigned, size>& passIndexMapping,
    const double& passTotalDistortion,
    const double& passChiralDistortion
  ) : indexMapping(passIndexMapping),
      totalDistortion(passTotalDistortion),
      chiralDistortion(passChiralDistortion)
  {}
};

template<class SymmetryClassFrom, class SymmetryClassTo>
constexpr auto ligandGainMappings() {
  static_assert(
    SymmetryClassTo::size == SymmetryClassFrom::size + 1,
    "Ligand gain pathway calculation must involve symmetry size increase"
  );

  auto indexMapping  = iota<SymmetryClassTo::size>();

  ArrayType<
    ArrayType<unsigned, SymmetryClassTo::size>,
    0
  > bestMappings;
  double lowestAngleDistortion = 100;

  // No way to do nextPermutation? need swap

}

template<typename SymmetryClass>
using IndicesList = ArrayType<unsigned, SymmetryClass::size>;

template<typename SymmetryClass, size_t size>
constexpr std::enable_if_t<
  size == 0,
  std::pair<
    ArrayType<IndicesList<SymmetryClass>, 1>,
    ArrayType<unsigned, 1>
  >
> collapseChain(
  const ArrayType<IndicesList<SymmetryClass>, size>& chainStructures __attribute__((unused)),
  const ArrayType<unsigned, size>& chain __attribute__((unused))
) {
  return {};
}

template<typename SymmetryClass, size_t size>
constexpr auto collapseChain(
  const ArrayType<IndicesList<SymmetryClass>, size>& chainStructures,
  const ArrayType<unsigned, size>& chain
) {
  if(size == 1 || chain.back() < SymmetryClass::rotations.size() - 1) {
    // Increment last item in chain via pop and push
    const auto lastItem = chain.back();
    const auto intermediateChain = ConstexprMagic::arrayPop(chain);

    return std::make_pair(
      chainStructures,
      ConstexprMagic::arrayPush(intermediateChain, lastItem + 1)
    );
  }

  return collapseChain<SymmetryClass>(
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
  const ArrayType<IndicesList<SymmetryClass>, chainLength>& chainStructures,
  const ArrayType<unsigned, chainLength>& chain
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

constexpr unsigned subtractOne(const unsigned& x) {
  return x - 1;
}

template<typename SymmetryClass>
constexpr auto generateAllChains() {
  constexpr auto symmetryRotationMultiplicites = rotationMultiplicities<SymmetryClass>();

  constexpr auto multiplicitiesMinusOne = ConstexprMagic::map(
    symmetryRotationMultiplicites, 
    subtractOne
  );

  constexpr auto chainLength = ConstexprMagic::sum(multiplicitiesMinusOne);

  constexpr auto nCombinations = ConstexprMagic::Math::factorial(chainLength) / ConstexprMagic::sum(
    ConstexprMagic::map(
      multiplicitiesMinusOne,
      ConstexprMagic::Math::factorial<unsigned>
    )
  );

  ArrayType<
    ArrayType<unsigned, chainLength>,
    nCombinations
  > returnArray {};

  auto runningPermutation = ArrayType<unsigned, chainLength> {};
  unsigned counter = 0;
  for(unsigned i = 0; i < SymmetryClass::rotations.size(); ++i) {
    for(unsigned j = 0; j < multiplicitiesMinusOne.at(i); ++j) {
      runningPermutation.at(counter + j) = i;
    }
    counter += multiplicitiesMinusOne.at(i);
  }

  unsigned permutationIndex = 0;
  do {
    returnArray.at(permutationIndex) = runningPermutation;
    permutationIndex += 1;
  } while(ConstexprMagic::inPlaceNextPermutation(runningPermutation));

  return returnArray;
}

constexpr auto linearChains = generateAllChains<data::Linear>();

static_assert(linearChains.size() == 1, "Unexpected result");
static_assert(linearChains.at(0).size() == 1, "Unexpected result");

constexpr auto tetrahedralChains = generateAllChains<data::Tetrahedral>();

template<typename SymmetryClass>
constexpr auto generateAllRotations() {
  constexpr auto chains = generateAllChains<SymmetryClass>();
}

/*! It works in principle, but template instantiation depth is a serious issue!
 * Finding all 5040 rotations of SquareAntiprismatic is pretty daunting as at
 * least that many different instantiations of generateAllRotationsImpl have to
 * be made, compiled and evaluated at compile-time. Perhaps it is better to
 * avoid set generation and checking here.
 */
template<typename SymmetryClass>
constexpr auto generateAllRotations(
  const ArrayType<unsigned, SymmetryClass::size>& indices
) {
  ConstexprMagic::Set<
    IndicesList<SymmetryClass>
  > setInstance;

  return generateAllRotationsImpl<SymmetryClass>(
    setInstance,
    ArrayType<IndicesList<SymmetryClass>, 1> {{indices}},
    ArrayType<unsigned, 1> {0}
  );
}

/*constexpr auto maxAsymTetrRots = generateAllRotations<data::Linear>(
  ArrayType<unsigned, 2> {{0, 1}}
);*/

} // namespace Symmetry
