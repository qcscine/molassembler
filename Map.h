#ifndef INCLUDE_CONSTEXPR_MAGIC_MAP_H
#define INCLUDE_CONSTEXPR_MAGIC_MAP_H

#include <array>

namespace ConstexprMagic { 

// This has NONE of the desired properties of a map. Lookup is linear!
template<typename KeyType, typename MappedType, size_t size = 0>
class Map {
private:
  std::array<KeyType, size> _keys;
  std::array<MappedType, size> _items;

public:
  constexpr Map() {}
  constexpr Map(
    std::array<KeyType, size>&& keys,
    std::array<MappedType, size>&& items
  ) : _keys(keys),
      _items(items)
  {}

  constexpr Map<KeyType, MappedType, (size + 1)> insert(
    KeyType key,
    MappedType item
  ) const {
    return Map<KeyType, MappedType, (size + 1)>(
      arrayPush(_keys, key),
      arrayPush(_items, item)
    );
  }

  constexpr MappedType at(const KeyType& key) const {
    for(unsigned i = 0; i < size; ++i) {
      if(_keys.at(i) == key) {
        return _keys.at(i);
      }
    }

    // If not found, returns default-constructed MappedType
    return {};
  }

  constexpr bool contains(const KeyType& key) const {
    for(unsigned i = 0; i < size; ++i) {
      if(_keys.at(i) == key) {
        return true;
      }
    }

    return false;
  }
};

/*constexpr auto testMap = Map<unsigned, double>();
constexpr auto withAnElement = testMap.insert(5, 1.3);
constexpr auto getElement = testMap.at(5);*/

} // namespace ConstexprMagic

#endif
