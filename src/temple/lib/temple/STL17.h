#ifndef INCLUDE_TEMPLE_STL_17
#define INCLUDE_TEMPLE_STL_17

#include <type_traits>
#include <cassert>
#include <functional>

namespace temple {

namespace stl17 {

// From cppreference, possible C++17 clamp implementation
template<class T, class Compare>
constexpr const T& clamp( const T& v, const T& lo, const T& hi, Compare comp ) {
    return assert( !comp(hi, lo) ),
        comp(v, lo) ? lo : comp(hi, v) ? hi : v;
}

template<class T>
constexpr const T& clamp( const T& v, const T& lo, const T& hi ) {
    return clamp( v, lo, hi, std::less<T>() );
}

template <class T>
constexpr std::add_const_t<T>& as_const(T& t) noexcept {
    return t;
}

} // namespace stl17

} // namespace temple

#endif
