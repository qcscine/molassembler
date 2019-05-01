/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Useful functions for dealing with STL containers / types.
 */

#ifndef INCLUDE_MOLASSEMBLER_STDLIB_TYPE_ALGORITHMS_H
#define INCLUDE_MOLASSEMBLER_STDLIB_TYPE_ALGORITHMS_H

#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <map>

/* TODO
 * - Figure out if this is used somewhere and get rid of all of these algorithms
 */

namespace Scine {

namespace molassembler {

namespace StdlibTypeAlgorithms {

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

void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

} // namespace StdlibTypeAlgorithms

} // namespace molassembler

} // namespace Scine

#endif
