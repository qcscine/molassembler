/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Base 64 encoding and decoding between strings and vectors of uint_8
 */

#ifndef INCLUDE_BASE_64_ENCODING_H
#define INCLUDE_BASE_64_ENCODING_H

#include <string>
#include <vector>

namespace Scine {
namespace base64 {

/** @brief Encode binary data as a string
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param inputBuffer binary data
 */
std::string encode(const std::vector<std::uint8_t>& inputBuffer);

/** @brief Decode base 64 string data to binary
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param input base 64 string
 */
std::vector<std::uint8_t> decode(const std::string& input);

} // namespace base64
} // namespace Scine

#endif
