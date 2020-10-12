/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief String-specific algorithms
 */

#ifndef INCLUDE_TEMPLE_STRING_ALGORITHMS_H
#define INCLUDE_TEMPLE_STRING_ALGORITHMS_H

#include <string>
#include <sstream>
#include <vector>

namespace Scine {
namespace Molassembler {
namespace Temple {

inline std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> elements;

  std::stringstream ss;
  ss.str(s);
  std::string item;

  while(std::getline(ss, item, delimiter)) {
    elements.push_back(item);
  }

  return elements;
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
