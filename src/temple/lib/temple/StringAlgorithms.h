#ifndef INCLUDE_TEMPLE_STRING_ALGORITHMS_H
#define INCLUDE_TEMPLE_STRING_ALGORITHMS_H

#include <string>
#include <sstream>
#include <vector>

namespace temple {

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

} // namespace temple

#endif
