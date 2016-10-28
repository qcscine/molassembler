#ifndef INCLUDE_MOLLIB_VERSION_H
#define INCLUDE_MOLLIB_VERSION_H

#include <string>

namespace Version {

unsigned short major = 0;
unsigned short minor = 1;

std::string String() {
  return std::to_string(major)+"."+std::to_string(minor);
}

}

#endif
