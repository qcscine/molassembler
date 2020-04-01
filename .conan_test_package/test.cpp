#include "Utils/Geometry/ElementInfo.h"

using namespace Scine::Utils;

int main() {
  if(ElementInfo::element(1) == ElementType::H) {
    return 0;
  }

  return 1;
}
