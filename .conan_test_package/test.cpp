#include "molassembler/Types.h"

using namespace Scine::Molassembler;

int main() {
  BondIndex f {4, 3};

  if(f.first == 3) {
    return 0;
  }

  return 1;
}
