#ifndef INCLUDE_TESTING_HELPER_MAKE_REPEATED_ELEMENT_COLLECTION
#define INCLUDE_TESTING_HELPER_MAKE_REPEATED_ELEMENT_COLLECTION

#include "Types/ElementTypeCollection.h"

Delib::ElementTypeCollection makeRepeatedElementCollection(
  const Delib::ElementType& elementType,
  const unsigned& repeat
) {
  Delib::ElementTypeCollection returnCollection;

  for(unsigned i = 0; i < repeat; i++) {
    returnCollection.push_back(elementType);
  }

  return returnCollection;
}

Delib::ElementTypeCollection makeIncrementedElementCollection(
  const unsigned& size
) {
  Delib::ElementTypeCollection returnCollection;

  for(unsigned i = 0; i < size; i++) {
    returnCollection.push_back(Delib::ElementType(i + 1));
  }

  return returnCollection;
}

#endif
