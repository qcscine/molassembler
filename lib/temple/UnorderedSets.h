#ifndef INCLUDE_TEMPLATE_MAGIC_UNORDERED_SET_ALGORITHMS_H
#define INCLUDE_TEMPLATE_MAGIC_UNORDERED_SET_ALGORITHMS_H

/*! @file
 *
 * Set-level algorithms for unordered sets.
 */

namespace temple {

template<class UnorderedSetContainer>
UnorderedSetContainer unorderedSetSymmetricDifference(
  const UnorderedSetContainer& a,
  const UnorderedSetContainer& b
) {
  UnorderedSetContainer symmetricDifference;

  for(const auto& element : a) {
    if(b.count(element) == 0) {
      symmetricDifference.insert(element);
    }
  }

  for(const auto& element : b) {
    if(a.count(element) == 0) {
      symmetricDifference.insert(element);
    }
  }

  return symmetricDifference;
}

template<class UnorderedSetContainer>
UnorderedSetContainer unorderedSetDifference(
  const UnorderedSetContainer& a,
  const UnorderedSetContainer& b
) {
  UnorderedSetContainer difference;

  for(const auto& element : a) {
    if(b.count(element) == 0) {
      difference.insert(element);
    }
  }

  return difference;
}

template<class UnorderedSetContainer>
UnorderedSetContainer unorderedSetIntersection(
  const UnorderedSetContainer& a,
  const UnorderedSetContainer& b
) {
  UnorderedSetContainer intersection;

  for(const auto& element : a) {
    if(b.count(element) == 1) {
      intersection.insert(element);
    }
  }

  return intersection;
}

template<class UnorderedSetContainer>
UnorderedSetContainer unorderedSetUnion(
  const UnorderedSetContainer& a,
  const UnorderedSetContainer& b
) {
  UnorderedSetContainer unionSet = a;

  for(const auto& element : b) {
    if(unionSet.count(element) == 0) {
      unionSet.insert(element);
    }
  }

  return unionSet;
}


} // namespace temple

#endif
