#ifndef INCLUDE_DG_SEQUENCE_SET
#define INCLUDE_DG_SEQUENCE_SET

#include "common_typedefs.h"

#include <set>
#include <array>

/* NOTES
 * - As soon as larger molecules need to be treated, add a bloom filter for set
 *   contains checking speed improvement. Set membership via count is O(log N),
 *   while using e.g. a Bloom filter would be O(1). Consider also using a hash
 *   table or so
 * - This implementation requires a lot of copying / reversing, which could
 *   be avoided, maybe by e.g. using forward and backwards comparison operators?
 *   Unsure if comparing twice is cheaper than reversing the trial sequence
 * - using different types, i.e. vector and array, causes a lot of
 *   type-conversion copying! Custom comparison operators could get rid of that
 */

namespace MoleculeManip {

namespace DistanceGeometry {

template<int size>
void standardize(std::array<AtomIndexType, size>& sequence) {
  if(sequence.front() > sequence.back()) {
    std::reverse(
      sequence.begin(),
      sequence.end()
    );
  }

  return sequence;
}

template<int sequenceLength>
class SequenceSet {
public:
  using SequenceSetType = std::set<
    std::array<AtomIndexType, sequenceLength>
  >;

private:
/* State */
  SequenceSetType _sequences;

/* Static private members */
  static void _standardizeInplace(std::array<AtomIndexType, sequenceLength>& sequence) {
    if(sequence.front() > sequence.back()) {
      std::reverse(
        sequence.begin(),
        sequence.end()
      );
    }

    return sequence;
  }

  static std::array<AtomIndexType, sequenceLength> _standardize(const std::array<AtomIndexType, sequenceLength>& sequence) {
    if(sequence.front() < sequence.back()) {
      return sequence;
    }

    std::array<AtomIndexType, sequenceLength> reversed;

    for(unsigned i = 0; i < sequenceLength; i++) {
      reversed.at(i) = sequence.at(sequenceLength - i - 1);
    }

    return reversed;
  }

  static std::array<AtomIndexType, sequenceLength> _standardize(const std::vector<AtomIndexType>& sequence) {
    assert(sequence.size() == sequenceLength);

    if(sequence.front() < sequence.back()) {
      std::array<AtomIndexType, sequenceLength> arraySequence;

      for(unsigned i = 0; i < sequenceLength; i++) {
        arraySequence.at(i) = sequence.at(i);
      }

      return arraySequence;
    }

    std::array<AtomIndexType, sequenceLength> reversed;

    for(unsigned i = 0; i < sequenceLength; i++) {
      reversed.at(i) = sequence.at(sequenceLength - i - 1);
    }

    return reversed;
  }

public:
/* Constructors */
  SequenceSet() {}

/* Modification */
  void insert(const std::array<AtomIndexType, sequenceLength>& sequence) {
    _sequences.insert(
      _standardize(sequence)
    );
  }

  void insert(const std::vector<AtomIndexType>& sequence) {
    _sequences.insert(
      _standardize(sequence)
    );
  }

/* Information */
  bool contains(const std::array<AtomIndexType, sequenceLength>& sequence) const {
    return _sequences.count(
      _standardize(sequence)
    ) > 0;
  }

  bool contains(const std::vector<AtomIndexType>& sequence) const {
    return _sequences.count(
      _standardize(sequence)
    ) > 0;
  }

  unsigned size() const {
    return _sequences.size();
  }

  typename SequenceSetType::const_iterator begin() const {
    return _sequences.cbegin();
  }

  typename SequenceSetType::const_iterator end() const {
    return _sequences.cend();
  }
};

} // namespace MoleculeManip

} // namespace DistanceGeometry

#endif
