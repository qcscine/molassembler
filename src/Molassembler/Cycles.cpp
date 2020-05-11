/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Cycles.h"

#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Graph.h"

#include "Molassembler/Temple/TinySet.h"
#include "Molassembler/Temple/Functional.h"

#include "RingDecomposerLib.h"
#include "boost/variant.hpp"

namespace Scine {
namespace Molassembler {
namespace {

inline std::vector<unsigned> intersect(
  const std::vector<unsigned>& a,
  const std::vector<unsigned>& b
) {
  assert(std::is_sorted(std::begin(a), std::end(a)));
  assert(std::is_sorted(std::begin(b), std::end(b)));
  std::vector<unsigned> intersection;
  intersection.reserve(
    std::min(a.size(), b.size())
  );

  std::set_intersection(
    std::begin(a),
    std::end(a),
    std::begin(b),
    std::end(b),
    std::back_inserter(intersection)
  );

  return intersection;
}

} // namespace

/* Cycles member Type declarations */
struct Cycles::RdlDataPtrs {
  //! Raw pointers to RDL_graph object
  RDL_graph* graphPtr;

  //! Raw pointer to calculated graph cycle data
  RDL_data* dataPtr;

  RdlDataPtrs() = delete;
  RdlDataPtrs(const PrivateGraph& sourceGraph, bool ignoreEtaBonds);

  RdlDataPtrs(const RdlDataPtrs& other) = delete;
  RdlDataPtrs(RdlDataPtrs&& other) = delete;
  RdlDataPtrs& operator = (const RdlDataPtrs& other) = delete;
  RdlDataPtrs& operator = (RdlDataPtrs&& other) = delete;
  ~RdlDataPtrs();

  bool bondExists(const BondIndex& bond) const;
};

Cycles::Cycles(const Graph& sourceGraph, const bool ignoreEtaBonds)
  : Cycles {sourceGraph.inner(), ignoreEtaBonds}
{}

Cycles::Cycles(const PrivateGraph& sourceGraph, const bool ignoreEtaBonds)
  : rdlPtr_(std::make_shared<RdlDataPtrs>(sourceGraph, ignoreEtaBonds))
{
  unsigned U = RDL_getNofURF(rdlPtr_->dataPtr);
  for(unsigned i = 0; i < U; ++i) {
    RDL_edge* edgeArray;
    unsigned nEdges = RDL_getEdgesForURF(rdlPtr_->dataPtr, i, &edgeArray);
    if(nEdges == RDL_INVALID_RESULT) {
      throw std::logic_error("RDL failure");
    }

    for(unsigned j = 0; j < nEdges; ++j) {
      const BondIndex bond { edgeArray[j][0], edgeArray[j][1] };
      urfMap_[bond].push_back(i);
    }
  }

  // Sort the URFs for set intersections
  for(auto& entry : urfMap_) {
    std::sort(
      std::begin(entry.second),
      std::end(entry.second)
    );
  }
}

unsigned Cycles::numCycleFamilies() const {
  return RDL_getNofURF(rdlPtr_->dataPtr);
}

unsigned Cycles::numCycleFamilies(const AtomIndex index) const {
  return RDL_getNofURFContainingNode(rdlPtr_->dataPtr, index);
}

//! Returns the number of relevant cycles (RCs)
unsigned Cycles::numRelevantCycles() const {
  return RDL_getNofRC(rdlPtr_->dataPtr);
}

unsigned Cycles::numRelevantCycles(const AtomIndex index) const {
  return RDL_getNofRCFContainingNode(rdlPtr_->dataPtr, index);
}

Cycles::AllCyclesIterator Cycles::begin() const {
  return {rdlPtr_};
}

Cycles::AllCyclesIterator Cycles::end() const {
  return {
    rdlPtr_,
    numRelevantCycles()
  };
}

/* Cycles::RdlCyclePtrs */
struct Cycles::RdlCyclePtrs {
  RDL_cycleIterator* cycleIterPtr = nullptr;
  RDL_cycle* cyclePtr = nullptr;
  unsigned rCycleIndex = 0;
  std::vector<BondIndex> bonds;

  RdlCyclePtrs() = delete;

  //! Construct cycle ptrs just from cycle data (get all RCycles)
  RdlCyclePtrs(const RdlDataPtrs& dataPtrs) {
    cycleIterPtr = RDL_getRCyclesIterator(dataPtrs.dataPtr);
    initializeCycle();
  }

  //! Construct cycle ptrs from a particular URF ID
  RdlCyclePtrs(const RdlDataPtrs& dataPtrs, const unsigned URFID) {
    cycleIterPtr = RDL_getRCyclesForURFIterator(dataPtrs.dataPtr, URFID);
    initializeCycle();
  }

  RdlCyclePtrs(const RdlCyclePtrs& other) = delete;
  RdlCyclePtrs(RdlCyclePtrs&& other) = delete;
  RdlCyclePtrs& operator = (const RdlCyclePtrs& other) = delete;
  RdlCyclePtrs& operator = (RdlCyclePtrs&& other) = delete;

  bool atEnd() const {
    return RDL_cycleIteratorAtEnd(cycleIterPtr) != 0;
  }

  bool operator == (const RdlCyclePtrs& other) const {
    if(cyclePtr == nullptr && other.cyclePtr == nullptr) {
      return true;
    }

    if(cyclePtr != nullptr && other.cyclePtr != nullptr) {
      return (
        std::tie(cyclePtr->weight, cyclePtr->rcf, cyclePtr->urf)
        == std::tie(other.cyclePtr->weight, other.cyclePtr->rcf, other.cyclePtr->urf)
      );
    }

    // Mixed case remains
    return false;
  }

  bool operator != (const RdlCyclePtrs& other) const {
    return !(*this == other);
  }


  ~RdlCyclePtrs() {
    if(cyclePtr != nullptr) {
      RDL_deleteCycle(cyclePtr);
      cyclePtr = nullptr;
    }

    RDL_deleteCycleIterator(cycleIterPtr);
  }

  void populateBondsForCurrentCycle() {
    assert(cyclePtr != nullptr);
    bonds.clear();

    const unsigned size = cyclePtr->weight;
    bonds.reserve(size);

    for(unsigned i = 0; i < size; ++i) {
      bonds.emplace_back(
        cyclePtr->edges[i][0],
        cyclePtr->edges[i][1]
      );
    }
  }

  void initializeCycle() {
    if(!atEnd()) {
      cyclePtr = RDL_cycleIteratorGetCycle(cycleIterPtr);
      populateBondsForCurrentCycle();
    }
  }

  /*!
   * @brief Advance internal iterator and cycle state
   *
   * Frees the memory for the current cycle and advances the iterator state.
   * If the iterator is now not at the end of all relevant cycles, then
   * the next cycle is allocated and constructed into cyclePtr. Otherwise,
   * cyclePtr is a nullptr.
   */
  void advance() {
    assert(!atEnd());

    RDL_deleteCycle(cyclePtr);
    cyclePtr = nullptr;

    RDL_cycleIteratorNext(cycleIterPtr);
    ++rCycleIndex;

    if(!atEnd()) {
      cyclePtr = RDL_cycleIteratorGetCycle(cycleIterPtr);
      populateBondsForCurrentCycle();
    }
  }

  void matchCycleState(const RdlCyclePtrs& other) {
    assert(rCycleIndex == 0);
    for(unsigned i = 0; i < other.rCycleIndex; ++i) {
      advance();
    }

    assert(rCycleIndex == other.rCycleIndex);
  }
};

namespace {

template<typename Arg>
IteratorRange<Cycles::UrfIdsCycleIterator>
makeURFIDsCycleIterator(const std::shared_ptr<Cycles::RdlDataPtrs>& dataPtr, Arg&& arg) {
  Cycles::UrfIdsCycleIterator begin {arg, dataPtr};
  Cycles::UrfIdsCycleIterator end = begin;
  end.advanceToEnd();

  return {
    std::move(begin),
    std::move(end)
  };
}

} // namespace

IteratorRange<Cycles::UrfIdsCycleIterator>
Cycles::containing(const AtomIndex atom) const {
  return makeURFIDsCycleIterator(rdlPtr_, atom);
}

IteratorRange<Cycles::UrfIdsCycleIterator>
Cycles::containing(const BondIndex& bond) const {
  return containing(std::vector<BondIndex> {bond});
}

IteratorRange<Cycles::UrfIdsCycleIterator>
Cycles::containing(const std::vector<BondIndex>& bonds) const {
  auto fetchBondURFs = [&](const BondIndex& bond) -> std::vector<unsigned> {
    auto findIter = urfMap_.find(bond);
    if(findIter == std::end(urfMap_)) {
      return {};
    }

    return findIter->second;
  };

  std::vector<unsigned> urfs = fetchBondURFs(bonds.front());
  for(unsigned i = 1; i < bonds.size() && !urfs.empty(); ++i) {
    urfs = intersect(urfs, fetchBondURFs(bonds.at(i)));
  }

  Cycles::UrfIdsCycleIterator begin {bonds, std::move(urfs), rdlPtr_};
  Cycles::UrfIdsCycleIterator end = begin;
  end.advanceToEnd();

  return {
    std::move(begin),
    std::move(end)
  };
}

RDL_data* Cycles::dataPtr() const {
  return rdlPtr_->dataPtr;
}

bool Cycles::operator == (const Cycles& other) const {
  return rdlPtr_ == other.rdlPtr_;
}

bool Cycles::operator != (const Cycles& other) const {
  return !(*this == other);
}

/* Cycles::RdlDataPtrs */
Cycles::RdlDataPtrs::RdlDataPtrs(
  const PrivateGraph& sourceGraph,
  const bool ignoreEtaBonds
) {
  // Initialize a new graph
  graphPtr = RDL_initNewGraph(sourceGraph.N());

  if(ignoreEtaBonds) {
    for(const auto edge : sourceGraph.edges()) {
      // Copy vertices (without bond type information)
      if(sourceGraph.bondType(edge) != BondType::Eta) {
        auto edgeAddResult = RDL_addUEdge(
          graphPtr,
          sourceGraph.source(edge),
          sourceGraph.target(edge)
        );

        if(edgeAddResult == RDL_INVALID_RESULT || edgeAddResult == RDL_DUPLICATE_EDGE) {
          throw std::runtime_error("RDL add edge failed!");
        }
      }
    }
  } else {
    for(const auto edge : sourceGraph.edges()) {
      auto edgeAddResult = RDL_addUEdge(
        graphPtr,
        sourceGraph.source(edge),
        sourceGraph.target(edge)
      );

      if(edgeAddResult == RDL_INVALID_RESULT || edgeAddResult == RDL_DUPLICATE_EDGE) {
        throw std::runtime_error("RDL add edge failed!");
      }
    }
  }

  // Calculate, and assure yourself of success
  dataPtr = RDL_calculate(graphPtr);
  assert(dataPtr != nullptr);
}

Cycles::RdlDataPtrs::~RdlDataPtrs() {
  // Awfully enough, calling this frees both dataPtr and graphPtr
  RDL_deleteData(dataPtr);
}

bool Cycles::RdlDataPtrs::bondExists(const BondIndex& bond) const {
  const unsigned bondID = RDL_getEdgeId(dataPtr, bond.first, bond.second);
  return bondID != RDL_INVALID_RESULT;
}

Cycles::AllCyclesIterator::AllCyclesIterator(const AllCyclesIterator& other)
  : rdlPtr_(other.rdlPtr_),
    cyclePtr_(std::make_unique<RdlCyclePtrs>(*rdlPtr_))
{
  cyclePtr_->matchCycleState(*other.cyclePtr_);
}

Cycles::AllCyclesIterator::AllCyclesIterator(AllCyclesIterator&& other) noexcept
  : rdlPtr_(std::move(other.rdlPtr_)),
    cyclePtr_(std::move(other.cyclePtr_))
{}

Cycles::AllCyclesIterator& Cycles::AllCyclesIterator::operator = (const AllCyclesIterator& other) {
  rdlPtr_ = other.rdlPtr_;
  cyclePtr_ = std::make_unique<RdlCyclePtrs>(*rdlPtr_);
  cyclePtr_->matchCycleState(*other.cyclePtr_);
  return *this;
}

Cycles::AllCyclesIterator& Cycles::AllCyclesIterator::operator = (AllCyclesIterator&& other) noexcept {
  rdlPtr_ = std::move(other.rdlPtr_);
  cyclePtr_ = std::move(other.cyclePtr_);
  return *this;
}

Cycles::AllCyclesIterator::~AllCyclesIterator() = default;

Cycles::AllCyclesIterator::AllCyclesIterator(
  const std::shared_ptr<RdlDataPtrs>& dataPtr,
  unsigned rCycleIndex
) : rdlPtr_ {dataPtr} {
  cyclePtr_ = std::make_unique<RdlCyclePtrs>(*dataPtr);

  if(rCycleIndex != 0) {
    // Constructor to advance to specific index
    while(cyclePtr_->rCycleIndex < rCycleIndex) {
      ++(*this);
    }

    assert(cyclePtr_->rCycleIndex == rCycleIndex);
  }
}

Cycles::AllCyclesIterator& Cycles::AllCyclesIterator::operator ++ () {
  if(!cyclePtr_->atEnd()) {
    cyclePtr_->advance();
  } else {
    throw std::logic_error("Advancing Cycles::AllCyclesIterator past end");
  }

  return *this;
}

Cycles::AllCyclesIterator Cycles::AllCyclesIterator::operator ++ (int) {
  AllCyclesIterator retval = *this;
  ++(*this);
  return retval;
}

Cycles::AllCyclesIterator::value_type Cycles::AllCyclesIterator::operator * () const {
  if(cyclePtr_->cyclePtr == nullptr) {
    throw std::range_error("Dereferencing Cycles::AllCyclesIterator end iterator");
  }

  assert(!cyclePtr_->bonds.empty());
  return cyclePtr_->bonds;
}

Cycles::AllCyclesIterator::pointer Cycles::AllCyclesIterator::operator -> () const {
  return &(this->operator *());
}

//! Must be constructed from same Cycles base and at same RC
bool Cycles::AllCyclesIterator::operator == (const Cycles::AllCyclesIterator& other) const {
  return (
    rdlPtr_ == other.rdlPtr_
    && cyclePtr_->rCycleIndex == other.cyclePtr_->rCycleIndex
  );
}

bool Cycles::AllCyclesIterator::operator != (const Cycles::AllCyclesIterator& other) const {
  return !(*this == other);
}

struct Cycles::UrfIdsCycleIterator::UrfHelper {
  using VariantType = boost::variant<AtomIndex, std::vector<BondIndex>>;

  template<typename RdlFunc, typename ... Args>
  static std::vector<unsigned> getURFsHelper(
    const RdlDataPtrs& dataPtrs,
    RdlFunc&& RDL_function,
    Args ... args
  ) {
    unsigned* idsPtr;

    // Call the RDL function, allocating memory at idsPtr
    const unsigned idCount = RDL_function(
      dataPtrs.dataPtr,
      args ...,
      &idsPtr
    );

    if(idCount == RDL_INVALID_RESULT) {
      throw std::logic_error("RDL returned invalid result to URF ID call");
    }

    std::vector<unsigned> ids(idCount);
    for(unsigned i = 0; i < idCount; ++i) {
      ids[i] = idsPtr[i];
    }

    // Free memory allocated by RDL
    free(idsPtr);
    return ids;
  }

  static std::vector<unsigned> getURFs(
    AtomIndex atom,
    const RdlDataPtrs& dataPtrs
  ) {
    return getURFsHelper(
      dataPtrs,
      &RDL_getURFsContainingNode,
      atom
    );
  }

  static std::vector<unsigned> getURFs(
    const BondIndex& bond,
    const RdlDataPtrs& dataPtrs
  ) {
    assert(dataPtrs.bondExists(bond));

    return getURFsHelper(
      dataPtrs,
      &RDL_getURFsContainingEdge,
      bond.first,
      bond.second
    );
  }

  static std::vector<unsigned> getURFs(
    const std::vector<BondIndex>& bonds,
    const RdlDataPtrs& dataPtrs
  ) {
    assert(!bonds.empty());
    std::vector<unsigned> idsIntersection = getURFs(bonds.front(), dataPtrs);

    if(idsIntersection.empty()) {
      return idsIntersection;
    }

    Temple::InPlace::sort(idsIntersection);

    // Continually intersect idsIntersection with the URF IDs of the next bond
    for(unsigned i = 1; i < bonds.size(); ++i) {
      auto newIDs = getURFs(bonds[i], dataPtrs);
      Temple::InPlace::sort(newIDs);

      idsIntersection = intersect(idsIntersection, newIDs);

      if(idsIntersection.empty()) {
        break;
      }
    }

    return idsIntersection;
  }

  struct CycleSatisfiesSoughtVisitor : public boost::static_visitor<bool> {
    const std::vector<BondIndex>& cycleReference;

    CycleSatisfiesSoughtVisitor(const std::vector<BondIndex>& cycleEdges)
      : cycleReference(cycleEdges) {}


    bool operator() (const AtomIndex atom) const {
      return Temple::any_of(
        cycleReference,
        [&atom](const BondIndex& bond) -> bool {
          return bond.contains(atom);
        }
      );
    }

    bool operator() (const std::vector<BondIndex>& bonds) const {
      return Temple::all_of(
        bonds,
        Temple::makeContainsPredicate(cycleReference)
      );
    }
  };

  bool cycleSatisfiesSoughtConditions(const std::vector<BondIndex>& bonds) {
    CycleSatisfiesSoughtVisitor visitor {bonds};
    return boost::apply_visitor(visitor, soughtVariant_);
  }

  template<typename Arg>
  UrfHelper(Arg&& arg, const RdlDataPtrs& dataPtrs)
    : soughtVariant_(arg),
      ids(getURFs(arg, dataPtrs)),
      idsIdx(0)
  {}

  template<typename Arg>
  UrfHelper(Arg&& arg, std::vector<unsigned> urfs)
    : soughtVariant_(arg),
      ids(std::move(urfs)),
      idsIdx(0)
  {}

  VariantType soughtVariant_;
  std::vector<unsigned> ids;
  unsigned idsIdx;
};

Cycles::UrfIdsCycleIterator::UrfIdsCycleIterator(
  AtomIndex soughtIndex,
  const std::shared_ptr<RdlDataPtrs>& dataPtr
) : rdlPtr_(dataPtr),
    urfsPtr_(std::make_unique<UrfHelper>(soughtIndex, *dataPtr)),
    cyclePtr_()
{
  initializeCyclesFromURFID_();
}

Cycles::UrfIdsCycleIterator::UrfIdsCycleIterator(
  const BondIndex& soughtBond,
  std::vector<unsigned> urfs,
  const std::shared_ptr<RdlDataPtrs>& dataPtr
) : rdlPtr_(dataPtr),
    urfsPtr_(
      std::make_unique<UrfHelper>(
        std::vector<BondIndex> {soughtBond},
        std::move(urfs)
      )
    ),
    cyclePtr_()
{
  initializeCyclesFromURFID_();
}

Cycles::UrfIdsCycleIterator::UrfIdsCycleIterator(
  const std::vector<BondIndex>& soughtBonds,
  std::vector<unsigned> urfs,
  const std::shared_ptr<RdlDataPtrs>& dataPtr
) : rdlPtr_(dataPtr),
    urfsPtr_(
      std::make_unique<UrfHelper>(soughtBonds, std::move(urfs))
    ),
    cyclePtr_()
{
  initializeCyclesFromURFID_();
}

Cycles::UrfIdsCycleIterator::UrfIdsCycleIterator(const UrfIdsCycleIterator& other)
  : rdlPtr_(other.rdlPtr_),
    urfsPtr_(std::make_unique<UrfHelper>(*other.urfsPtr_)),
    cyclePtr_()
{
  matchCycleState_(other);
}

Cycles::UrfIdsCycleIterator::UrfIdsCycleIterator(UrfIdsCycleIterator&& other) noexcept
  : rdlPtr_(std::move(other.rdlPtr_)),
    urfsPtr_(std::move(other.urfsPtr_)),
    cyclePtr_(std::move(other.cyclePtr_))
{}

Cycles::UrfIdsCycleIterator& Cycles::UrfIdsCycleIterator::operator = (const UrfIdsCycleIterator& other) {
  rdlPtr_ = other.rdlPtr_;
  *urfsPtr_ = *other.urfsPtr_;
  matchCycleState_(other);

  return *this;
}

void Cycles::UrfIdsCycleIterator::matchCycleState_(const UrfIdsCycleIterator& other) {
  if(urfsPtr_->idsIdx < urfsPtr_->ids.size()) {
    initializeCyclesFromURFID_();
  }

  if(other.cyclePtr_ && other.cyclePtr_->cyclePtr != nullptr) {
    while(
      cyclePtr_->cyclePtr != nullptr
      && *cyclePtr_ != *other.cyclePtr_
    ) {
      cyclePtr_->advance();
    }

    if(cyclePtr_->cyclePtr == nullptr) {
      throw std::runtime_error("Could not match state in copy!");
    }
  }
}

Cycles::UrfIdsCycleIterator& Cycles::UrfIdsCycleIterator::operator = (UrfIdsCycleIterator&& other) noexcept {
  rdlPtr_ = std::move(other.rdlPtr_);
  urfsPtr_ = std::move(other.urfsPtr_);
  cyclePtr_ = std::move(other.cyclePtr_);

  return *this;
}

Cycles::UrfIdsCycleIterator::~UrfIdsCycleIterator() = default;

void Cycles::UrfIdsCycleIterator::initializeCyclesFromURFID_() {
  assert(urfsPtr_);
  if(!urfsPtr_->ids.empty()) {
    assert(urfsPtr_->idsIdx < urfsPtr_->ids.size());
    cyclePtr_ = std::make_unique<RdlCyclePtrs>(
      *rdlPtr_,
      urfsPtr_->ids[urfsPtr_->idsIdx]
    );

    if(
      !cyclePtr_->atEnd()
      && !urfsPtr_->cycleSatisfiesSoughtConditions(cyclePtr_->bonds)
    ) {
      advanceToNextPermissibleCycle_();
    }

    /* This is odd, but there does NOT have to be at least one RC in each URF
     * that fulfills the sought conditions
     */
  }
}

void Cycles::UrfIdsCycleIterator::advanceToNextPermissibleCycle_() {
  assert(!cyclePtr_->atEnd());

  do {
    cyclePtr_->advance();
  } while(
    !cyclePtr_->atEnd()
    && !urfsPtr_->cycleSatisfiesSoughtConditions(cyclePtr_->bonds)
  );
}

Cycles::UrfIdsCycleIterator& Cycles::UrfIdsCycleIterator::operator ++ () {
  assert(!cyclePtr_->atEnd());

  // Advance the cycle iterator if it's not at the end yet
  if(!cyclePtr_->atEnd()) {
    advanceToNextPermissibleCycle_();
  }

  /* If the cycle iterator is at the end now, advance to the next URF ID from
   * the list if possible. It may be strange, but it is not necessary
   * that a URF ID yields any RCs that fulfill the sought conditions.
   */
  while(cyclePtr_->atEnd()) {
    ++urfsPtr_->idsIdx;
    if(urfsPtr_->idsIdx < urfsPtr_->ids.size()) {
      initializeCyclesFromURFID_();
    } else {
      cyclePtr_ = {};
      break;
    }
  }

  return *this;
}

Cycles::UrfIdsCycleIterator Cycles::UrfIdsCycleIterator::operator ++ (int) {
  UrfIdsCycleIterator thisCopy = *this;
  ++(*this);
  return thisCopy;
}

Cycles::UrfIdsCycleIterator::value_type Cycles::UrfIdsCycleIterator::operator * () const {
  if(cyclePtr_->cyclePtr == nullptr) {
    throw std::range_error("Dereferencing Cycles::UrfIdsCycleIterator end iterator");
  }

  assert(urfsPtr_->cycleSatisfiesSoughtConditions(cyclePtr_->bonds));
  assert(!cyclePtr_->bonds.empty());
  return cyclePtr_->bonds;
}

Cycles::UrfIdsCycleIterator::pointer Cycles::UrfIdsCycleIterator::operator -> () const {
  return &(this->operator *());
}

bool Cycles::UrfIdsCycleIterator::operator == (const UrfIdsCycleIterator& other) const {
  // Compare address of shared_ptr and the URF id index
  if(
    std::tie(rdlPtr_, urfsPtr_->idsIdx)
    != std::tie(other.rdlPtr_, other.urfsPtr_->idsIdx)
  ) {
    return false;
  }

  // cyclePtr_ is nullable (end iterator)
  if(cyclePtr_.operator bool() xor other.cyclePtr_.operator bool()) {
    // Heterogeneous case, one of both is an end iterator
    return false;
  }

  // Homogeneous case 1 (both are the end iterator)
  if(!cyclePtr_ && !other.cyclePtr_) {
    return true;
  }

  // Homogeneous case 2 (neither are the end iterator)
  return cyclePtr_->rCycleIndex == other.cyclePtr_->rCycleIndex;
}

bool Cycles::UrfIdsCycleIterator::operator != (const UrfIdsCycleIterator& other) const {
  return !(*this == other);
}

void Cycles::UrfIdsCycleIterator::advanceToEnd() {
  cyclePtr_.reset();
  urfsPtr_->idsIdx = urfsPtr_->ids.size();
}

boost::optional<unsigned> smallestCycleContaining(AtomIndex atom, const Cycles& cycles) {
  auto iteratorPair = cycles.containing(atom);

  auto minIter = std::min_element(
    iteratorPair.first,
    iteratorPair.second,
    [](const std::vector<BondIndex>& cycleA, const std::vector<BondIndex>& cycleB) -> bool {
      return cycleA.size() < cycleB.size();
    }
  );

  if(minIter == iteratorPair.second) {
    return boost::none;
  }

  return (*minIter).size();
}

std::unordered_map<AtomIndex, unsigned> makeSmallestCycleMap(const Cycles& cycleData) {
  std::unordered_map<AtomIndex, unsigned> smallestCycle;

  for(const auto& cycleEdges : cycleData) {
    const unsigned cycleSize = cycleEdges.size();

    for(const BondIndex& bond : cycleEdges) {
      for(const AtomIndex index: bond) {
        auto findIter = smallestCycle.find(index);

        if(findIter != std::end(smallestCycle)) {
          if(cycleSize < findIter->second) {
            findIter->second = cycleSize;
          }
        } else {
          smallestCycle.emplace(index, cycleSize);
        }
      }
    }
  }

  return smallestCycle;
}

std::vector<AtomIndex> makeRingIndexSequence(
  std::vector<BondIndex> edgeDescriptors
) {

  auto firstEdgeIter = edgeDescriptors.begin();
  // Initialize with first edge descriptor's indices
  std::vector<AtomIndex> indexSequence {
    firstEdgeIter->first,
    firstEdgeIter->second
  };

  indexSequence.reserve(edgeDescriptors.size());

  edgeDescriptors.erase(firstEdgeIter);
  // firstEdgeIter is now invalid!

  while(!edgeDescriptors.empty()) {
    for(
      auto edgeIter = edgeDescriptors.begin();
      edgeIter != edgeDescriptors.end();
      ++edgeIter
    ) {
      if(edgeIter->first == indexSequence.back()) {
        indexSequence.emplace_back(
          edgeIter->second
        );
        edgeDescriptors.erase(edgeIter);
        break;
      }

      if(edgeIter->second == indexSequence.back()) {
        indexSequence.emplace_back(
          edgeIter->first
        );
        edgeDescriptors.erase(edgeIter);
        break;
      }
    }
  }
  /* Now indexSequence should contain the entire sequence, but the first
   * vertex index occurs twice, once at the front and once at the back.
   *
   * Remove the repetition and then return.
   */
  indexSequence.pop_back();

  return indexSequence;
}

std::vector<AtomIndex> centralizeRingIndexSequence(
  std::vector<AtomIndex> ringIndexSequence,
  const AtomIndex center
) {
  auto findIter = std::find(
    std::begin(ringIndexSequence),
    std::end(ringIndexSequence),
    center
  );

  assert(findIter != std::end(ringIndexSequence));

  std::rotate(
    std::begin(ringIndexSequence),
    findIter,
    std::end(ringIndexSequence)
  );

  return ringIndexSequence;
}


unsigned countPlanarityEnforcingBonds(
  const std::vector<BondIndex>& edgeSet,
  const Graph& graph
) {
  return std::accumulate(
    std::begin(edgeSet),
    std::end(edgeSet),
    0u,
    [&graph](const unsigned carry, const BondIndex& edge) {
      if(graph.bondType(edge) == BondType::Double) {
        return carry + 1;
      }

      return carry;
    }
  );
}

} // namespace Molassembler
} // namespace Scine
