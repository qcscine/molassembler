/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 *
 * The implementation of this consumer-facing type is a little complicated, so
 * it warrants an explanation at length. Since molassembler mostly pImpls its
 * public types, this class was another opportunity to hide dependencies from
 * library consumers. But this one is a little tricky, because the boost graph
 * dependency brings with it vertex and edge descriptors which need to be freely
 * distributable to work with the graph itself.
 *
 * In order to hide those types, alternate vertex and edge descriptors are
 * defined in Types.h (AtomIndex and BondIndex) and the conversion to the
 * underlying boost graph descriptors is implemented in Graph/Bridge.h.
 *
 * The iterators giving easy access to atoms, bonds, adjacents and incidents
 * have to be declared in the header, these are all template specializations of
 * the InnerBasedIterator facade, which are themselves pImpl-ed. Their
 * implementations are in OuterGraphIterators.cpp.
 *
 * The InnerGraph type is the Impl struct of OuterGraph, except it's not a
 * OuterGraph-local type, but a free type that many implementation details of
 * other molassembler components expect to be passed.
 */

#include "molassembler/OuterGraph.h"

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Typenames.h"

#include "molassembler/Cycles.h"
#include "molassembler/Graph/Bridge.h"
#include "molassembler/Modeling/BondDistance.h"

namespace Scine {
namespace molassembler {

static_assert(
  std::is_same<AtomIndex, InnerGraph::Vertex>::value,
  "Atomic index types are mismatched!"
);

OuterGraph::OuterGraph(OuterGraph&& other) noexcept = default;
OuterGraph& OuterGraph::operator = (OuterGraph&& other) noexcept = default;
OuterGraph::OuterGraph(const OuterGraph& other) : _innerPtr(
  std::make_unique<InnerGraph>(*other._innerPtr)
) {}
OuterGraph& OuterGraph::operator = (const OuterGraph& other) {
  *_innerPtr = *other._innerPtr;
  return *this;
}
OuterGraph::~OuterGraph() = default;

OuterGraph::OuterGraph() : _innerPtr(
  std::make_unique<InnerGraph>()
) {}

OuterGraph::OuterGraph(InnerGraph&& inner) : _innerPtr(
  std::make_unique<InnerGraph>(std::move(inner))
) {}

Utils::ElementTypeCollection OuterGraph::elementCollection() const {
  const AtomIndex size = N();

  Utils::ElementTypeCollection elements;
  elements.reserve(size);

  for(AtomIndex i = 0; i < size; ++i) {
    elements.push_back(elementType(i));
  }

  return elements;
}

Utils::ElementType OuterGraph::elementType(const AtomIndex a) const {
  return inner().elementType(a);
}

AtomIndex OuterGraph::N() const {
  return inner().N();
}

unsigned OuterGraph::B() const {
  return inner().B();
}

//! Determine which vertices belong to which side of a bridge edge
std::pair<
  std::vector<AtomIndex>,
  std::vector<AtomIndex>
> OuterGraph::splitAlongBridge(BondIndex bridge) const {
  return inner().splitAlongBridge(
    toInner(bridge, inner())
  );
}

bool OuterGraph::adjacent(const AtomIndex a, const AtomIndex b) const {
  return static_cast<bool>(
    inner().edgeOption(a, b)
  );
}

std::vector<AtomIndex> OuterGraph::atomsOfElement(const Utils::ElementType e) const {
  std::vector<AtomIndex> matches;
  for(AtomIndex i : boost::make_iterator_range(inner().vertices())) {
    if(inner().elementType(i) == e) {
      matches.push_back(i);
    }
  }
  return matches;
}

boost::optional<BondIndex> OuterGraph::bond(const AtomIndex a, const AtomIndex b) const {
  if(auto edgeOption = inner().edgeOption(a, b)) {
    return toOuter(edgeOption.value(), inner());
  }

  return boost::none;
}

Utils::BondOrderCollection OuterGraph::bondOrders() const {
  Utils::BondOrderCollection BOs(inner().N());

  for(const auto edge : boost::make_iterator_range(inner().edges())) {
    BOs.setOrder(
      inner().source(edge),
      inner().target(edge),
      Bond::bondOrderMap.at(
        static_cast<std::underlying_type_t<BondType>>(
          inner().bondType(edge)
        )
      )
    );
  }

  return BOs;
}

BondType OuterGraph::bondType(const BondIndex& edge) const {
  return inner().bondType(toInner(edge, inner()));
}

bool OuterGraph::canRemove(const AtomIndex a) const {
  return inner().canRemove(a);
}

bool OuterGraph::canRemove(const BondIndex& edge) const {
  return inner().canRemove(
    toInner(edge, inner())
  );
}

const Cycles& OuterGraph::cycles() const {
  return inner().cycles();
}

InnerGraph& OuterGraph::inner() {
  return *_innerPtr;
}

const InnerGraph& OuterGraph::inner() const {
  return *_innerPtr;
}

unsigned OuterGraph::degree(const AtomIndex a) const {
  return inner().degree(a);
}

OuterGraph::Range<OuterGraph::AtomIterator> OuterGraph::atoms() const {
  return {
    AtomIterator(inner(), true),
    AtomIterator(inner(), false)
  };
}

OuterGraph::Range<OuterGraph::BondIterator> OuterGraph::bonds() const {
  return {
    BondIterator(inner(), true),
    BondIterator(inner(), false)
  };
}

OuterGraph::Range<OuterGraph::AdjacencyIterator> OuterGraph::adjacents(const AtomIndex a) const {
  return {
    AdjacencyIterator(a, inner(), true),
    AdjacencyIterator(a, inner(), false)
  };
}

OuterGraph::Range<OuterGraph::IncidentEdgesIterator> OuterGraph::bonds(const AtomIndex a) const {
  return {
    IncidentEdgesIterator(a, inner(), true),
    IncidentEdgesIterator(a, inner(), false)
  };
}

} // namespace molassembler
} // namespace Scine
