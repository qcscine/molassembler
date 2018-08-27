#include "molassembler/OuterGraph.h"

#include "Delib/ElementTypeCollection.h"

#include "molassembler/Graph/Bridge.h"

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

Delib::ElementTypeCollection OuterGraph::elementCollection() const {
  const AtomIndex size = N();

  Delib::ElementTypeCollection elements;
  elements.reserve(size);

  for(AtomIndex i = 0; i < size; ++i) {
    elements.push_back(
      elementType(i)
    );
  }

  return elements;
}

Delib::ElementType OuterGraph::elementType(const AtomIndex a) const {
  return inner().elementType(a);
}

AtomIndex OuterGraph::N() const {
  return inner().N();
}

unsigned OuterGraph::B() const {
  return inner().B();
}

bool OuterGraph::adjacent(const AtomIndex a, const AtomIndex b) const {
  return static_cast<bool>(
    inner().edgeOption(a, b)
  );
}

boost::optional<BondIndex> OuterGraph::bond(const AtomIndex a, const AtomIndex b) const {
  if(auto edgeOption = inner().edgeOption(a, b)) {
    return toOuter(edgeOption.value(), inner());
  }

  return boost::none;
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
