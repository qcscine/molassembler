#include "molassembler/Layers/Outer.h"

#include "molassembler/Layers/Bridge.h"

namespace molassembler {

static_assert(
  std::is_same<OuterGraph::AtomIndex, InnerGraph::Vertex>::value,
  "Atomic index types are mismatched!"
);

Delib::ElementType OuterGraph::elementType(const OuterGraph::AtomIndex a) const {
  return inner().elementType(a);
}

BondType OuterGraph::bondType(const OuterGraph::BondIndex& a) const {
  return inner().bondType(toInner(a, inner()));
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
