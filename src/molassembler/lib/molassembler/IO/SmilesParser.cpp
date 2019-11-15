#include "molassembler/IO/SmilesParser.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix/fusion/at.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include "molassembler/Molecule.h"
#include "molassembler/Graph/InnerGraph.h"
#include "shapes/Shapes.h"
#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Geometry/ElementInfo.h"

#include <stack>
#include <iostream>

/* TODO
 * - Piece together the molecules by adding the cycle closing bonds when paired
 * - Collect the passed atom data and use to apply valence rules to fill up
 *   with hydrogen atoms
 * - Collect the passed bond data and use it to figure out some Kekule
 *   variation of double bonds in aromatic cycles or at least to set triangle
 *   shapes
 * - Stereostuff
 *   - Double bond stereo
 *   - @/@@
 */

namespace Scine {
namespace molassembler {
namespace IO {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

struct ElementData {
  unsigned Z = 0;
  bool aromatic = false;

  ElementData() = default;
  ElementData(Utils::ElementType e) : Z(Utils::ElementInfo::Z(e)) {}

  static ElementData aromaticElement(Utils::ElementType e) {
    ElementData d(e);
    d.aromatic = true;
    return d;
  }
};

struct ChiralData {
  Shapes::Shape shape;
  unsigned chiralIndex;
};

struct AtomData {
  unsigned A = 0;
  ElementData partialElement;
  boost::optional<ChiralData> chiralOptional;
  boost::optional<unsigned> hCount;
  boost::optional<int> chargeOptional;

  Utils::ElementType getElement() const {
    if(partialElement.Z == 0) {
      return Utils::ElementType::none;
    }

    if(A == 0) {
      return Utils::ElementInfo::element(partialElement.Z);
    }

    return Utils::ElementInfo::isotope(partialElement.Z, A);
  }
};

struct BondData {
  BondType type = BondType::Single;
  boost::optional<unsigned> ezStereo;
  boost::optional<unsigned> ringNumber;
};

} // namespace IO
} // namespace molassembler
} // namespace Scine


BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::ElementData,
  (unsigned, Z),
  (bool, aromatic)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::ChiralData,
  (Scine::Shapes::Shape, shape),
  (unsigned, chiralIndex)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::AtomData,
  (unsigned, A),
  (Scine::molassembler::IO::ElementData, partialElement),
  (boost::optional<Scine::molassembler::IO::ChiralData>, chiralOptional),
  (boost::optional<unsigned>, hCount),
  (boost::optional<int>, chargeOptional)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::BondData,
  (Scine::molassembler::BondType, type),
  (boost::optional<unsigned>, ezStereo)
  (boost::optional<unsigned>, ringNumber)
)

namespace Scine {
namespace molassembler {
namespace IO {
namespace symbols {

struct organic_aliphatic_element_ : qi::symbols<char, ElementData> {
  organic_aliphatic_element_() {
    add
      ("B",  ElementData(Utils::ElementType::B))
      ("C",  ElementData(Utils::ElementType::C))
      ("N",  ElementData(Utils::ElementType::N))
      ("O",  ElementData(Utils::ElementType::O))
      ("S",  ElementData(Utils::ElementType::S))
      ("P",  ElementData(Utils::ElementType::P))
      ("F",  ElementData(Utils::ElementType::F))
      ("Cl", ElementData(Utils::ElementType::Cl))
      ("Br", ElementData(Utils::ElementType::Br))
      ("I",  ElementData(Utils::ElementType::I));
  }
} organic_aliphatic_element;

//! Organic aromatic elements (for alternative in atom)
struct organic_aromatic_element_ : qi::symbols<char, ElementData> {
  organic_aromatic_element_() {
    add
      ("b", ElementData::aromaticElement(Utils::ElementType::B))
      ("c", ElementData::aromaticElement(Utils::ElementType::C))
      ("n", ElementData::aromaticElement(Utils::ElementType::N))
      ("o", ElementData::aromaticElement(Utils::ElementType::O))
      ("s", ElementData::aromaticElement(Utils::ElementType::S))
      ("p", ElementData::aromaticElement(Utils::ElementType::P));
  }
} organic_aromatic_element;

//! Aromatic symbols (for use in atom bracket)
struct aromatic_symbols_ : qi::symbols<char, ElementData> {
  aromatic_symbols_() {
    add
      ("b",  ElementData::aromaticElement(Utils::ElementType::B))
      ("c",  ElementData::aromaticElement(Utils::ElementType::C))
      ("n",  ElementData::aromaticElement(Utils::ElementType::N))
      ("o",  ElementData::aromaticElement(Utils::ElementType::O))
      ("s",  ElementData::aromaticElement(Utils::ElementType::S))
      ("p",  ElementData::aromaticElement(Utils::ElementType::P))
      ("se", ElementData::aromaticElement(Utils::ElementType::Se))
      ("as", ElementData::aromaticElement(Utils::ElementType::As));
  }
} aromatic_symbols;

//! All element symbol strings (for use in atom bracket)
struct element_symbols_ : qi::symbols<char, ElementData> {
  element_symbols_() {
    // All symbols from Z = 0 to 110
    for(unsigned i = 1; i < 110; ++i) {
      auto element = Utils::ElementInfo::element(i);
      add(Utils::ElementInfo::symbol(element), ElementData(element));
    }
  }
} element_symbols;

struct chiral_subset_ : qi::symbols<char, ChiralData> {
  chiral_subset_() {
    add
      ("@", {Shapes::Shape::Tetrahedron, 0})
      ("@@", {Shapes::Shape::Tetrahedron, 1})
      ("@TH1", {Shapes::Shape::Tetrahedron, 0})
      ("@TH2", {Shapes::Shape::Tetrahedron, 1})
      ("@AL1", {Shapes::Shape::Tetrahedron, 0})
      ("@AL2", {Shapes::Shape::Tetrahedron, 1})
      ("@SP1", {Shapes::Shape::Square, 0})
      ("@SP2", {Shapes::Shape::Square, 1})
      ("@SP3", {Shapes::Shape::Square, 2});
  }
} chiral_subset;

struct bond_ : qi::symbols<char, BondData> {
  bond_() {
    add
      ("-", {BondType::Single, boost::none, boost::none})
      ("=", {BondType::Double, boost::none, boost::none})
      ("#", {BondType::Triple, boost::none, boost::none})
      ("$", {BondType::Quadruple, boost::none, boost::none})
      (":", {BondType::Single, boost::none, boost::none})  // TODO actually aromatic
      ("/", {BondType::Single, 0, boost::none})
      ("\\", {BondType::Single, 1, boost::none});
  }
} bond;

} // namespace symbols

struct MoleculeBuilder {
  // On atom addition
  void addAtom(const AtomData& atom) {
    InnerGraph::Vertex newVertex = graph.addVertex(atom.getElement());

    if(lastBondData.which() == 0) {
      auto data = boost::get<SimpleLastBondData>(lastBondData);
      if(data == SimpleLastBondData::Single) {
        assert(!vertexStack.empty());
        graph.addEdge(
          newVertex,
          vertexStack.top(),
          BondType::Single
        );
      }
    } else {
      assert(!vertexStack.empty());
      auto data = boost::get<BondData>(lastBondData);
      graph.addEdge(
        newVertex,
        vertexStack.top(),
        data.type
      );
    }

    if(vertexStack.empty()) {
      vertexStack.push(newVertex);
    } else {
      vertexStack.top() = newVertex;
    }

    lastBondData = SimpleLastBondData::Single;
  }

  void addRingClosure(const BondData& bond) {
    ringClosures.emplace_back(
      vertexStack.top(),
      bond
    );
  }

  // Trigger on branch open
  void branchOpen() {
    vertexStack.push(vertexStack.top());
  }

  void branchClose() {
    assert(!vertexStack.empty());
    vertexStack.pop();
  }

  // Trigger on dot bond parse
  void setNextAtomUnbonded() {
    lastBondData = SimpleLastBondData::Unbonded;
  }

  // Triggered on non-default bond information after atom addition
  void setNextAtomBondInformation(const BondData& bond) {
    lastBondData = bond;
  }

  // Interpret the graph as possibly distinct molecules
  std::vector<Molecule> interpret() {
    std::vector<unsigned> componentMap;
    const unsigned M = graph.connectedComponents(componentMap);

    std::vector<InnerGraph> precursors;
    precursors.resize(M);

    const unsigned N = graph.N();

    std::vector<InnerGraph::Vertex> indexInComponentMap(N);
    // Copy vertices
    for(unsigned i = 0; i < N; ++i) {
      auto& precursor = precursors.at(componentMap.at(i));
      InnerGraph::Vertex newIndex = precursor.addVertex(graph.elementType(i));
      indexInComponentMap.at(i) = newIndex;
    }

    // Copy edges
    for(const InnerGraph::Edge& edge : boost::make_iterator_range(graph.edges())) {
      const InnerGraph::Vertex source = graph.source(edge);
      const InnerGraph::Vertex target = graph.target(edge);

      // Both vertices must be in the same component
      auto& precursor = precursors.at(componentMap.at(source));

      precursor.addEdge(
        indexInComponentMap.at(source),
        indexInComponentMap.at(target),
        graph.bondType(edge)
      );
    }

    std::vector<Molecule> molecules;
    molecules.reserve(M);
    for(auto&& precursor : precursors) {
      molecules.emplace_back(
        OuterGraph(std::move(precursor))
      );
    }
    return molecules;
  }

  enum class SimpleLastBondData {
    Unbonded,
    Single
  };

  boost::variant<SimpleLastBondData, BondData> lastBondData = SimpleLastBondData::Unbonded;
  InnerGraph graph;
  std::stack<InnerGraph::Vertex> vertexStack;
  std::vector<
    std::pair<AtomIndex, BondData>
  > ringClosures;
};

template<typename Iterator>
struct openSMILES : qi::grammar<Iterator> {
  template<typename F, typename ... Args>
  auto bind(F&& f, Args ... args) {
    return boost::phoenix::bind(std::forward<F>(f), std::forward<Args>(args) ...);
  }

  template<unsigned count>
  auto digits() {
    return qi::uint_parser<unsigned, 10, count, count>();
  }

  template<unsigned lower, unsigned upper>
  auto digits() {
    static_assert(lower <= upper, "Reversed range!");
    return qi::uint_parser<unsigned, 10, lower, upper>();
  }

  openSMILES() : openSMILES::base_type(start) {
    namespace phoenix = boost::phoenix;
    using phoenix::at_c;
    using qi::_1;
    using qi::_val;
    using qi::lit;
    using qi::uint_;
    using qi::eps;

    // builder fns
    auto addAtom = bind([this](const AtomData& x) { builder.addAtom(x); }, qi::_1);
    auto addRingClosure = bind([this](const BondData& x) { builder.addRingClosure(x); }, qi::_1);
    auto branchOpen = bind([this]() { builder.branchOpen(); });
    auto branchClose = bind([this]() { builder.branchClose(); });
    auto setNextAtomUnbonded = bind([this]() { builder.setNextAtomUnbonded(); });
    auto setNextAtomBondInformation = bind([this](const BondData& x) { builder.setNextAtomBondInformation(x); }, qi::_1);

    // chiral ::= lots of cases (see chiral_subset and the @TB(num) @OH(num) cases here
    chiral = (
      symbols::chiral_subset[_val = _1]
      | (
        lit("@TB")[at_c<0>(_val) = Shapes::Shape::TrigonalBipyramid]
        >> digits<1, 2>()[at_c<1>(_val) = _1]
      )
      | (
        lit("@OH")[at_c<0>(_val) = Shapes::Shape::Octahedron]
        >> digits<1, 2>()[at_c<1>(_val) = _1]
      )
    );

    /* hcount is defined as
     *
     * hcount ::= H | H digit
     *
     * But we could possibly want more hydrogens than 9 for inorganic cases,
     * so we accept up to two digits
     */
    hcount = eps[_val = 0] >> lit('H')[_val = 1] >> -(digits<1, 2>()[_val = _1]);

    // charge ::= `-` num | `+` num | `--` | `++`
    charge = (
      (lit('-')[_val = -1] >> -(digits<1, 2>()[_val = -_1]))
      | (lit('+')[_val = +1] >> -(digits<1, 2>()[_val = +_1]))
      | lit("--")[_val = -2]
      | lit("++")[_val = +2]
    );

    // atom_class ::= `:` num
    atom_class = lit(':') >> uint_;

    // bracket_atom ::= `[` isotope? symbol chiral? hcount? charge? class? `]`
    bracket_atom = (
      lit('[') > ( // Opening bracket must be followed by a match of the rest
        -digits<1, 3>()[at_c<0>(_val) = _1]
        >> ( // symbol
          symbols::aromatic_symbols[at_c<1>(_val) = _1]
          | symbols::element_symbols[at_c<1>(_val) = _1]
          | lit('*')
        )
        >> -chiral[at_c<2>(_val) = _1]
        >> -hcount[at_c<3>(_val) = _1]
        >> -charge[at_c<4>(_val) = _1]
        >> -atom_class
        >> lit(']')
      )
    );

    atom = (
      bracket_atom[_val = _1]
      | symbols::organic_aliphatic_element[at_c<1>(_val) = _1]
      | symbols::organic_aromatic_element[at_c<1>(_val) = _1]
      | lit('*')
    );

    /* Now for the bond stuff */
    // Bond order and stereo
    bond = symbols::bond[_val = _1];
    // Bond info and a ring closure
    ringbond = -bond[_val = _1] >> (
      digits<1>()[at_c<2>(_val) = _1]
      | (lit("%") >> digits<2>()[at_c<2>(_val) = _1])
    );
    // Branching atom: atom and any ring bonds with their entire branches
    branched_atom = atom[addAtom] >> *(ringbond[addRingClosure]) >> *branch;
    /* Branch is defined in the spec as:
     *
     * branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
     *
     * but we do two transformations:
     * 1. '(' (bond | dot)? chain ')' (? is optional)
     * 2. '(' (bond | dot | eps) chain ')' (eps == epsilon is empty string)
     *
     * We need to trigger some addition of the bond order to the builder even
     * if no information is matched (empty string last, representing the
     * default bond order single), hence we need epsilon to place that semantic
     * action
     *
     */
    branch = (
      lit('(')[branchOpen] > ( // Opening bracket must be followed by match of the rest
        (bond[setNextAtomBondInformation] | dot[setNextAtomUnbonded] | eps)
        >> chain
        >> lit(')')[branchClose]
      )
    );
    /* Chain is defined in the spec as:
     * chain ::= (
     *   branched_atom
     *   | chain branched_atom
     *   | chain bond branched_atom
     *   | chain dot branched_atom
     * )
     *
     * But we want to add bonds and atoms eagerly before the point of recursion,
     * not at the end when the recursive descent collapses.
     *
     * So, equivalently (I think), reversing order:
     * chain ::= (
     *   branched_atom
     *   | branched_atom chain
     *   | branched_atom bond chain
     *   | branched_atom dot chain
     * )
     *
     * and also equivalently (I think):
     * chain ::= branched_atom (eps | chain | bond chain | dot chain)
     * chain ::= branched_atom (chain | bond chain | dot chain)?
     * chain ::= branched_atom ((bond | dot | eps) chain)?
     * chain ::= branched_atom (bond | dot | eps) chain?
     *
     * We need the eps as before in branch to trigger some bond order
     * information push semantic action.
     */
    chain = (
      branched_atom
      >> (bond[setNextAtomBondInformation] | dot[setNextAtomUnbonded] | eps)
      >> -chain
    );
    // Just a dot. Molecule separator. Will have its own little action later.
    dot = lit('.');

    start = chain;

    /* Error handling */
    qi::on_error<qi::fail>(bracket_atom,
      std::cout << phoenix::val("Expected atom symbol and ']' after opening of atom bracket '[' here: \"")
        << phoenix::construct<std::string>(qi::_3, qi::_2)
        << phoenix::val("\"\n")
    );

    qi::on_error<qi::fail>(branch,
      std::cout << phoenix::val("Expected branch continuation and ')' after branch opening '(' here: \"")
        << phoenix::construct<std::string>(qi::_3, qi::_2)
        << phoenix::val("\"\n")
    );
  }

  MoleculeBuilder builder;

  // Everything needed for atom
  qi::rule<Iterator, ChiralData()> chiral;
  qi::rule<Iterator, unsigned()> hcount;
  qi::rule<Iterator, int()> charge;
  qi::rule<Iterator> atom_class;
  qi::rule<Iterator, AtomData()> bracket_atom;
  qi::rule<Iterator, AtomData()> atom;

  // Everything needed for bonds
  qi::rule<Iterator, BondData()> bond;
  qi::rule<Iterator, BondData()> ringbond;
  qi::rule<Iterator> branched_atom;
  qi::rule<Iterator> branch;
  qi::rule<Iterator> chain;
  qi::rule<Iterator> dot;

  qi::rule<Iterator> start;
};

std::vector<Molecule> parseSmiles(const std::string& smiles) {
  auto iter = std::begin(smiles);
  const auto end = std::end(smiles);

  using Iterator = std::string::const_iterator;
  using Parser = openSMILES<Iterator>;

  Parser parser;
  bool result = qi::parse(iter, end, parser);

  if(result && iter == end) {
    return parser.builder.interpret();
  }

  throw std::logic_error("Failed to parse SMILES");
}

} // namespace IO
} // namespace molassembler
} // namespace Scine
