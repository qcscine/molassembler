/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/IO/SmilesParser.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix/fusion/at.hpp>

#include "Molassembler/IO/SmilesParseData.h"
#include "Molassembler/IO/SmilesMoleculeBuilder.h"
#include "Molassembler/Molecule.h"

#include <iostream>

/* TODO
 * - Allenes
 * - Collect the passed bond data and use it to figure out some Kekule
 *   variation of double bonds in aromatic cycles or at least to set triangle
 *   shapes when elements are aromatic
 *   The spec has the following to say regarding aromaticity: In an aromatic
 *   system, all of the aromatic atoms must be sp 2 hybridized, and the number
 *   of π electrons must meet Huckel’s 4n+2 criterion When parsing a SMILES, a
 *   parser must note the aromatic designation of each atom on input, then when
 *   the parsing is complete, the SMILES software must verify that electrons
 *   can be assigned without violating the valence rules, consistent with the
 *   sp 2 markings, the specified or implied hydrogens, external bonds, and
 *   charges on the atoms.
 *
 *   BUT it also says that "aromatic" bond types are also allowed in
 *   antiaromatic systems such as cyclobutadiene!
 * - Stereostuff
 *   - @/@@
 *     SMILES use the order in which atoms are specified together with @
 *     (anticlockwise looking along the first to the center) to define order.
 *     (Watch out for the special case involving h-count hydrogens instead of
 *     explicit ones, then they are first in the ordering).
 */

namespace Scine {
namespace Molassembler {
namespace IO {
namespace {

std::string abbreviate(std::string x, const unsigned len) {
  assert(len > 3);

  if(x.size() > len) {
    x.resize(len - 3);
    x += "...";
  }

  return x;
}

} // namespace

namespace qi = boost::spirit::qi;

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
} const organic_aliphatic_element;

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
} const organic_aromatic_element;

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
} const aromatic_symbols;

//! All element symbol strings (for use in atom bracket)
struct element_symbols_ : qi::symbols<char, ElementData> {
  element_symbols_() {
    // All symbols from Z = 0 to 110
    for(unsigned i = 1; i < 110; ++i) {
      auto element = Utils::ElementInfo::element(i);
      add(Utils::ElementInfo::symbol(element), ElementData(element));
    }
  }
} const element_symbols;

struct chiral_subset_ : qi::symbols<char, ChiralData> {
  chiral_subset_() {
    add
      ("@", {Shapes::Shape::Tetrahedron, 1})
      ("@@", {Shapes::Shape::Tetrahedron, 2})
      ("@TH1", {Shapes::Shape::Tetrahedron, 1}) // Same as @
      ("@TH2", {Shapes::Shape::Tetrahedron, 2}) // Same as @@
      // ("@AL1", {Shapes::Shape::Tetrahedron, 1}) // Allene
      // ("@AL2", {Shapes::Shape::Tetrahedron, 2}) // Allene
      ("@SP1", {Shapes::Shape::Square, 1})
      ("@SP2", {Shapes::Shape::Square, 2})
      ("@SP3", {Shapes::Shape::Square, 3});
  }
} const chiral_subset;

struct bond_ : qi::symbols<char, BondData> {
  bond_() {
    add
      ("-", {SmilesBondType::Single, boost::none})
      ("=", {SmilesBondType::Double, boost::none})
      ("#", {SmilesBondType::Triple, boost::none})
      ("$", {SmilesBondType::Quadruple, boost::none})
      (":", {SmilesBondType::Aromatic, boost::none})
      ("/", {SmilesBondType::Forward, boost::none})
      ("\\", {SmilesBondType::Backward, boost::none});
  }
} const bond;

} // namespace symbols


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
      (
        lit("@TB")[at_c<0>(_val) = Shapes::Shape::TrigonalBipyramid]
        >> digits<1, 2>()[at_c<1>(_val) = _1]
      )
      | (
        lit("@OH")[at_c<0>(_val) = Shapes::Shape::Octahedron]
        >> digits<1, 2>()[at_c<1>(_val) = _1]
      )
      | symbols::chiral_subset[_val = _1]
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
      lit('[')[at_c<5>(_val) = true] > ( // Opening bracket must be followed by a match of the rest
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
      digits<1>()[at_c<1>(_val) = _1]
      | (lit("%") >> digits<2>()[at_c<1>(_val) = _1])
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
      phoenix::ref(error) = (
        phoenix::val("Expected atom symbol and ']' after atom bracket '[' here: \"")
        + phoenix::bind(
          &abbreviate,
          phoenix::construct<std::string>(qi::_3, qi::_2),
          20
        )
        + phoenix::val("\"\n")
      )
    );

    qi::on_error<qi::fail>(branch,
      phoenix::ref(error) = (
        phoenix::val("Expected branch continuation and ')' after '(' here: \"")
        + phoenix::bind(
          &abbreviate,
          phoenix::construct<std::string>(qi::_3, qi::_2),
          20
        )
        + phoenix::val("\"\n")
      )
    );
  }

  MoleculeBuilder builder;
  std::string error;

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

namespace Experimental {

std::vector<Molecule> parseSmiles(const std::string& smiles) {
  auto iter = std::begin(smiles);
  const auto end = std::end(smiles);

  using Iterator = std::string::const_iterator;
  using Parser = openSMILES<Iterator>;

  Parser parser;
  bool result = qi::parse(iter, end, parser);

  if(result && iter == end) {
    return parser.builder.interpret(smiles);
  }

  throw std::runtime_error("Parsing failure: " + parser.error);
}

Molecule parseSmilesSingleMolecule(const std::string& smiles) {
  auto results = parseSmiles(smiles);

  if(results.size() > 1) {
    throw std::logic_error("Passed smiles strings string contains multiple molecules");
  }

  return results.front();
}

} // namespace Experimental
} // namespace IO
} // namespace Molassembler
} // namespace Scine
