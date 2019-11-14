#include "molassembler/IO/SmilesParser.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/phoenix/fusion/at.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include "molassembler/Molecule.h"
#include "shapes/Shapes.h"
#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Geometry/ElementInfo.h"

/* TODO
 * - uint_ to digit changes (otherwise ringbond doesn't make much sense)
 */

namespace Scine {
namespace molassembler {
namespace IO {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

struct ChiralData {
  Shapes::Shape shape;
  unsigned chiralIndex;
};

struct AtomData {
  unsigned Z = 0;
  unsigned A = 0;
  boost::optional<ChiralData> chiralOptional;
  boost::optional<unsigned> hCount;
  boost::optional<int> chargeOptional;

  AtomData() = default;
  AtomData(Utils::ElementType e) : Z(Utils::ElementInfo::Z(e)), A(Utils::ElementInfo::A(e)) {}
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
  Scine::molassembler::IO::ChiralData,
  (Scine::Shapes::Shape, shape),
  (unsigned, chiralIndex)
);

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::AtomData,
  (unsigned, Z),
  (unsigned, A),
  (boost::optional<Scine::molassembler::IO::ChiralData>, chiralOptional),
  (boost::optional<unsigned>, hCount),
  (boost::optional<int>, chargeOptional)
);

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::BondData,
  (Scine::molassembler::BondType, type),
  (boost::optional<unsigned>, ezStereo)
  (boost::optional<unsigned>, ringNumber)
);

namespace Scine {
namespace molassembler {
namespace IO {
namespace symbols {

struct organic_aliphatic_element_ : qi::symbols<char, AtomData> {
  organic_aliphatic_element_() {
    add
      ("B",  AtomData(Utils::ElementType::B))
      ("C",  AtomData(Utils::ElementType::C))
      ("N",  AtomData(Utils::ElementType::N))
      ("O",  AtomData(Utils::ElementType::O))
      ("S",  AtomData(Utils::ElementType::S))
      ("P",  AtomData(Utils::ElementType::P))
      ("F",  AtomData(Utils::ElementType::F))
      ("Cl", AtomData(Utils::ElementType::Cl))
      ("Br", AtomData(Utils::ElementType::Br))
      ("I",  AtomData(Utils::ElementType::I));
  }
} organic_aliphatic_element;

//! Organic aromatic elements (for alternative in atom)
struct organic_aromatic_element_ : qi::symbols<char, AtomData> {
  // TODO AtomData needs to be able to hold that these are aromatic for use later
  organic_aromatic_element_() {
    add
      ("b", AtomData(Utils::ElementType::B))
      ("c", AtomData(Utils::ElementType::C))
      ("n", AtomData(Utils::ElementType::N))
      ("o", AtomData(Utils::ElementType::O))
      ("s", AtomData(Utils::ElementType::S))
      ("p", AtomData(Utils::ElementType::P));
  }
} organic_aromatic_element;

//! Aromatic symbols (for use in atom bracket)
struct aromatic_symbols_ : qi::symbols<char, unsigned> {
  aromatic_symbols_() {
    add
      ("b",  Utils::ElementInfo::Z(Utils::ElementType::B))
      ("c",  Utils::ElementInfo::Z(Utils::ElementType::C))
      ("n",  Utils::ElementInfo::Z(Utils::ElementType::N))
      ("o",  Utils::ElementInfo::Z(Utils::ElementType::O))
      ("s",  Utils::ElementInfo::Z(Utils::ElementType::S))
      ("p",  Utils::ElementInfo::Z(Utils::ElementType::P))
      ("se", Utils::ElementInfo::Z(Utils::ElementType::Se))
      ("as", Utils::ElementInfo::Z(Utils::ElementType::As));
  }
} aromatic_symbols;

//! All element symbol strings (for use in atom bracket)
struct element_symbols_ : qi::symbols<char, AtomData> {
  element_symbols_() {
    // All symbols from Z = 0 to 110
    for(unsigned i = 1; i < 110; ++i) {
      auto element = Utils::ElementInfo::element(i);
      add(Utils::ElementInfo::symbol(element), AtomData(element));
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

struct MoleculeData {};

template<typename Iterator>
struct openSMILES : qi::grammar<Iterator, MoleculeData()> {
  openSMILES() : openSMILES::base_type(start) {
    namespace phoenix = boost::phoenix;
    using phoenix::at_c;
    using qi::_1;
    using qi::_val;
    using qi::lit;
    using qi::uint_;
    using qi::eps;

    // chiral ::= lots of cases (see chiral_subset and the @TB(num) @OH(num) cases here
    chiral = (
      symbols::chiral_subset[_val = _1]
      | (
        lit("@TB")[at_c<0>(_val) = Shapes::Shape::TrigonalBipyramid]
        >> uint_[at_c<1>(_val) = _1]
      )
      | (
        lit("@OH")[at_c<0>(_val) = Shapes::Shape::Octahedron]
        >> uint_[at_c<1>(_val) = _1]
      )
    );

    // hcount ::= H | H digit
    hcount = eps[_val = 0] >> lit('H')[_val = 1] >> -(uint_[_val = _1]);

    // charge ::= `-` num | `+` num | `--` | `++`
    charge = (
      (lit('-')[_val = -1] >> -(uint_[_val = -_1]))
      | (lit('+')[_val = +1] >> -(uint_[_val = +_1]))
      | lit("--")[_val = -2]
      | lit("++")[_val = +2]
    );

    atom_class = lit(':') >> uint_;

    bracket_atom = (
      lit('[')
      >> -uint_[at_c<1>(_val) = _1]
      >> ( // symbol
        symbols::aromatic_symbols[at_c<0>(_val) = _1]
        | symbols::element_symbols[at_c<0>(_val) = _1]
        | lit('*')
      )
      >> -chiral[at_c<2>(_val) = _1]
      >> -hcount[at_c<3>(_val) = _1]
      >> -charge[at_c<4>(_val) = _1]
      >> -atom_class
      >> lit(']')
    );

    atom = (
      bracket_atom[_val = _1]
      | symbols::organic_aliphatic_element[_val = _1]
      | symbols::organic_aromatic_element[_val = _1]
      | lit('*')
    );

    /* Now for the bond stuff
     *
     * TODO YIKES
     */
    // Bond order and stereo
    bond = symbols::bond[_val = _1];
    // Bond info and a ring closure
    ringbond = (
      (-bond[_val = _1] >> uint_[at_c<3>(_val) = _1])
      | (-bond[_val = _1] >> lit("%") >> uint_[at_c<3>(_val) = _1])
    );
    // Branching atom: atom and any ring bonds with their entire branches
    branched_atom = atom >> *ringbond >> *branch;
    branch = (
      (lit('(') >> chain >> lit(')'))
      | (lit('(') >> bond >> chain >> lit(')'))
      | (lit('(') >> dot >> chain >> lit(')'))
    );
    chain = (
      branched_atom
      | (chain >> branched_atom)
      | (chain >> bond >> branched_atom)
      | (chain >> dot >> branched_atom)
    );
    // Just a dot. Molecule separator
    dot = lit('.');

    start = chain;
  }

  // Everything needed for atom
  qi::rule<Iterator, ChiralData()> chiral;
  qi::rule<Iterator, unsigned()> hcount;
  qi::rule<Iterator, int()> charge;
  qi::rule<Iterator> atom_class;
  qi::rule<Iterator, AtomData()> bracket_atom;
  qi::rule<Iterator, AtomData()> atom;

  // Everything needed for bonds
  qi::rule<Iterator, BondData()> bond;
  qi::rule<Iterator> ringbond;
  qi::rule<Iterator> branched_atom;
  qi::rule<Iterator> branch;
  qi::rule<Iterator> chain;
  qi::rule<Iterator> dot;

  qi::rule<Iterator, MoleculeData()> start;
};

Molecule parseSMILES(const std::string& smiles) {
}

} // namespace IO
} // namespace molassembler
} // namespace Scine
