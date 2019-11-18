#include "molassembler/IO/SmilesParser.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix/fusion/at.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include "molassembler/StereopermutatorList.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Molecule.h"
#include "shapes/Shapes.h"
#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Geometry/ElementInfo.h"

#include "temple/Optionals.h"

#include <stack>
#include <iostream>

/* TODO
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
 *     explicit ones, then they are first in the ordering). I think it might be
 *     possible to just figure out whether the ranking of the substituents in
 *     order is an even or odd permutation as proxy for figuring out the
 *     stereopermutation.
 *   - Square planar centers
 *     This is more complicated than tetrahedral owing to the fact that there
 *     are more possible stereopermutations. The specified index corresponds to
 *     a 2D character's shape that is made if you follow the atoms as ordered
 *     in the smiles: One for U, two for 4, three for Z. You can relate these
 *     shapes to the shape vertex indinces for square planar (ccw in order,
 *     so one is just 0123, 3210, or any rotation thereof)
 *   - Trigonal biypramidal centers
 *     The trigonal bipyramidal spec is that each number corresponds to a pair
 *     of apical atoms and an ccw/cw sense of the equatorial atoms viewed
 *     along the axis. The ordering is arbitrarily weird (without explanation).
 *     I think the same applies as before to the shape vertex ordering choice
 *     (4 and 5 are the ordered top and bottom pair forming the axis, and then
 *     the rest can only be rotations of 012 (ccw) or 021 (cw))
 *   - Octahedral centers
 *     Combination of the earlier ideas (square planar shapes combined with a
 *     specified axis and an ordering criterion, have another look at this
 *     when done with the rest)
 */

namespace Scine {
namespace molassembler {
namespace IO {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

struct TrigonalBipyramidSpec {
  std::uint8_t top;
  std::uint8_t bottom;
  bool clockwise;
};

// These are the TBxx numbers as specified
constexpr std::array<TrigonalBipyramidSpec, 20> trigonalBipyramidStereoSpec {{
  {0, 4, false},
  {0, 4, true},
  {0, 3, false},
  {0, 3, true},
  {0, 2, false},
  {0, 2, true},
  {0, 1, false},
  {0, 1, true},
  {1, 4, false}, // Weird order starts below here (this is on purpose)
  {1, 3, false},
  {1, 4, true},
  {1, 3, true},
  {1, 2, false},
  {1, 2, true},
  {2, 4, false},
  {2, 3, false},
  {3, 4, false},
  {3, 4, true},
  {2, 3, true},
  {2, 4, true}
}};

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
  bool atomBracket = false;

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
  enum class StereoMarker {Forward, Backward};

  boost::optional<BondType> type;
  boost::optional<StereoMarker> ezStereo;
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
  (boost::optional<int>, chargeOptional),
  (bool, atomBracket)
)

BOOST_FUSION_ADAPT_STRUCT(
  Scine::molassembler::IO::BondData,
  (boost::optional<Scine::molassembler::BondType>, type),
  (boost::optional<Scine::molassembler::IO::BondData::StereoMarker>, ezStereo)
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
      ("@TH1", {Shapes::Shape::Tetrahedron, 0}) // Same as @
      ("@TH2", {Shapes::Shape::Tetrahedron, 1}) // Same as @@
      ("@AL1", {Shapes::Shape::Tetrahedron, 0}) // Allene
      ("@AL2", {Shapes::Shape::Tetrahedron, 1}) // Allene
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
      ("/", {BondType::Single, BondData::StereoMarker::Forward, boost::none})
      ("\\", {BondType::Single, BondData::StereoMarker::Backward, boost::none});
  }
} bond;

} // namespace symbols

struct MoleculeBuilder {
  static bool isValenceFillElement(Utils::ElementType e) {
    const unsigned Z = Utils::ElementInfo::Z(e);
    if(5 <= Z && Z <= 9) {
      // B, C, N, O, F
      return true;
    }

    if(15 <= Z && Z <= 17) {
      // P, S, Cl
      return true;
    }

    if(Z == 35 || Z == 53) {
      // Br, I
      return true;
    }

    return false;
  }

  static unsigned valenceFillElementImplicitHydrogenCount(
    const int valence,
    Utils::ElementType e
  ) {
    assert(valence >= 0);
    assert(isValenceFillElement(e));

    /* Quoting from the spec:
     *
     * The implicit hydrogen count is determined by summing the bond orders of
     * the bonds connected to the atom. If that sum is equal to a known valence
     * for the element or is greater than any known valence then the implicit
     * hydrogen count is 0. Otherwise the implicit hydrogen count is the
     * difference between that sum and the next highest known valence.
     */

    switch(Utils::ElementInfo::Z(e)) {
      case 5: return std::max(0, 3 - valence); // B
      case 6: return std::max(0, 4 - valence); // C
      case 7: { // N
        return std::min(
          std::max(0, 3 - valence),
          std::max(0, 5 - valence)
        );
      }
      case 8: return std::max(0, 2 - valence); // O
      case 15: { // P
        return std::min(
          std::max(0, 3 - valence),
          std::max(0, 5 - valence)
        );
      }
      case 16: { // S
        return std::min({
          std::max(0, 2 - valence),
          std::max(0, 4 - valence),
          std::max(0, 6 - valence)
        });
      }
      default: return std::max(0, 1 - valence); // F, Cl, Br, I are the remaining cases
    }
  }

  static BondType mutualBondType(
    const boost::optional<BondType>& a,
    const boost::optional<BondType>& b
  ) {
    /* Ensure that the specified bond order matches. Both bond orders
     * are optionals. If one of both is specified, that one's type is used.
     * If neither is specified, use a single bond. If both are specified,
     * their type is used if it matches, otherwise throw.
     */
    if(!a && !b) {
      return BondType::Single;
    }

    if(a && !b) {
      return a.value();
    }

    if(!a && b) {
      return b.value();
    }

    if(a.value() != b.value()) {
      throw std::runtime_error("Mismatched ring closing bond order");
    }

    return a.value();
  }

  // On atom addition
  void addAtom(const AtomData& atom) {
    InnerGraph::Vertex newVertex = graph.addVertex(atom.getElement());

    if(atom.partialElement.Z == 1 && atom.hCount && atom.hCount.value() != 0) {
      throw std::runtime_error("Hydrogen atoms cannot have hydrogen counts!");
    }

    vertexData.push_back(atom);

    if(lastBondData.which() == 0) {
      auto data = boost::get<SimpleLastBondData>(lastBondData);
      if(data == SimpleLastBondData::Unspecified) {
        assert(!vertexStack.empty());
        graph.addEdge(
          vertexStack.top(),
          newVertex,
          BondType::Single
        );
      }
    } else {
      assert(!vertexStack.empty());
      auto data = boost::get<BondData>(lastBondData);
      graph.addEdge(
        vertexStack.top(),
        newVertex,
        data.type.value_or(BondType::Single)
      );

      // Store stereo-marked bonds for later
      if(data.ezStereo) {
        stereoMarkedBonds.emplace_back(
          vertexStack.top(),
          newVertex,
          data.ezStereo.value()
        );
      }
    }

    if(vertexStack.empty()) {
      vertexStack.push(newVertex);
    } else {
      vertexStack.top() = newVertex;
    }

    lastBondData = SimpleLastBondData::Unspecified;
  }

  void addRingClosure(const BondData& bond) {
    assert(bond.ringNumber);
    const unsigned key = bond.ringNumber.value();
    auto findIter = ringClosures.find(key);
    if(findIter == std::end(ringClosures)) {
      // Add the entry to the map for later
      ringClosures.emplace_hint(
        findIter,
        std::piecewise_construct,
        std::make_tuple(key),
        std::make_tuple(vertexStack.top(), bond.type)
      );
    } else {
      // Add the edge to the graph now and remove the map entry
      InnerGraph::Vertex a = findIter->second.first;
      InnerGraph::Vertex b = vertexStack.top();

      if(a == b) {
        throw std::runtime_error("Same-atom ring-closing bond!");
      }

      if(graph.edgeOption(a, b) != boost::none) {
        throw std::runtime_error("Ring closing bond already exists!");
      }

      // Ensure the specified bond types match (this fn can throw)
      const BondType type = mutualBondType(
        findIter->second.second,
        bond.type
      );

      graph.addEdge(a, b, type);

      // Remove the entry from the map
      ringClosures.erase(findIter);
    }
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

  void setBondStereo(
    std::vector<Molecule>& molecules,
    std::vector<unsigned> componentMap,
    std::vector<InnerGraph::Vertex> indexInComponentMap
  ) {
    /* Setting the bond stereo from the forward and backward markers is tricky
     * for several reasons.
     *
     * 1. The "/" and "\" markers indicate the "up" or "down" positioning
     * relative to the carbon atom, which may freely occur before or after the
     * marked atom.
     *
     * 2. All we have from the parsing is the sequence in which the
     * stereomarkers were discovered. During parsing, none of these are
     * directly interpreted. We have to make sure they make sense in the first
     * place (no two markers on one side of the double bond may indicate the
     * same relative positioning).
     *
     * That said, I think there are some patterns we can exploit to structure
     * the mess.
     *
     * F/C=C/F is trans, and so is C(\F)=C/F
     *
     * The only side for which there is freedom of reordering for the
     * stereomarkers is left of the bond. I.e. F/C=(\F)C is invalid. So the only
     * side of the double bond for which we have to figure out the relative
     * ordering is left. We can exploit the double bond that must be present in
     * the graph:
     *
     * F/C=C/F -> AB forward, CD forward and BC double bonded
     * C(\F)=C/F -> AB backward, CD forward and AC double bonded
     *
     * It is legal to mark the second substituent at either side too: We can
     * express this molecule also as:
     *
     * [H]\C(\F)=C/F(\[H]) -> AB backward, BC backward, DE forward, DF backward, B = D
     * C(\F)(/[H])=C/F(\[H]) -> AB backward, AC forward, DE forward, DF backward, A = D
     *
     * Unfortunately, all we have is the ordered sequence of marked bonds, so
     * we have to write a state machine that can deal with that. Index
     * repetitions and changes can be our guides to deciding when we have
     * crossed sides of the bond.
     */
    struct BondStereo {
      boost::optional<InnerGraph::Vertex> left;
      InnerGraph::Vertex right;
      boost::optional<InnerGraph::Vertex> upOfLeft;
      boost::optional<InnerGraph::Vertex> downOfLeft;
      boost::optional<InnerGraph::Vertex> upOfRight;
      boost::optional<InnerGraph::Vertex> downOfRight;

      unsigned findAssignment(
        BondStereopermutator stereopermutator,
        const Molecule& mol,
        const std::vector<InnerGraph::Vertex>& indexInComponentMap
      ) const {
        auto first = mol.stereopermutators().option(stereopermutator.edge().first).value();
        auto second = mol.stereopermutators().option(stereopermutator.edge().second).value();

        if(first.centralIndex() == indexInComponentMap.at(right)) {
          std::swap(first, second);
        }

        auto getSiteIndexLeft = [&](const InnerGraph::Vertex i) -> unsigned {
          return first.getRanking().getSiteIndexOf(indexInComponentMap.at(i));
        };
        auto getSiteIndexRight = [&](const InnerGraph::Vertex i) -> unsigned {
          return second.getRanking().getSiteIndexOf(indexInComponentMap.at(i));
        };

        assert(stereopermutator.numAssignments() == 2);
        for(unsigned i = 0; i < 2; ++i) {
          stereopermutator.assign(i);

          auto upOfLeftSiteIndex = temple::optionals::map(upOfLeft, getSiteIndexLeft);
          auto downOfLeftSiteIndex = temple::optionals::map(downOfLeft, getSiteIndexLeft);
          auto upOfRightSiteIndex = temple::optionals::map(upOfRight, getSiteIndexRight);
          auto downOfRightSiteIndex = temple::optionals::map(downOfRight, getSiteIndexRight);

          if(upOfLeftSiteIndex) {
            if(upOfRightSiteIndex) {
              if(std::fabs(stereopermutator.dihedral(first, *upOfLeftSiteIndex, second, *upOfRightSiteIndex)) > 1e-3) {
                continue;
              }
            }
            if(downOfRightSiteIndex) {
              if(std::fabs(stereopermutator.dihedral(first, *upOfLeftSiteIndex, second, *downOfRightSiteIndex) - M_PI) > 1e-3) {
                continue;
              }
            }
          }
          if(downOfLeftSiteIndex) {
            if(upOfRightSiteIndex) {
              if(std::fabs(stereopermutator.dihedral(first, *downOfLeftSiteIndex, second, *upOfRightSiteIndex) - M_PI) > 1e-3) {
                continue;
              }
            }
            if(downOfRightSiteIndex) {
              if(std::fabs(stereopermutator.dihedral(first, *downOfLeftSiteIndex, second, *downOfRightSiteIndex)) > 1e-3) {
                continue;
              }
            }
          }

          return i;
        }

        throw std::logic_error("Failed to find matching stereopermutation for BondStereo state.");
      }
    };

    auto first = [](const auto& tup) { return std::get<0>(tup); };
    auto second = [](const auto& tup) { return std::get<1>(tup); };
    auto marker = [](const auto& tup) { return std::get<2>(tup); };

    using Iterator = std::vector<StereoMarkedBondTuple>::const_iterator;
    Iterator start = std::cbegin(stereoMarkedBonds);
    const Iterator end = std::cend(stereoMarkedBonds);

    while(start != end) {
      BondStereo state;

      InnerGraph::Vertex A = first(*start);
      InnerGraph::Vertex B = second(*start);

      // We assume that all vertices are in the same component
      Molecule& mol = molecules.at(componentMap.at(A));

      std::vector<Iterator> leftMarkers {start};
      std::vector<Iterator> rightMarkers;

      /* Iff the first two marked bonds have an overlapping atom index, then
       * they are on the same side of the bond. The overlapping bond must
       * then be the left atom.
       */
      Iterator explorer = start + 1;
      if(explorer == end) {
        throw std::runtime_error("Missing right side of stereo-marked double bond");
      }

      // Check for second marker left of bond
      {
        InnerGraph::Vertex X = first(*explorer);

        if(A == X) {
          // Two markers left of bond, C(\F)(/[H]) pattern
          state.left = A;
          leftMarkers.push_back(explorer);
          ++explorer;
        } else if(B == X) {
          // Two markers left of bond, [H]\C(\F) pattern
          state.left = B;
          leftMarkers.push_back(explorer);
          ++explorer;
        }
      }

      // Now we have ensured the explorer is right of the bond
      if(explorer == end) {
        throw std::runtime_error("Missing right side of stereo-marked double bond");
      }

      auto bondTypeOption = [&](const InnerGraph::Vertex a, const InnerGraph::Vertex b) {
        return temple::optionals::map(
          mol.graph().bond(
            indexInComponentMap.at(a),
            indexInComponentMap.at(b)
          ),
          [&](const BondIndex& bond) -> BondType {
            return mol.graph().bondType(bond);
          }
        );
      };

      // Establish the right atom
      {
        rightMarkers.push_back(explorer);
        state.right = first(*explorer);

        // Establish the left atom if it is unknown
        if(!state.left) {
          if(bondTypeOption(A, state.right) == BondType::Double) {
            state.left = A;
          } else if(bondTypeOption(B, state.right) == BondType::Double) {
            state.left = B;
          } else {
            throw std::runtime_error("Right side of marked double bond expected, got unrelated bond");
          }
        }
      }

      /* Check for an additional right marker and place explorer at end of
       * relevant marked bonds
       */
      ++explorer;
      if(explorer != end && first(*explorer) == state.right) {
        rightMarkers.push_back(explorer);
        ++explorer;
      }

      assert(state.left);
      /* Now process the collected markers for directionality */
      for(const Iterator& leftMarker : leftMarkers) {
        /* Four cases:
         * - first of the marker is left and forward: C(/F) -> second is up
         * - first of the marker is left and backward: C(\F) -> second is down
         * - second of the marker is left and forward: F/C -> first is down
         * - second of the marker is left and backward: F\C -> first is up
         */
        bool firstIsLeft = (first(*leftMarker) == state.left.value());
        bool markerIsForward = (marker(*leftMarker) == BondData::StereoMarker::Forward);
        bool up = (firstIsLeft == markerIsForward);
        InnerGraph::Vertex which = (firstIsLeft ? second(*leftMarker) : first(*leftMarker));

        if(up) {
          if(state.upOfLeft) {
            throw std::runtime_error("Both markers left of double bond indicate 'up' directionality");
          }
          state.upOfLeft = which;
        } else {
          if(state.downOfLeft) {
            throw std::runtime_error("Both markers left of double bond indicate 'down' directionality");
          }
          state.downOfLeft = which;
        }
      }
      for(const Iterator& rightMarker : rightMarkers) {
        assert(first(*rightMarker) == state.right);
        if(marker(*rightMarker) == BondData::StereoMarker::Forward) {
          if(state.upOfRight) {
            throw std::runtime_error("Both markers right of double bond indicate 'up' directionality");
          }
          state.upOfRight = second(*rightMarker);
        } else { // Backward
          if(state.downOfRight) {
            throw std::runtime_error("Both markers right of double bond indicate 'down' directionality");
          }
          state.downOfRight = second(*rightMarker);
        }
      }

      /* Add the information to the molecular graph */
      auto molBondOption = mol.graph().bond(
        indexInComponentMap.at(state.left.value()),
        indexInComponentMap.at(state.right)
      );
      assert(molBondOption);

      if(auto stereopermutatorOption = mol.stereopermutators().option(molBondOption.value())) {
        if(stereopermutatorOption->numAssignments() == 2) {
          mol.assignStereopermutator(
            molBondOption.value(),
            state.findAssignment(
              *stereopermutatorOption,
              mol,
              indexInComponentMap
            )
          );
        } else {
          std::cerr << "Warning: Smiles contains stereo markers for non-stereogenic double bond\n";
        }
      } else {
        std::cerr << "Warning: Smiles contains stereo markers for non-stereogenic double bond\n";
      }

      // Advance the start iterator
      start = explorer;
    }
  }

  // Interpret the graph as possibly distinct molecules
  std::vector<Molecule> interpret() {
    if(!ringClosures.empty()) {
      throw std::runtime_error("Unmatched ring closure markers remain!");
    }

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

    /* Copy edges into the separate components */
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

    /* Valence fill organic subset in each precursor */
    assert(vertexData.size() == N);
    for(unsigned i = 0; i < N; ++i) {
      const AtomData& data = vertexData.at(i);
      auto& precursor = precursors.at(componentMap.at(i));
      InnerGraph::Vertex vertexInPrecursor = indexInComponentMap.at(i);

      if(data.hCount) {
        // Fill with specified number of hydrogen atoms
        for(unsigned j = 0; j < data.hCount.value(); ++j) {
          InnerGraph::Vertex newHydrogenVertex = precursor.addVertex(Utils::ElementType::H);
          precursor.addEdge(vertexInPrecursor, newHydrogenVertex, BondType::Single);
        }
      } else if(!data.atomBracket && isValenceFillElement(precursor.elementType(vertexInPrecursor))) {
        // Figure out current valence.
        int currentValence = 0;
        for(
          const InnerGraph::Edge edge :
          boost::make_iterator_range(precursor.edges(vertexInPrecursor))
        ) {
          currentValence += Bond::bondOrderMap.at(
            static_cast<unsigned>(
              precursor.bondType(edge)
            )
          );
        }

        const unsigned fillCount = valenceFillElementImplicitHydrogenCount(
          currentValence,
          precursor.elementType(vertexInPrecursor)
        );

        for(unsigned j = 0; j < fillCount; ++j) {
          InnerGraph::Vertex newHydrogenVertex = precursor.addVertex(Utils::ElementType::H);
          precursor.addEdge(vertexInPrecursor, newHydrogenVertex, BondType::Single);
        }
      }
    }

    /* Convert the graphs to molecules */
    std::vector<Molecule> molecules;
    molecules.reserve(M);
    for(auto&& precursor : precursors) {
      molecules.emplace_back(
        OuterGraph(std::move(precursor))
      );

      /* TODO set shapes here */
    }

    /* Set bond stereo */
    setBondStereo(
      molecules,
      componentMap,
      indexInComponentMap
    );

    return molecules;
  }

  enum class SimpleLastBondData {
    Unbonded,
    Unspecified
  };

  //! State for last stored bond data
  boost::variant<SimpleLastBondData, BondData> lastBondData = SimpleLastBondData::Unbonded;

  //! Possibly disconnected tracking graph
  InnerGraph graph;

  //! State to track the vertex a new vertex is bound to
  std::stack<InnerGraph::Vertex> vertexStack;

  //! Storage for bonds marked with stereo indicators ("/" and "\")
  using StereoMarkedBondTuple = std::tuple<InnerGraph::Vertex, InnerGraph::Vertex, BondData::StereoMarker>;
  std::vector<StereoMarkedBondTuple> stereoMarkedBonds;

  //! Storage for ring closure bond indicators
  std::unordered_map<
    unsigned,
    std::pair<InnerGraph::Vertex, boost::optional<BondType>>
  > ringClosures;

  //! AtomData for each created vertex
  std::vector<AtomData> vertexData;
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

  throw std::runtime_error("Failed to parse SMILES");
}

} // namespace IO
} // namespace molassembler
} // namespace Scine
