/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Serialization.h"

#include "shapes/Data.h"
#include "boost/range/adaptor/map.hpp"
#include "Utils/Typenames.h"
#include "json/json.hpp"

#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/IO/Base64.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Version.h"

#include "temple/Functional.h"

namespace nlohmann {

template<>
struct adl_serializer<Scine::molassembler::BondStereopermutator::Alignment> {
  using Type = Scine::molassembler::BondStereopermutator::Alignment;
  using Underlying = std::underlying_type_t<Type>;

  static void to_json(json& j, const Type& value) {
    j = static_cast<Underlying>(value);
  }

  static void from_json(const json& j, Type& value) {
    value = static_cast<Type>(
      j.get<Underlying>()
    );
  }
};

template<>
struct adl_serializer<Scine::molassembler::AtomEnvironmentComponents> {
  using Type = Scine::molassembler::AtomEnvironmentComponents;
  using Underlying = std::underlying_type_t<Type>;

  static void to_json(json& j, const Type& value) {
    j = static_cast<Underlying>(value);
  }

  static void from_json(const json& j, Type& value) {
    value = static_cast<Type>(
      j.get<Underlying>()
    );
  }
};

template<>
struct adl_serializer<Scine::Utils::ElementType> {
  /* Static casting element types around here is fine. It's not particularly
   * human readable for monoisotopic elements or isotopes, but certainly
   * doesn't violate correctness in any fashion.
   */
  using Type = Scine::Utils::ElementType;
  using Underlying = std::underlying_type<Scine::Utils::ElementType>::type;

  static void to_json(json& j, const Type& value) {
    j = static_cast<Underlying>(value);
  }

  static void from_json(const json& j, Type& value) {
    value = static_cast<Scine::Utils::ElementType>(
      j.get<Underlying>()
    );
  }
};

template<>
struct adl_serializer<Scine::Utils::ElementTypeCollection> {
  using Type = Scine::Utils::ElementTypeCollection;

  static void to_json(json& j, const Type& value) {
    j = json::array();
    for(const Scine::Utils::ElementType elementType : value) {
      j.push_back(elementType);
    }
  }

  static void from_json(const json& j, Type& value) {
    for(const auto& elementJSON : j) {
      value.push_back(
        elementJSON.get<Scine::Utils::ElementType>()
      );
    }
  }
};

template<>
struct adl_serializer<Scine::molassembler::LinkInformation> {
  using Type = Scine::molassembler::LinkInformation;

  static void to_json(json& j, const Type& link) {
    j["p"] = json::array();
    j["p"].push_back(link.indexPair.first);
    j["p"].push_back(link.indexPair.second);

    j["seq"] = link.cycleSequence;
  }

  static void from_json(const json& j, Type& link) {
    link.indexPair = {
      j["p"].at(0),
      j["p"].at(1)
    };

    link.cycleSequence.reserve(j["seq"].size());
    for(const auto& sequenceElementJSON : j["seq"]) {
      link.cycleSequence.push_back(sequenceElementJSON);
    }

    /* Compatibility fix: old formats repeated the central index in the cycle
     * sequence at the back.
     */
    if(link.cycleSequence.front() == link.cycleSequence.back()) {
      link.cycleSequence.pop_back();
    }
  }
};

template<>
struct adl_serializer<Scine::molassembler::RankingInformation> {
  using Type = Scine::molassembler::RankingInformation;

  static void to_json(json& j, const Type& ranking) {
    j["s"] = ranking.substituentRanking;
    // Omit links member if the list is empty (common)
    if(!ranking.links.empty()) {
      j["lnk"] = ranking.links;
    }
    j["l"] = ranking.sites;
    j["lr"] = ranking.siteRanking;
  }

  static void from_json(const json& j, Type& ranking) {
    ranking.substituentRanking.reserve(j["s"].size());

    for(const auto& listJSON : j["s"]) {
      std::vector<Scine::molassembler::AtomIndex> subGroup;
      for(const auto& listElementJSON : listJSON) {
        subGroup.push_back(listElementJSON);
      }

      ranking.substituentRanking.push_back(std::move(subGroup));
    }

    if(j.count("lnk") > 0) {
      ranking.links.reserve(j["lnk"].size());

      for(const auto& listJSON : j["lnk"]) {
        ranking.links.push_back(listJSON);
      }
    }

    ranking.sites.reserve(j["l"].size());
    for(const auto& listJSON : j["l"]) {
      std::vector<Scine::molassembler::AtomIndex> siteConstitutingAtoms;
      for(const auto& listElementJSON : listJSON) {
        siteConstitutingAtoms.push_back(listElementJSON);
      }
      ranking.sites.push_back(siteConstitutingAtoms);
    }

    ranking.siteRanking.reserve(j["lr"].size());
    for(const auto& listJSON : j["lr"]) {
      std::vector<unsigned> equalSiteIndices;
      for(const auto& listElementJSON : listJSON) {
        equalSiteIndices.push_back(listElementJSON);
      }
      ranking.siteRanking.push_back(equalSiteIndices);
    }
  }
};

template<>
struct adl_serializer<Scine::molassembler::Graph> {
  using Type = Scine::molassembler::Graph;

  static void to_json(json& j, const Type& graph) {
    const Scine::molassembler::PrivateGraph& inner = graph.inner();

    j["Z"] = json::array();
    auto& elements = j["Z"];

    for(const auto vertexIndex : inner.vertices()) {
      elements.push_back(
        inner.elementType(vertexIndex)
      );
    }


    j["E"] = json::array();
    auto& edges = j["E"];

    for(
      const Scine::molassembler::PrivateGraph::Edge& edgeDescriptor :
      inner.edges()
    ) {
      json e = json::array();
      e.push_back(
        static_cast<int>(
          inner.source(edgeDescriptor)
        )
      );
      e.push_back(
        static_cast<int>(
          inner.target(edgeDescriptor)
        )
      );
      e.push_back(
        static_cast<int>(
          inner.bondType(edgeDescriptor)
        )
      );

      edges.push_back(std::move(e));
    }
  }

  static void from_json(const json& j, Type& graph) {
    const unsigned N = j["Z"].size();

    Scine::molassembler::PrivateGraph inner (N);

    for(unsigned i = 0; i < N; ++i) {
      inner.elementType(i) = j["Z"].at(i);
    }

    for(const auto& edgeJSON : j["E"]) {
      inner.addEdge(
        edgeJSON.at(0),
        edgeJSON.at(1),
        static_cast<Scine::molassembler::BondType>(
          edgeJSON.at(2)
        )
      );
    }

    graph = Type {std::move(inner)};
  }
};

} // namespace nlohmann

namespace Scine {

namespace molassembler {

// Top level keys into molecule json object, for easier reading
constexpr const char* atomStereopermutatorKey = "a";
constexpr const char* bondStereopermutatorKey = "b";
constexpr const char* graphKey = "g";
constexpr const char* versionKey = "v";
constexpr const char* canonicalKey = "c";

namespace detail {

void sortAtomStereopermutatorsByCentralIndex(nlohmann::json& m) {
  temple::inplace::sort(
    m.at(atomStereopermutatorKey),
    [&](const nlohmann::json& lhs, const nlohmann::json& rhs) -> bool {
      assert(lhs.count("c") > 0 && rhs.count("c") > 0);
      return (
        lhs["c"].template get<unsigned>()
        < rhs["c"].template get<unsigned>()
      );
    }
  );
}

void standardizeAtomStereopermutatorRepresentation(nlohmann::json& permutator) {
  using json = nlohmann::json;

  assert(permutator.count("r") > 0);
  json& ranking = permutator.at("r");
  assert(ranking.is_object());

  // Sort sub-lists of sorted substituents
  assert(ranking.count("s") > 0);
  json& sortedSubstituents = ranking.at("s");
  assert(sortedSubstituents.is_array());
  for(json& subList : sortedSubstituents) {
    assert(subList.is_array());
    temple::inplace::sort(subList);
  }

  // Sort ligands' sub-lists
  assert(ranking.count("l") > 0);
  json& ligands = ranking.at("l");
  assert(ligands.is_array());
  for(json& subList : ligands) {
    assert(subList.is_array());
    temple::inplace::sort(subList);
  }

  /* Sort ligands lists lexicographically (and their ranking
   * simultaneously), since they are index-connected
   */
  assert(ranking.count("lr") > 0);
  json& ligandsRanking = ranking.at("lr");
  assert(ligandsRanking.is_array());

  json unsortedLigands = ligands;
  temple::inplace::sort(ligands);

  auto newLigandIndex = [&](unsigned oldLigandIndex) -> unsigned {
    const auto& ligandIndices = unsortedLigands.at(oldLigandIndex);
    auto findIter = temple::find(ligands, ligandIndices);
    assert(findIter != std::end(ligands));
    return findIter - std::begin(ligands);
  };

  // Changing ligands means ligandsRanking has to be adapted
  for(json& equallyRankedLigandsList : ligandsRanking) {
    for(json& oldLigandIndexJSON : equallyRankedLigandsList) {
      oldLigandIndexJSON = newLigandIndex(oldLigandIndexJSON.get<unsigned>());
    }

    // Sort the sub list, too
    temple::inplace::sort(equallyRankedLigandsList);
  }

  // Changing ligands also means that links' ligand indices have to be adapted
  if(ranking.count("lnk") > 0) {
    json& links = ranking.at("lnk");
    assert(links.is_array() && !links.empty());
    for(json& link : links) {
      assert(link.count("p") > 0);
      json& ligandIndexPair = link.at("p");
      assert(ligandIndexPair.is_array() && ligandIndexPair.size() == 2);

      ligandIndexPair[0] = newLigandIndex(ligandIndexPair[0].get<unsigned>());
      ligandIndexPair[1] = newLigandIndex(ligandIndexPair[1].get<unsigned>());

      // Reorder the pair
      if(ligandIndexPair[0] > ligandIndexPair[1]) {
        std::swap(ligandIndexPair[0], ligandIndexPair[1]);
      }
    }

    // Sort links by their index pair
    temple::inplace::sort(
      links,
      [](const json& lhs, const json& rhs) -> bool {
        assert(lhs.is_object() && lhs.count("p") > 0);
        assert(rhs.is_object() && rhs.count("p") > 0);

        return lhs["p"] < rhs["p"];
      }
    );
  }
}

void sortBondStereopermutatorsByEdge(nlohmann::json& m) {
  temple::inplace::sort(
    m.at(bondStereopermutatorKey),
    [&](const nlohmann::json& lhs, const nlohmann::json& rhs) -> bool {
      assert(lhs.count("e") > 0 && rhs.count("e") > 0);
      assert(lhs["e"].is_array() && rhs["e"].is_array());
      assert(lhs["e"].size() == 2 && rhs["e"].size() == 2);
      return (
        std::tie(lhs["e"][0], lhs["e"][1])
        < std::tie(rhs["e"][0], rhs["e"][1])
      );
    }
  );
}

void standardizeJSON(nlohmann::json& m) {
  using json = nlohmann::json;
  /* Ensure that m["a"], which is the list of atom stereopermutator objects,
   * is sorted by the central index
   */
  if(m.count(atomStereopermutatorKey) > 0) {
    sortAtomStereopermutatorsByCentralIndex(m);

    // Ranking objects of AtomStereopermutators have lots of notational freedom
    for(json& permutator : m.at(atomStereopermutatorKey)) {
      standardizeAtomStereopermutatorRepresentation(permutator);
    }
  }

  /* Ensure that m["b"], which is the list of bond stereopermutator objects,
   * is sorted by their placement edges
   */
  if(m.count(bondStereopermutatorKey) > 0) {
    sortBondStereopermutatorsByEdge(m);
  }

  /* Ensure that each graph edge is ordered, i.e. the first vertex index is
   * smaller than the second
   */
  for(json& edge : m.at(graphKey).at("E")) {
    assert(edge.is_array() && edge.size() == 3);
    if(edge[0] > edge[1]) {
      std::swap(edge[0], edge[1]);
    }
  }

  /* Sort the list of edges */
  temple::inplace::sort(m.at(graphKey).at("E"));
}

nlohmann::json serialize(const Molecule& molecule) {
  using json = nlohmann::json;

  json m;

  // Add a version tag. Always serialize to the newest version information
  m[versionKey] = {version::major, version::minor, version::patch};

  m[graphKey] = molecule.graph();

  const auto& stereopermutators = molecule.stereopermutators();

  /* Manual conversion of stereopermutators because these need the graph to be
   * constructed
   */
  if(stereopermutators.A() > 0) {
    // Atom stereopermutators first
    m["a"] = json::array();
    for(const auto& stereopermutator : stereopermutators.atomStereopermutators()) {
      json s;

      s["c"] = stereopermutator.centralIndex();
      s["s"] = shapes::nameIndex(stereopermutator.getShape());
      s["r"] = stereopermutator.getRanking();

      if(stereopermutator.assigned()) {
        s["a"] = stereopermutator.assigned().value();
      }

      m["a"].push_back(std::move(s));
    }
  }

  if(stereopermutators.B() > 0) {
    m["b"] = json::array();
    for(const auto& stereopermutator : stereopermutators.bondStereopermutators()) {
      json s;

      s["e"] = {stereopermutator.edge().first, stereopermutator.edge().second};

      if(stereopermutator.assigned()) {
        s["a"] = stereopermutator.assigned().value();
      }

      if(stereopermutator.alignment() == BondStereopermutator::Alignment::Staggered) {
        s["al"] = stereopermutator.alignment();
      }

      m["b"].push_back(std::move(s));
    }
  }

  if(auto canonicalComponentsOption = molecule.canonicalComponents()) {
    m[canonicalKey] = canonicalComponentsOption.value();
  }

  return m;
}

Molecule deserialize(const nlohmann::json& m) {
  /* In the future, we will have to check the version information to ensure the
   * chosen deserialization algorithm is compatible (in case breaking changes
   * are made). As long as we do not increment the version past major = 1,
   * though, this algorithm does not require any changes.
   */
  // Ensure the future hasn't come yet
  static_assert(
    version::major <= 1,
    "Incrementing the major version past 1 requires that you look at JSON deserialization!"
  );

  /* Look at the version information to determine if there is a suitable
   * deserialization algorithm or not
   */
  // std::vector<unsigned> version = m[versionKey];

  Graph graph = m.at(graphKey);

  StereopermutatorList stereopermutators;

  // Atom stereopermutators
  if(m.count(atomStereopermutatorKey) > 0) {
    for(const auto& j : m[atomStereopermutatorKey]) {
      shapes::Shape shape = shapes::allShapes.at(j["s"]);
      AtomIndex centralIndex = j["c"];

      auto stereopermutator = AtomStereopermutator {
        graph,
        shape,
        centralIndex,
        j["r"].get<RankingInformation>()
      };

      // Assign if present
      if(j.count("a") > 0) {
        stereopermutator.assign(
          static_cast<unsigned>(j["a"])
        );
      }

      stereopermutators.add(std::move(stereopermutator));
    }
  }

  if(m.count(bondStereopermutatorKey) > 0) {
    // Bond stereopermutators
    for(const auto& j : m[bondStereopermutatorKey]) {
      AtomIndex a = j["e"].at(0);
      AtomIndex b = j["e"].at(1);

      assert(stereopermutators.option(a) && stereopermutators.option(b));

      BondIndex molEdge {a, b};

      auto alignment = BondStereopermutator::Alignment::Eclipsed;
      if(j.count("al") > 0) {
        alignment = j["al"].get<BondStereopermutator::Alignment>();
      }

      auto stereopermutator = BondStereopermutator {
        graph.inner(),
        stereopermutators,
        molEdge,
        alignment
      };

      // Assign if present
      if(j.count("a") > 0) {
        stereopermutator.assign(
          static_cast<unsigned>(j["a"])
        );
      }

      stereopermutators.add(std::move(stereopermutator));
    }
  }

  boost::optional<AtomEnvironmentComponents> canonicalComponentsOption;
  if(m.count(canonicalKey) > 0) {
    canonicalComponentsOption = m[canonicalKey];
  }

  return Molecule {graph, stereopermutators, canonicalComponentsOption};
}

} // namespace detail

struct JsonSerialization::Impl {
  explicit Impl(const std::string& jsonString) : serialization(nlohmann::json::parse(jsonString)) {}
  explicit Impl(const Molecule& molecule) : serialization(detail::serialize(molecule)) {
    // If the molecule is fully canonical, make serializations canonical too
    if(molecule.canonicalComponents() == AtomEnvironmentComponents::All) {
      detail::standardizeJSON(serialization);
    }
  }

  Impl(
    const BinaryType& binary,
    const BinaryFormat format
  ) {
    if(format == BinaryFormat::CBOR) {
      serialization = nlohmann::json::from_cbor(binary);
    } else if(format == BinaryFormat::BSON) {
      serialization = nlohmann::json::from_bson(binary);
    } else if(format == BinaryFormat::MsgPack) {
      serialization = nlohmann::json::from_msgpack(binary);
    } else if(format == BinaryFormat::UBJSON) {
      serialization = nlohmann::json::from_ubjson(binary);
    }
  }

  operator std::string() const {
    return serialization.dump();
  }

  operator Molecule() const {
    return detail::deserialize(serialization);
  }

  BinaryType toBinary(const BinaryFormat format) {
    if(format == BinaryFormat::CBOR) {
      return nlohmann::json::to_cbor(serialization);
    }

    if(format == BinaryFormat::BSON) {
      return nlohmann::json::to_bson(serialization);
    }

    if(format == BinaryFormat::MsgPack) {
      return nlohmann::json::to_msgpack(serialization);
    }

    if(format == BinaryFormat::UBJSON) {
      return nlohmann::json::to_ubjson(serialization);
    }

    return {};
  }

  void standardize() {
    if(serialization.count(canonicalKey) == 0) {
      throw std::logic_error("Molecule is not canonical. Standardizing the JSON representation does not make sense.");
    }

    if(serialization.at(canonicalKey) != AtomEnvironmentComponents::All) {
      throw std::logic_error("Molecule is not fully canonical. Standardizing the JSON representation does not make sense.");
    }

    detail::standardizeJSON(serialization);
  }

  nlohmann::json serialization;
};

/* JsonSerialization implementation */

std::string JsonSerialization::base64Encode(const BinaryType& binary) {
  return base64::encode(binary);
}

JsonSerialization::BinaryType JsonSerialization::base64Decode(const std::string& binary) {
  return base64::decode(binary);
}

JsonSerialization::JsonSerialization(JsonSerialization&& other) noexcept = default;
JsonSerialization& JsonSerialization::operator = (JsonSerialization&& other) noexcept = default;

JsonSerialization::JsonSerialization(const JsonSerialization& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}

JsonSerialization& JsonSerialization::operator = (const JsonSerialization& other) {
  *_pImpl = *other._pImpl;
  return *this;
}

JsonSerialization::~JsonSerialization() = default;

JsonSerialization::JsonSerialization(const std::string& jsonString)
  : _pImpl(std::make_unique<Impl>(jsonString)) {}

JsonSerialization::JsonSerialization(const Molecule& molecule)
  : _pImpl(std::make_unique<Impl>(molecule)) {}

JsonSerialization::JsonSerialization(const BinaryType& binary, const BinaryFormat format)
  : _pImpl(std::make_unique<Impl>(binary, format)) {}

JsonSerialization::operator std::string() const {
  return _pImpl->operator std::string();
}

JsonSerialization::operator Molecule() const {
  return _pImpl->operator Molecule();
}

JsonSerialization::BinaryType JsonSerialization::toBinary(const BinaryFormat format) const {
  return _pImpl->toBinary(format);
}

JsonSerialization& JsonSerialization::standardize() {
  _pImpl->standardize();
  return *this;
}

} // namespace molassembler

} // namespace Scine
