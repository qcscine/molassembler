/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Serialization.h"

#include "Molassembler/Shapes/Data.h"
#include "boost/range/adaptor/map.hpp"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"
#include "Utils/Typenames.h"
#include "nlohmann/json.hpp"

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/IO/Base64.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/RankingInformation.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Version.h"

#include "Molassembler/Temple/Functional.h"

namespace nlohmann {

template<>
struct adl_serializer<Scine::Molassembler::SiteIndex> {
  using Type = Scine::Molassembler::SiteIndex;

  static void to_json(json& j, const Type& value) {
    j = static_cast<unsigned>(value);
  }

  static void from_json(const json& j, Type& value) {
    value = Type(j.get<unsigned>());
  }
};

template<>
struct adl_serializer<Scine::Molassembler::BondStereopermutator::Alignment> {
  using Type = Scine::Molassembler::BondStereopermutator::Alignment;
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
struct adl_serializer<Scine::Molassembler::AtomEnvironmentComponents> {
  using Type = Scine::Molassembler::AtomEnvironmentComponents;
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
struct adl_serializer<Scine::Molassembler::RankingInformation::Link> {
  using Type = Scine::Molassembler::RankingInformation::Link;

  static void to_json(json& j, const Type& link) {
    j["p"] = json::array();
    j["p"].push_back(link.sites.first);
    j["p"].push_back(link.sites.second);

    j["seq"] = link.cycleSequence;
  }

  static void from_json(const json& j, Type& link) {
    link.sites = {
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
struct adl_serializer<Scine::Molassembler::RankingInformation> {
  using Type = Scine::Molassembler::RankingInformation;

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
      std::vector<Scine::Molassembler::AtomIndex> subGroup;
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
      std::vector<Scine::Molassembler::AtomIndex> siteConstitutingAtoms;
      for(const auto& listElementJSON : listJSON) {
        siteConstitutingAtoms.push_back(listElementJSON);
      }
      ranking.sites.push_back(siteConstitutingAtoms);
    }

    ranking.siteRanking.reserve(j["lr"].size());
    for(const auto& listJSON : j["lr"]) {
      std::vector<Scine::Molassembler::SiteIndex> equalSiteIndices;
      for(const auto& listElementJSON : listJSON) {
        equalSiteIndices.push_back(listElementJSON);
      }
      ranking.siteRanking.push_back(equalSiteIndices);
    }
  }
};

template<>
struct adl_serializer<Scine::Molassembler::Graph> {
  using Type = Scine::Molassembler::Graph;

  static void to_json(json& j, const Type& graph) {
    const Scine::Molassembler::PrivateGraph& inner = graph.inner();

    j["Z"] = json::array();
    auto& elements = j["Z"];

    for(const auto vertexIndex : inner.vertices()) {
      elements.push_back(inner.elementType(vertexIndex));
    }


    j["E"] = json::array();
    auto& edges = j["E"];

    for(
      const Scine::Molassembler::PrivateGraph::Edge& edgeDescriptor :
      inner.edges()
    ) {
      json e = json::array();
      e.push_back(
        static_cast<int>(inner.source(edgeDescriptor))
      );
      e.push_back(
        static_cast<int>(inner.target(edgeDescriptor))
      );
      e.push_back(
        static_cast<int>(inner.bondType(edgeDescriptor))
      );

      edges.push_back(std::move(e));
    }
  }

  static void from_json(const json& j, Type& graph) {
    const unsigned N = j["Z"].size();

    Scine::Molassembler::PrivateGraph inner (N);

    for(unsigned i = 0; i < N; ++i) {
      inner.elementType(i) = j["Z"].at(i);
    }

    for(const auto& edgeJSON : j["E"]) {
      inner.addEdge(
        edgeJSON.at(0),
        edgeJSON.at(1),
        static_cast<Scine::Molassembler::BondType>(edgeJSON.at(2))
      );
    }

    graph = Type {std::move(inner)};
  }
};

} // namespace nlohmann

namespace Scine {

namespace Molassembler {

// Top level keys into molecule json object, for easier reading
constexpr const char* atomStereopermutatorKey = "a";
constexpr const char* bondStereopermutatorKey = "b";
constexpr const char* graphKey = "g";
constexpr const char* versionKey = "v";
constexpr const char* canonicalKey = "c";

namespace {

void sortAtomStereopermutatorsByCentralIndex(nlohmann::json& m) {
  Temple::sort(
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
    Temple::sort(subList);
  }

  // Sort sites' sub-lists
  assert(ranking.count("l") > 0);
  json& sites = ranking.at("l");
  assert(sites.is_array());
  for(json& subList : sites) {
    assert(subList.is_array());
    Temple::sort(subList);
  }

  /* Sort site index lists lexicographically (and their ranking
   * simultaneously), since they are index-connected
   */
  assert(ranking.count("lr") > 0);
  json& sitesRanking = ranking.at("lr");
  assert(sitesRanking.is_array());

  json unsortedSites = sites;
  Temple::sort(sites);

  auto newSiteIndex = [&](unsigned oldSiteIndex) -> unsigned {
    const auto& siteIndices = unsortedSites.at(oldSiteIndex);
    auto findIter = Temple::find(sites, siteIndices);
    assert(findIter != std::end(sites));
    return findIter - std::begin(sites);
  };

  // Changing sites means sitesRanking has to be adapted
  for(json& equallyRankedLigandsList : sitesRanking) {
    for(json& oldLigandIndexJSON : equallyRankedLigandsList) {
      oldLigandIndexJSON = newSiteIndex(oldLigandIndexJSON.get<unsigned>());
    }

    // Sort the sub list, too
    Temple::sort(equallyRankedLigandsList);
  }

  // Changing sites also means that links' site indices have to be adapted
  if(ranking.count("lnk") > 0) {
    json& links = ranking.at("lnk");
    assert(links.is_array() && !links.empty());
    for(json& link : links) {
      assert(link.count("p") > 0);
      json& siteIndexPair = link.at("p");
      assert(siteIndexPair.is_array() && siteIndexPair.size() == 2);

      siteIndexPair[0] = newSiteIndex(siteIndexPair[0].get<unsigned>());
      siteIndexPair[1] = newSiteIndex(siteIndexPair[1].get<unsigned>());

      // Reorder the pair
      if(siteIndexPair[0] > siteIndexPair[1]) {
        std::swap(siteIndexPair[0], siteIndexPair[1]);
      }
    }

    // Sort links by their index pair
    Temple::sort(
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
  Temple::sort(
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
  Temple::sort(m.at(graphKey).at("E"));
}

nlohmann::json serialize(const Molecule& molecule) {
  using json = nlohmann::json;

  json m;

  // Add a version tag. Always serialize to the newest version information
  m[versionKey] = {Version::major, Version::minor, Version::patch};

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

      s["c"] = stereopermutator.placement();
      s["s"] = Shapes::nameIndex(stereopermutator.getShape());
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

      s["e"] = {stereopermutator.placement().first, stereopermutator.placement().second};

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
    Version::major <= 1,
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
      Shapes::Shape shape = Shapes::allShapes.at(j["s"]);
      AtomIndex placement = j["c"];

      AtomStereopermutator stereopermutator {
        placement,
        shape,
        j["r"].get<RankingInformation>(),
        Stereopermutators::Feasible::Functor(graph),
        AtomStereopermutator::thermalizationFunctor(graph)
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

      BondStereopermutator stereopermutator {
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

} // namespace

struct JsonSerialization::Impl {
  explicit Impl(const std::string& jsonString) : serialization(nlohmann::json::parse(jsonString)) {}
  explicit Impl(const Molecule& molecule) : serialization(serialize(molecule)) {
    // If the molecule is fully canonical, make serializations canonical too
    if(molecule.canonicalComponents() == AtomEnvironmentComponents::All) {
      standardizeJSON(serialization);
    }
  }

  Impl(const BinaryType& binary, const BinaryFormat format)
    : serialization(fromBinary(binary, format)) {}

  operator std::string() const {
    return serialization.dump();
  }

  operator Molecule() const {
    return deserialize(serialization);
  }

  static nlohmann::json fromBinary(const BinaryType& binary, const BinaryFormat format) {
    switch(format) {
      case BinaryFormat::CBOR: return nlohmann::json::from_cbor(binary);
      case BinaryFormat::BSON: return nlohmann::json::from_bson(binary);
      case BinaryFormat::MsgPack: return nlohmann::json::from_msgpack(binary);
      case BinaryFormat::UBJSON: return nlohmann::json::from_ubjson(binary);
      default: throw std::runtime_error("Unknown binary format!");
    }
  }

  BinaryType toBinary(const BinaryFormat format) const {
    switch(format) {
      case BinaryFormat::CBOR: return nlohmann::json::to_cbor(serialization);
      case BinaryFormat::BSON: return nlohmann::json::to_bson(serialization);
      case BinaryFormat::MsgPack: return nlohmann::json::to_msgpack(serialization);
      case BinaryFormat::UBJSON: return nlohmann::json::to_ubjson(serialization);
      default: return {};
    }
  }

  void standardize() {
    if(serialization.count(canonicalKey) == 0) {
      throw std::logic_error("Molecule is not canonical. Standardizing the JSON representation does not make sense.");
    }

    if(serialization.at(canonicalKey) != AtomEnvironmentComponents::All) {
      throw std::logic_error("Molecule is not fully canonical. Standardizing the JSON representation does not make sense.");
    }

    standardizeJSON(serialization);
  }

  nlohmann::json serialization;
};

/* JsonSerialization implementation */

std::string JsonSerialization::base64Encode(const BinaryType& binary) {
  return base64::encode(binary);
}

JsonSerialization::BinaryType JsonSerialization::base64Decode(const std::string& base64String) {
  return base64::decode(base64String);
}

JsonSerialization::JsonSerialization(JsonSerialization&& other) noexcept = default;
JsonSerialization& JsonSerialization::operator = (JsonSerialization&& other) noexcept = default;

JsonSerialization::JsonSerialization(const JsonSerialization& other) : pImpl_(
  std::make_unique<Impl>(*other.pImpl_)
) {}

JsonSerialization& JsonSerialization::operator = (const JsonSerialization& other) {
  *pImpl_ = *other.pImpl_;
  return *this;
}

JsonSerialization::~JsonSerialization() = default;

JsonSerialization::JsonSerialization(const std::string& jsonString)
  : pImpl_(std::make_unique<Impl>(jsonString)) {}

JsonSerialization::JsonSerialization(const Molecule& molecule)
  : pImpl_(std::make_unique<Impl>(molecule)) {}

JsonSerialization::JsonSerialization(const BinaryType& binary, const BinaryFormat format)
  : pImpl_(std::make_unique<Impl>(binary, format)) {}

JsonSerialization::operator std::string() const {
  return pImpl_->operator std::string();
}

JsonSerialization::operator Molecule() const {
  return pImpl_->operator Molecule();
}

JsonSerialization::BinaryType JsonSerialization::toBinary(const BinaryFormat format) const {
  return pImpl_->toBinary(format);
}

JsonSerialization& JsonSerialization::standardize() {
  pImpl_->standardize();
  return *this;
}

} // namespace Molassembler

} // namespace Scine
