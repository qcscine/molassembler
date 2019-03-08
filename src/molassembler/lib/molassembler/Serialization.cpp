/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Serialization.h"

#include "chemical_symmetries/Symmetries.h"
#include "boost/range/adaptor/map.hpp"
#include "Utils/Typenames.h"
#include "json/json.hpp"

#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/IO/Base64.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Version.h"

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
  }
};

template<>
struct adl_serializer<Scine::molassembler::RankingInformation> {
  using Type = Scine::molassembler::RankingInformation;

  static void to_json(json& j, const Type& ranking) {
    j["s"] = ranking.sortedSubstituents;
    // Omit links member if the list is empty (common)
    if(!ranking.links.empty()) {
      j["lnk"] = ranking.links;
    }
    j["l"] = ranking.ligands;
    j["lr"] = ranking.ligandsRanking;
  }

  static void from_json(const json& j, Type& ranking) {
    ranking.sortedSubstituents.reserve(j["s"].size());

    for(const auto& listJSON : j["s"]) {
      std::vector<Scine::molassembler::AtomIndex> subGroup;
      for(const auto& listElementJSON : listJSON) {
        subGroup.push_back(listElementJSON);
      }

      ranking.sortedSubstituents.push_back(std::move(subGroup));
    }

    if(j.count("lnk") > 0) {
      ranking.links.reserve(j["lnk"].size());

      for(const auto& listJSON : j["lnk"]) {
        ranking.links.push_back(listJSON);
      }
    }

    ranking.ligands.reserve(j["l"].size());
    for(const auto& listJSON : j["l"]) {
      std::vector<Scine::molassembler::AtomIndex> ligandConstitutingAtoms;
      for(const auto& listElementJSON : listJSON) {
        ligandConstitutingAtoms.push_back(listElementJSON);
      }
      ranking.ligands.push_back(ligandConstitutingAtoms);
    }

    ranking.ligandsRanking.reserve(j["lr"].size());
    for(const auto& listJSON : j["lr"]) {
      std::vector<unsigned> equalLigandIndices;
      for(const auto& listElementJSON : listJSON) {
        equalLigandIndices.push_back(listElementJSON);
      }
      ranking.ligandsRanking.push_back(equalLigandIndices);
    }
  }
};

template<>
struct adl_serializer<Scine::molassembler::OuterGraph> {
  using Type = Scine::molassembler::OuterGraph;

  static void to_json(json& j, const Type& graph) {
    const Scine::molassembler::InnerGraph& inner = graph.inner();

    j["Z"] = json::array();
    auto& elements = j["Z"];

    for(const auto vertexIndex : boost::make_iterator_range(inner.vertices())) {
      elements.push_back(
        inner.elementType(vertexIndex)
      );
    }


    j["E"] = json::array();
    auto& edges = j["E"];

    for(
      const Scine::molassembler::InnerGraph::Edge& edgeDescriptor :
      boost::make_iterator_range(inner.edges())
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

    Scine::molassembler::InnerGraph inner (N);

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

namespace detail {

nlohmann::json serialize(const Molecule& molecule) {
  using json = nlohmann::json;

  json m;

  // Add a version tag. Always serialize to the newest version information
  m["v"] = {version::major, version::minor, version::patch};

  m["g"] = molecule.graph();

  const auto& stereopermutators = molecule.stereopermutators();

  // Manual conversion of stereopermutators
  // Atom stereopermutators first
  m["a"] = json::array();
  for(const auto& stereopermutator : stereopermutators.atomStereopermutators()) {
    json s;

    s["c"] = stereopermutator.centralIndex();
    s["s"] = Symmetry::nameIndex(stereopermutator.getSymmetry());
    s["r"] = stereopermutator.getRanking();

    if(stereopermutator.assigned()) {
      s["a"] = stereopermutator.assigned().value();
    }

    m["a"].push_back(std::move(s));
  }

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

  m["c"] = molecule.canonicalComponents();

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
  // std::vector<unsigned> version = m["v"];

  OuterGraph graph = m["g"];

  StereopermutatorList stereopermutators;

  // Atom stereopermutators
  for(const auto& j : m["a"]) {
    Symmetry::Name symmetry = Symmetry::allNames.at(j["s"]);
    AtomIndex centralIndex = j["c"];

    auto stereopermutator = AtomStereopermutator {
      graph,
      symmetry,
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

  // Bond stereopermutators
  for(const auto& j : m["b"]) {
    AtomIndex a = j["e"].at(0);
    AtomIndex b = j["e"].at(1);

    auto aStereopermutatorOption = stereopermutators.option(a);
    auto bStereopermutatorOption = stereopermutators.option(b);

    assert(aStereopermutatorOption && bStereopermutatorOption);

    BondIndex molEdge {a, b};

    auto alignment = BondStereopermutator::Alignment::Eclipsed;
    if(j.count("al") > 0) {
      alignment = j["al"].get<BondStereopermutator::Alignment>();
    }

    auto stereopermutator = BondStereopermutator {
      aStereopermutatorOption.value(),
      bStereopermutatorOption.value(),
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

  AtomEnvironmentComponents canonicalComponents = m["c"];

  return Molecule {graph, stereopermutators, canonicalComponents};
}

} // namespace detail

std::string toJSON(const Molecule& molecule) {
  return detail::serialize(molecule).dump();
}

std::vector<std::uint8_t> toCBOR(const Molecule& molecule) {
  return nlohmann::json::to_cbor(
    detail::serialize(molecule)
  );
}

std::string toBase64EncodedCBOR(const Molecule& molecule) {
  return base64::encode(
    toCBOR(molecule)
  );
}

Molecule fromJSON(const std::string& serializedMolecule) {
  return detail::deserialize(
    nlohmann::json::parse(serializedMolecule)
  );
}

Molecule fromCBOR(const std::vector<std::uint8_t>& cbor) {
  return detail::deserialize(
    nlohmann::json::from_cbor(cbor)
  );
}

Molecule fromBase64EncodedCBOR(const std::string& base64EncodedCBOR) {
  return fromCBOR(
    base64::decode(base64EncodedCBOR)
  );
}

} // namespace molassembler

} // namespace Scine
