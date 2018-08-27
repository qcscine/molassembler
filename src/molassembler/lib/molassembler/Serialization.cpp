#include "molassembler/Serialization.h"

#include "chemical_symmetries/Symmetries.h"
#include "boost/range/adaptor/map.hpp"
#include "Delib/ElementTypeCollection.h"
#include "json/json.hpp"

#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/IO/Base64.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/StereocenterList.h"

namespace nlohmann {

template<>
struct adl_serializer<Delib::ElementType> {
  using Type = Delib::ElementType;
  using Underlying = std::underlying_type<Delib::ElementType>::type;

  static void to_json(json& j, const Type& value) {
    j = static_cast<Underlying>(value);
  }

  static void from_json(const json& j, Type& value) {
    value = static_cast<Delib::ElementType>(
      j.get<Underlying>()
    );
  }
};

template<>
struct adl_serializer<Delib::ElementTypeCollection> {
  using Type = Delib::ElementTypeCollection;

  static void to_json(json& j, const Type& value) {
    j = json::array();
    for(const Delib::ElementType elementType : value) {
      j.push_back(elementType);
    }
  }

  static void from_json(const json& j, Type& value) {
    for(const auto& elementJSON : j) {
      value.push_back(
        elementJSON.get<Delib::ElementType>()
      );
    }
  }
};

template<>
struct adl_serializer<molassembler::LinkInformation> {
  using Type = molassembler::LinkInformation;

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
struct adl_serializer<molassembler::RankingInformation> {
  using Type = molassembler::RankingInformation;

  static void to_json(json& j, const Type& ranking) {
    j["sorted"] = ranking.sortedSubstituents;
    j["links"] = ranking.links;
    j["lig"] = ranking.ligands;
    j["ligRank"] = ranking.ligandsRanking;
  }

  static void from_json(const json& j, Type& ranking) {
    ranking.sortedSubstituents.reserve(j["sorted"].size());

    for(const auto& listJSON : j["sorted"]) {
      std::vector<molassembler::AtomIndex> subGroup;
      for(const auto& listElementJSON : listJSON) {
        subGroup.push_back(listElementJSON);
      }

      ranking.sortedSubstituents.push_back(std::move(subGroup));
    }

    ranking.links.reserve(j["links"].size());

    for(const auto& listJSON : j["links"]) {
      ranking.links.push_back(listJSON);
    }

    ranking.ligands.reserve(j["lig"].size());
    for(const auto& listJSON : j["lig"]) {
      std::vector<molassembler::AtomIndex> ligandConstitutingAtoms;
      for(const auto& listElementJSON : listJSON) {
        ligandConstitutingAtoms.push_back(listElementJSON);
      }
      ranking.ligands.push_back(ligandConstitutingAtoms);
    }

    ranking.ligandsRanking.reserve(j["ligRank"].size());
    for(const auto& listJSON : j["ligRank"]) {
      std::vector<unsigned> equalLigandIndices;
      for(const auto& listElementJSON : listJSON) {
        equalLigandIndices.push_back(listElementJSON);
      }
      ranking.ligandsRanking.push_back(equalLigandIndices);
    }
  }
};

template<>
struct adl_serializer<molassembler::OuterGraph> {
  using Type = molassembler::OuterGraph;

  static void to_json(json& j, const Type& graph) {
    const molassembler::InnerGraph& inner = graph.inner();

    j["Z"] = json::array();
    auto& elements = j["Z"];

    for(
      const auto vertexIndex :
      boost::make_iterator_range(
        inner.vertices()
      )
    ) {
      elements.push_back(
        inner.elementType(vertexIndex)
      );
    }


    j["edges"] = json::array();
    auto& edges = j["edges"];

    for(
      const molassembler::InnerGraph::Edge& edgeDescriptor :
      boost::make_iterator_range(
        inner.edges()
      )
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
    unsigned N = j["Z"].size();

    molassembler::InnerGraph inner (N);

    for(unsigned i = 0; i < N; ++i) {
      inner.elementType(i) = j["Z"].at(i);
    }

    for(const auto& edgeJSON : j["edges"]) {
      inner.addEdge(
        edgeJSON.at(0),
        edgeJSON.at(1),
        static_cast<molassembler::BondType>(
          edgeJSON.at(2)
        )
      );
    }

    graph = Type {std::move(inner)};
  }
};

} // namespace nlohmann

namespace molassembler {

namespace detail {

nlohmann::json serialize(const Molecule& molecule) {
  using json = nlohmann::json;

  json m;

  m["g"] = molecule.graph();

  const auto& stereocenters = molecule.stereocenters();

  // Manual conversion of stereocenters
  // Atom stereocenters first
  m["a"] = json::array();
  for(const auto& stereocenter : stereocenters.atomStereocenters()) {
    json s;

    s["c"] = stereocenter.centralIndex();
    s["sym"] = Symmetry::name(stereocenter.getSymmetry());
    s["rank"] = stereocenter.getRanking();

    if(stereocenter.assigned()) {
      s["assign"] = stereocenter.assigned().value();
    }

    m["a"].push_back(std::move(s));
  }

  m["b"] = json::array();
  for(const auto& stereocenter : stereocenters.bondStereocenters()) {
    json s;

    s["edge"] = {
      stereocenter.edge().first,
      stereocenter.edge().second
    };

    if(stereocenter.assigned()) {
      s["assign"] = stereocenter.assigned().value();
    }

    m["b"].push_back(std::move(s));
  }

  return m;
}

Molecule deserialize(const nlohmann::json& m) {
  OuterGraph graph = m["g"];

  StereocenterList stereocenters;

  // Atom stereocenters
  for(const auto& j : m["a"]) {
    Symmetry::Name symmetry = Symmetry::nameFromString(j["sym"]);
    AtomIndex centralIndex = j["centers"].front();

    auto stereocenter = AtomStereocenter {
      graph,
      symmetry,
      centralIndex,
      j["rank"].get<RankingInformation>()
    };

    // Assign if present
    if(j.count("assign") > 0) {
      stereocenter.assign(
        static_cast<unsigned>(j["assign"])
      );
    }

    stereocenters.add(centralIndex, std::move(stereocenter));
  }

  // Bond stereocenters
  for(const auto& j : m["b"]) {
    AtomIndex a = j["edge"].at(0);
    AtomIndex b = j["edge"].at(1);

    auto aStereocenterOption = stereocenters.option(a);
    auto bStereocenterOption = stereocenters.option(b);

    assert(aStereocenterOption && bStereocenterOption);

    BondIndex molEdge {a, b};

    auto stereocenter = BondStereocenter {
      aStereocenterOption.value(),
      bStereocenterOption.value(),
      molEdge
    };

    // Assign if present
    if(j.count("assign") > 0) {
      stereocenter.assign(
        static_cast<unsigned>(j["assign"])
      );
    }

    stereocenters.add(molEdge, std::move(stereocenter));
  }

  return Molecule {
    graph,
    stereocenters
  };
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
    nlohmann::json::from_cbor(
      cbor
    )
  );
}

Molecule fromBase64EncodedCBOR(const std::string& base64EncodedCBOR) {
  return fromCBOR(
    base64::decode(base64EncodedCBOR)
  );
}

} // namespace molassembler
