#include "json/json.hpp"

#include "detail/Base64.h"
#include "Serialization.h"

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
struct adl_serializer<molassembler::GraphAlgorithms::LinkInformation> {
  using Type = molassembler::GraphAlgorithms::LinkInformation;

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
      std::vector<molassembler::AtomIndexType> subGroup;
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
      std::vector<molassembler::AtomIndexType> ligandConstitutingAtoms;
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
struct adl_serializer<molassembler::GraphType> {
  using Type = molassembler::GraphType;

  static void to_json(json& j, const Type& graph) {
    j["Z"] = json::array();
    auto& elements = j["Z"];

    for(
      const auto vertexIndex :
      RangeForTemporary<Type::vertex_iterator>(
        boost::vertices(graph)
      )
    ) {
      elements.push_back(
        graph[vertexIndex].elementType
      );
    }


    j["edges"] = json::array();
    auto& edges = j["edges"];

    for(
      const auto& edgeDescriptor :
      RangeForTemporary<Type::edge_iterator>(
        boost::edges(graph)
      )
    ) {
      json e = json::array();
      e.push_back(
        static_cast<int>(
          boost::source(edgeDescriptor, graph)
        )
      );
      e.push_back(
        static_cast<int>(
          boost::target(edgeDescriptor, graph)
        )
      );
      e.push_back(
        static_cast<int>(
          graph[edgeDescriptor].bondType
        )
      );

      edges.push_back(std::move(e));
    }
  }

  static void from_json(const json& j, Type& graph) {
    unsigned N = j["Z"].size();

    graph = Type (N);

    for(unsigned i = 0; i < N; ++i) {
      graph[i].elementType = j["Z"].at(i);
    }

    for(const auto& edgeJSON : j["edges"]) {
      auto addEdgePair = boost::add_edge(
        edgeJSON.at(0),
        edgeJSON.at(1),
        graph
      );

      assert(addEdgePair.second);

      graph[addEdgePair.first].bondType = static_cast<molassembler::BondType>(
        edgeJSON.at(2)
      );
    }
  }
};

} // namespace nlohmann

namespace molassembler {

namespace detail {

nlohmann::json serialize(const Molecule& molecule) {
  using json = nlohmann::json;

  json m;

  m["g"] = molecule.getGraph();

  // Manual conversion of stereocenters (need access to molecule members for deserialization)
  m["s"] = json::array();
  for(const auto& stereocenterPtr : molecule.getStereocenterList()) {
    json c;

    if(stereocenterPtr->type() == Stereocenters::Type::CNStereocenter) {
      auto cnPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        stereocenterPtr
      );

      c["sym"] = Symmetry::name(cnPtr->getSymmetry());
      c["rank"] = cnPtr->getRanking();
    } else if(stereocenterPtr->type() == Stereocenters::Type::EZStereocenter) {
      auto ezPtr = std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
        stereocenterPtr
      );
      c["lRank"] = ezPtr->getLeftRanking();
      c["rRank"] = ezPtr->getRightRanking();
    }

    if(stereocenterPtr->assigned()) {
      c["assign"] = stereocenterPtr->assigned().value();
    }

    c["centers"] = stereocenterPtr->involvedAtoms();

    m["s"].push_back(std::move(c));
  }

  return m;
}

Molecule deserialize(const nlohmann::json& m) {
  GraphType graph = m["g"];

  StereocenterList stereocenters;
  for(const auto& j : m["s"]) {
    if(j.count("sym") > 0) { // CNStereocenter
      Symmetry::Name symmetry = Symmetry::nameFromString(j["sym"]);
      AtomIndexType centralIndex = j["centers"].front();

      auto cnPtr = std::make_shared<Stereocenters::CNStereocenter>(
        graph,
        symmetry,
        centralIndex,
        j["rank"].get<RankingInformation>()
      );

      // Assign if present
      if(j.count("assign") > 0) {
        cnPtr->assign(
          static_cast<unsigned>(j["assign"])
        );
      }

      stereocenters.add(cnPtr);
    } else { // EZStereocenter
      AtomIndexType left = j["centers"].front();
      AtomIndexType right = j["centers"].back();

      auto ezPtr = std::make_shared<Stereocenters::EZStereocenter>(
        left,
        j["lRank"].get<RankingInformation>(),
        right,
        j["rRank"].get<RankingInformation>()
      );

      // Assign if present
      if(j.count("assign") > 0) {
        ezPtr->assign(
          static_cast<unsigned>(j["assign"])
        );
      }

      stereocenters.add(ezPtr);
    }
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
