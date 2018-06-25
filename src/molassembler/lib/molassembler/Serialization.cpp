#include "json/json.hpp"

#include "detail/Base64.h"
#include "Serialization.h"

namespace nlohmann {

template<>
struct adl_serializer<Delib::ElementType> {
  using Type = Delib::ElementType;

  static void to_json(json& j, const Type& value) {
    j = Delib::ElementInfo::symbol(value);
  }

  static void from_json(const json& j, Type& value) {
    value = Delib::ElementInfo::elementTypeForSymbol(
      j.get<std::string>()
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
    j["pair"] = json::array();
    j["pair"].push_back(link.indexPair.first);
    j["pair"].push_back(link.indexPair.second);

    j["sequence"] = link.cycleSequence;
  }

  static void from_json(const json& j, Type& link) {
    link.indexPair = {
      j["pair"].at(0),
      j["pair"].at(1)
    };

    link.cycleSequence.reserve(j["sequence"].size());
    for(const auto& sequenceElementJSON : j["sequence"]) {
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
    j["ligands"] = ranking.ligands;
    j["ligandsRanking"] = ranking.ligandsRanking;
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

    ranking.ligands.reserve(j["ligands"].size());
    for(const auto& listJSON : j["ligands"]) {
      std::vector<molassembler::AtomIndexType> ligandConstitutingAtoms;
      for(const auto& listElementJSON : listJSON) {
        ligandConstitutingAtoms.push_back(listElementJSON);
      }
      ranking.ligands.push_back(ligandConstitutingAtoms);
    }

    ranking.ligandsRanking.reserve(j["ligandsRanking"].size());
    for(const auto& listJSON : j["ligandsRanking"]) {
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
    j["elements"] = json::array();
    auto& elements = j["elements"];

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
      const auto edgeDescriptor :
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
    unsigned N = j["elements"].size();

    graph = Type (N);

    for(unsigned i = 0; i < N; ++i) {
      graph[i].elementType = j["elements"].at(i);
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

  m["graph"] = molecule.getGraph();

  // Manual conversion of stereocenters (need access to molecule members for deserialization)
  m["stereocenters"] = json::array();
  for(const auto& stereocenterPtr : molecule.getStereocenterList()) {
    json c;

    if(stereocenterPtr->type() == Stereocenters::Type::CNStereocenter) {
      c["type"] = "CN";
      auto cnPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        stereocenterPtr
      );

      c["symmetry"] = Symmetry::name(cnPtr->getSymmetry());
      c["ranking"] = cnPtr->getRanking();
    } else if(stereocenterPtr->type() == Stereocenters::Type::EZStereocenter) {
      c["type"] = "EZ";

      auto ezPtr = std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
        stereocenterPtr
      );
      c["leftRanking"] = ezPtr->getLeftRanking();
      c["rightRanking"] = ezPtr->getRightRanking();
    }

    if(stereocenterPtr->assigned()) {
      c["assignment"] = stereocenterPtr->assigned().value();
    }

    c["centers"] = stereocenterPtr->involvedAtoms();

    m["stereocenters"].push_back(std::move(c));
  }

  return m;
}

Molecule deserialize(const nlohmann::json& m) {
  GraphType graph = m["graph"];

  StereocenterList stereocenters;
  for(const auto& j : m["stereocenters"]) {
    if(j["type"] == "CN"s) {
      Symmetry::Name symmetry = Symmetry::nameFromString(j["symmetry"]);
      AtomIndexType centralIndex = j["centers"].front();

      auto cnPtr = std::make_shared<Stereocenters::CNStereocenter>(
        graph,
        symmetry,
        centralIndex,
        j["ranking"].get<RankingInformation>()
      );

      // Assign if present
      if(j.count("assignment") > 0) {
        cnPtr->assign(
          static_cast<unsigned>(j["assignment"])
        );
      }

      stereocenters.add(cnPtr);
    } else if(j["type"] == "EZ") {
      AtomIndexType left = j["centers"].front();
      AtomIndexType right = j["centers"].back();

      auto ezPtr = std::make_shared<Stereocenters::EZStereocenter>(
        left,
        j["leftRanking"].get<RankingInformation>(),
        right,
        j["rightRanking"].get<RankingInformation>()
      );

      // Assign if present
      if(j.count("assignment") > 0) {
        ezPtr->assign(
          static_cast<unsigned>(j["assignment"])
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
