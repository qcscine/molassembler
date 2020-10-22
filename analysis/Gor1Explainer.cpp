/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/graph/two_bit_color_map.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"

// NOTE: Include order is important here!
#include "Molassembler/DistanceGeometry/ImplicitBoundsGraphBoost.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"
#include "Molassembler/IO.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Graph/Gor1.h"
#include "boost/graph/graphviz.hpp"
#include "Molassembler/Temple/StringAlgorithms.h"

#include <random>

using namespace Scine;
using namespace Molassembler;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

namespace Scine {

namespace Molassembler {

namespace DistanceGeometry {

struct WriterState {
  using Graph = ImplicitBoundsGraph;
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

  boost::optional<Edge> highlightEdge = boost::none;
  boost::optional<Vertex> highlightVertex = boost::none;
};

class GraphvizWriter {
public:
  using Graph = ImplicitBoundsGraph;
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<Graph>::edge_descriptor;
  using ColorMapType = boost::two_bit_color_map<>;
  using Color = boost::color_traits<typename boost::property_traits<ColorMapType>::value_type>;

private:
  const Graph& graph_;
  const std::vector<double>& distances_;
  const std::vector<Vertex>& predecessors_;
  const ColorMapType& colors_;

public:
  GraphvizWriter(
    const Graph& graph,
    const std::vector<double>& distances,
    const std::vector<Vertex>& predecessors,
    const ColorMapType& colors
  ) : graph_ {graph},
      distances_ {distances},
      predecessors_ {predecessors},
      colors_ {colors}
  {}

  std::shared_ptr<WriterState> statePtr = std::make_shared<WriterState>();

  void write_clusters(std::ostream& os) {
    os << "\nsubgraph cluster_left {\n";
    for(Vertex i = 0; i < boost::num_vertices(graph_); i += 2) {
      os << "  " << i << "; ";
    }
    os << "\n}\nsubgraph cluster_right {\n";
    for(Vertex i = 1; i < boost::num_vertices(graph_); i += 2) {
      os << "  " << i << "; ";
    }
    os << "\n}\n";
  }

  void operator() (std::ostream& os) {
    os << R"(  graph [fontname = "Arial"];)" << "\n"
      << R"(  node [fontname = "Arial", shape = "circle", style = "filled"];)" << "\n"
      << R"(  edge [fontname = "Arial"];)" << "\n"
      << R"(  rankdir="LR";)" << nl;

    os << std::fixed << std::setprecision(2);
  }

  bool inShortestPathsTree(const Vertex& v) {
    return predecessors_.at(v) != v;
  }

  void operator() (std::ostream& os, const Vertex& v) {
    os << "[";


    auto v_color = boost::get(colors_, v);
    if(statePtr->highlightVertex && statePtr->highlightVertex.value() == v) {
      os << R"(fillcolor="tomato", fontcolor="white")";
      statePtr->highlightVertex = boost::none;
    } else if(v_color == Color::white()) {
      os << R"(fillcolor="white")";
    } else if(v_color == Color::gray()) {
      os << R"(fillcolor="gray60", fontcolor="white")";
    } else if(v_color == Color::black()) {
      os << R"(fillcolor="black", fontcolor="white")";
    }

    if(inShortestPathsTree(v)) {
      os << R"(, color="steelblue", penwidth=2)";
    }

    // Index and distance label
    os << R"(, label=")" << v << nl;

    if(distances_.at(v) < 10000) {
      os << distances_.at(v);
    } else {
      os << "inf";
    }
    os << "\"";

    os << "]";

    if(boost::num_vertices(graph_) <= 20 && v == boost::num_vertices(graph_) - 1) {
      write_clusters(os);
    }
  }

  bool inShortestPathsTree(const Vertex& source, const Vertex& target) {
    return (
      predecessors_.at(target) == source
      || (
        target % 2 == source % 2
        && predecessors_.at(source) == target
      )
    );
  }

  void operator() (std::ostream& os, const Edge& e) {
    os << "[";

    os << R"(label=")" << boost::get(boost::edge_weight, graph_)(e) << R"(")";

    auto source = boost::source(e, graph_);
    auto target = boost::target(e, graph_);

    bool isHighlightedEdge = false;
    if(statePtr->highlightEdge) {
      auto hSource = boost::source(statePtr->highlightEdge.value(), graph_);
      auto hTarget = boost::target(statePtr->highlightEdge.value(), graph_);

      if(
        (source == hSource && target == hTarget)
        || (source == hTarget && target == hSource)
      ) {
        os << R"(, color="tomato", penwidth=4)";
        statePtr->highlightEdge = boost::none;
        isHighlightedEdge = true;
      }
    }
    if(!isHighlightedEdge && inShortestPathsTree(source, target)) {
      os << R"(, color="steelblue", penwidth=3)";
    }

    if(source % 2 != target % 2 && !graph_.hasExplicit(e)) {
      os << R"(, style="dashed")";
    }

    if(source % 2 == target % 2) { // Within-group edges are bidirectional
      if(source < target) {
        // Show the edge
        os << R"(, dir="both")";
      } else {
        // Hide it
        os << R"(, constraint=false, style=invis)";
      }
    }

    os << "]";
  }
};

template<typename T>
void write_graphviz(std::ostream& os, const std::deque<T>& a, const std::string& name) {
  os << "digraph G {";
  os << R"(  graph [fontname = "Arial"];)" << nl
    << R"(  node [fontname = "Arial", shape = "circle", style = "filled"];)" << nl
    << R"(  edge [fontname = "Arial"];)" << nl;
  os << R"(  0 [shape="square", label=")" << name << R"("];)" << nl;

  unsigned previousNode = 0;
  for(const auto& v : a) {
    os << "  " << (previousNode + 1) << R"([label=")" << v << R"("];)" << nl;
    os << "  " << previousNode << " -> " << (previousNode + 1) << ";" << nl;

    ++previousNode;
  }

  os << "}";
}

template<typename Vertex>
void write_predecessor_graphviz(std::ostream& os, const std::vector<Vertex>& predecessors) {
  os << "digraph G {";
  os << R"(  graph [fontname = "Arial"];)" << nl
    << R"(  node [fontname = "Arial", shape="circle"];)" << nl;
  for(unsigned i = 0; i < predecessors.size(); ++i) {

    if(predecessors.at(i) != i) {
      os << "  " << i << R"([label=")" << i << R"("];)" << nl;
      os << "  " << predecessors.at(i) << " -> " << i << ";" << nl;
    } else {
      // Check if i is listed somewhere else in the predecessor list
      if(std::count(predecessors.begin(), predecessors.end(), i) > 1) {
        os << "  " << i << R"([label=")" << i << R"("];)" << nl;
      }
    }
  }
  os << "}";
}

class VisualizationVisitor {
public:
  using Graph = ImplicitBoundsGraph;
  using Vertex = typename boost::graph_traits<ImplicitBoundsGraph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<ImplicitBoundsGraph>::edge_descriptor;

private:
  const std::vector<Vertex>& predecessors_;

  GraphvizWriter& writer_;
  unsigned counter_ = 0;
  std::deque<Vertex> A_, B_;

public:
  VisualizationVisitor(
    const std::vector<Vertex>& predecessors,
    GraphvizWriter& writer
  ) : predecessors_ {predecessors},
      writer_ {writer}
  {}

  void examine_edge(const Edge& e, const Graph& g) {
    writer_.statePtr->highlightEdge = e;
    makeGraphvizFiles(g);
  }

  void relax_edge(const Edge& e, const Graph& g) {
    writer_.statePtr->highlightEdge = e;
    makeGraphvizFiles(g);
  }

  void a_push(const Vertex& v, const Graph& g) {
    writer_.statePtr->highlightVertex = v;
    A_.push_front(v);
    makeGraphvizFiles(g);
  }
  void a_pop(const Vertex& v, const Graph& g) {
    writer_.statePtr->highlightVertex = v;
    A_.pop_front();
    makeGraphvizFiles(g);
  }
  void b_push(const Vertex& v, const Graph& g) {
    writer_.statePtr->highlightVertex = v;
    B_.push_front(v);
    makeGraphvizFiles(g);
  }
  void b_pop(const Vertex& v, const Graph& g) {
    writer_.statePtr->highlightVertex = v;
    B_.pop_front();
    makeGraphvizFiles(g);
  }

  void mark_white(const Vertex& /* v */, const Graph& g) {
    makeGraphvizFiles(g);
  }
  void mark_gray(const Vertex& /* v */, const Graph& g) {
    makeGraphvizFiles(g);
  }
  void mark_black(const Vertex& /* v */, const Graph& g) {
    makeGraphvizFiles(g);
  }

  void makeGraphvizFiles(const Graph& g) {
    using namespace std::string_literals;

    const std::string prefix = "gor1-explainer-"s
      + std::to_string(counter_)
      + "-"s;

    std::ofstream graphFile(prefix + "0.dot"s);

    boost::write_graphviz(
      graphFile,
      g,
      writer_,
      writer_,
      writer_
    );

    graphFile.close();

    std::ofstream stackAFile(prefix + "1.dot"s);
    write_graphviz(stackAFile, A_, "A");
    stackAFile.close();

    std::ofstream stackBFile(prefix + "2.dot"s);
    write_graphviz(stackBFile, B_, "B");
    stackBFile.close();

    std::ofstream treeGraph(prefix + "3.dot"s);
    write_predecessor_graphviz(treeGraph, predecessors_);
    treeGraph.close();

    ++counter_;
  }
};

} // namespace DistanceGeometry

} // namespace Molassembler

} // namespace Scine

int main(int argc, char* argv[]) {
  using namespace std::string_literals;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    (
      "f",
      boost::program_options::value<std::string>(),
      "Read molecule to generate from file"
    )
    (
      "r",
      boost::program_options::value<unsigned>(),
      "Root index to run Gor from"
    )
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  unsigned root_vertex = 0;
  if(options_variables_map.count("r") > 0) {
    root_vertex = options_variables_map["r"].as<unsigned>();
  }

  if(options_variables_map.count("help") > 0) {
    std::cout << options_description << nl;
    return 0;
  }

  if(options_variables_map.count("f") > 0) {
    auto filename = options_variables_map["f"].as<std::string>();

    if(!boost::filesystem::exists(filename)) {
      std::cout << "The specified file could not be found!" << nl;
      return 1;
    }

    auto mol = IO::read(filename);

    DistanceGeometry::SpatialModel spatialModel {mol, DistanceGeometry::Configuration {}};

    DistanceGeometry::ImplicitBoundsGraph shortestPathsGraph {
      mol.graph().inner(),
      spatialModel.makePairwiseBounds()
    };

    /* Prep */
    using Vertex = typename boost::graph_traits<DistanceGeometry::ImplicitBoundsGraph>::vertex_descriptor;

    unsigned N = boost::num_vertices(shortestPathsGraph);
    std::vector<double> distances(N);
    std::vector<Vertex> predecessors(N);

    // Initialize every element to itself for interpretation purposes! Not necessary by default for the algorithm
    std::iota(
      predecessors.begin(),
      predecessors.end(),
      0
    );

    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      get(boost::vertex_index, shortestPathsGraph)
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      get(boost::vertex_index, shortestPathsGraph)
    );

    using ColorMapType = boost::two_bit_color_map<>;
    ColorMapType color_map {N};

    DistanceGeometry::GraphvizWriter writer {
      shortestPathsGraph,
      distances,
      predecessors,
      color_map
    };

    DistanceGeometry::VisualizationVisitor visitor {
      predecessors,
      writer
    };

    /* Execution */
    boost::gor1_simplified_shortest_paths(
      shortestPathsGraph,
      Vertex {root_vertex},
      predecessor_map,
      color_map,
      distance_map,
      visitor
    );

    std::cout << mol << std::endl;

    auto filepath = boost::filesystem::path(filename);
    std::string folderName = "gor1-"s + filepath.stem().string();

    // Create a new ranking-results folder if it doesn't exist
    if(!boost::filesystem::is_directory(folderName)) {
      boost::filesystem::create_directory(folderName);
    }

    const boost::regex fileFilterRegex {R"(gor1-explainer-[0-9]{1,}-[0-9]\.dot)"};

    boost::filesystem::directory_iterator endIter;

    std::map<
      unsigned, // step index
      unsigned // number of graphs - 1
    > numGraphsMap;

    for(
      boost::filesystem::directory_iterator iter {"."};
      iter != endIter;
      ++iter
    ) {
      if(!boost::filesystem::is_regular_file(iter->status())) {
        continue;
      }

      // Move any generated ranking-tree files to the new folder
      boost::smatch what;

      if(!boost::regex_match(iter -> path().filename().string(), what, fileFilterRegex)) {
        continue;
      }

      boost::filesystem::path newPath {"./"s + folderName};
      newPath /= iter->path().filename();

      auto splat = Temple::split(iter->path().filename().string(), '-');
      auto step = std::stoul(splat.at(2));
      auto graphIndex = std::stoul(splat.at(3));

      // Collect information on how many graphs were generated for each step
      if(numGraphsMap.count(step) == 0) {
        numGraphsMap.emplace(
          step,
          graphIndex
        );
      } else if(numGraphsMap.at(step) < graphIndex) {
        numGraphsMap.at(step) = graphIndex;
      }

      // Move the file
      boost::filesystem::rename(
        iter->path(),
        newPath
      );

    }

    // Write a bash file for the generation of combined graphs
    std::ofstream bashFile(folderName + "/create_graphs.sh"s);

    /* Explanation of combined bash command
     * 1. Concatenate streams of separately layouted graphviz graphs, in order
     * 2. Pipe that into gvpack, which combines the graphs into an ordered
     *    array (by which is read in first), into columns with only one row
     *    (redirect stderr to /dev/null, since the warning of node names being
     *    adapted is unneeded)
     * 3. Create the final layout of the combined graph using neato, preserving
     *    position attributes from gvpack (-n2) into an SVG file
     * 4. Write to a iteration-compatibly delimited (i.e. step 0 is 000)
     *    filename
     */

    for(const auto& iterPair : numGraphsMap) {
      const auto& step = iterPair.first;
      const auto& highestGraphIndex = iterPair.second;

      bashFile << "cat ";

      for(unsigned i = 0; i <= highestGraphIndex; ++i) {
        bashFile << "<(dot gor1-explainer-" << step << "-" << i << ".dot) ";
      }

      bashFile << "| gvpack -array_uc1 2>/dev/null "
        << "| neato -n2 -Tsvg "
        << "> gor1-explainer-"
        << std::setw(3) << std::setfill('0') << step
        << ".svg" << nl;

    }

    bashFile.close();
  }

  return 0;
}
