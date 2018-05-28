#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "boost/regex.hpp"

#include "temple/Containers.h"

#include "DistanceGeometry/ImplicitGraphBoost.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "IO.h"
#include "StdlibTypeAlgorithms.h"

#include "gor1/Gor1.h"

#include "boost/graph/graphviz.hpp"

#include <random>

using namespace molassembler;

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

/* TODO
 * - Alter to specialized SPG Gor1 variant
 */

namespace molassembler {

namespace DistanceGeometry {

struct WriterState {
  using Graph = ImplicitGraph;
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

  boost::optional<Edge> highlightEdge = boost::none;
  boost::optional<Vertex> highlightVertex = boost::none;
};

class GraphvizWriter {
public:
  using Graph = ImplicitGraph;
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<Graph>::edge_descriptor;
  using ColorMapType = boost::two_bit_color_map<>;
  using Color = boost::color_traits<typename boost::property_traits<ColorMapType>::value_type>;

private:
  const Graph& _graph;
  const std::vector<double>& _distances;
  const std::vector<Vertex>& _predecessors;
  const ColorMapType& _colors;

public:
  GraphvizWriter(
    const Graph& graph,
    const std::vector<double>& distances,
    const std::vector<Vertex>& predecessors,
    const ColorMapType& colors
  ) : _graph {graph},
      _distances {distances},
      _predecessors {predecessors},
      _colors {colors}
  {}

  std::shared_ptr<WriterState> statePtr {new WriterState};

  void write_clusters(std::ostream& os) {
    os << "\nsubgraph cluster_left {\n";
    for(Vertex i = 0; i < boost::num_vertices(_graph); i += 2) {
      os << "  " << i << "; ";
    }
    os << "\n}\nsubgraph cluster_right {\n";
    for(Vertex i = 1; i < boost::num_vertices(_graph); i += 2) {
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
    return _predecessors.at(v) != v;
  }

  void operator() (std::ostream& os, const Vertex& v) {
    os << "[";


    auto v_color = boost::get(_colors, v);
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

    if(_distances.at(v) < 10000) {
      os << _distances.at(v);
    } else {
      os << "inf";
    }
    os << "\"";

    os << "]";

    if(boost::num_vertices(_graph) <= 20 && v == boost::num_vertices(_graph) - 1) {
      write_clusters(os);
    }
  }

  bool inShortestPathsTree(const Vertex& source, const Vertex& target) {
    return (
      _predecessors.at(target) == source
      || (
        target % 2 == source % 2
        && _predecessors.at(source) == target
      )
    );
  }

  void operator() (std::ostream& os, const Edge& e) {
    os << "[";

    os << R"(label=")" << boost::get(boost::edge_weight, _graph)(e) << R"(")";

    auto source = boost::source(e, _graph);
    auto target = boost::target(e, _graph);

    bool isHighlightedEdge = false;
    if(statePtr->highlightEdge) {
      auto hSource = boost::source(statePtr->highlightEdge.value(), _graph);
      auto hTarget = boost::target(statePtr->highlightEdge.value(), _graph);

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

    if(source % 2 != target % 2 && !_graph.hasExplicit(e)) {
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
  using Graph = ImplicitGraph;
  using Vertex = typename boost::graph_traits<ImplicitGraph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<ImplicitGraph>::edge_descriptor;

private:
  const std::vector<Vertex>& _predecessors;

  GraphvizWriter& _writer;
  unsigned _counter = 0;
  std::deque<Vertex> _A, _B;

public:
  VisualizationVisitor(
    const std::vector<Vertex>& predecessors,
    GraphvizWriter& writer
  ) : _predecessors {predecessors},
      _writer {writer}
  {}

  void examine_edge(const Edge& e, const Graph& g) {
    _writer.statePtr->highlightEdge = e;
    makeGraphvizFiles(g);
  }

  void relax_edge(const Edge& e, const Graph& g) {
    _writer.statePtr->highlightEdge = e;
    makeGraphvizFiles(g);
  }

  void a_push(const Vertex& v, const Graph& g) {
    _writer.statePtr->highlightVertex = v;
    _A.push_front(v);
    makeGraphvizFiles(g);
  }
  void a_pop(const Vertex& v, const Graph& g) {
    _writer.statePtr->highlightVertex = v;
    _A.pop_front();
    makeGraphvizFiles(g);
  }
  void b_push(const Vertex& v, const Graph& g) {
    _writer.statePtr->highlightVertex = v;
    _B.push_front(v);
    makeGraphvizFiles(g);
  }
  void b_pop(const Vertex& v, const Graph& g) {
    _writer.statePtr->highlightVertex = v;
    _B.pop_front();
    makeGraphvizFiles(g);
  }

  void mark_white(const Vertex&, const Graph& g) {
    makeGraphvizFiles(g);
  }
  void mark_gray(const Vertex&, const Graph& g) {
    makeGraphvizFiles(g);
  }
  void mark_black(const Vertex&, const Graph& g) {
    makeGraphvizFiles(g);
  }

  void makeGraphvizFiles(const Graph& g) {
    using namespace std::string_literals;

    const std::string prefix = "gor1-explainer-"s
      + std::to_string(_counter)
      + "-"s;

    std::ofstream graphFile(prefix + "0.dot"s);

    boost::write_graphviz(
      graphFile,
      g,
      _writer,
      _writer,
      _writer
    );

    graphFile.close();

    std::ofstream stackAFile(prefix + "1.dot"s);
    write_graphviz(stackAFile, _A, "A");
    stackAFile.close();

    std::ofstream stackBFile(prefix + "2.dot"s);
    write_graphviz(stackBFile, _B, "B");
    stackBFile.close();

    std::ofstream treeGraph(prefix + "3.dot"s);
    write_predecessor_graphviz(treeGraph, _predecessors);
    treeGraph.close();

    ++_counter;
  }
};

} // namespace DistanceGeometry

} // namespace molassembler

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
  if(options_variables_map.count("r")) {
    root_vertex = options_variables_map["r"].as<unsigned>();
  }

  if(options_variables_map.count("help")) {
    std::cout << options_description << nl;
    return 0;
  }

  if(options_variables_map.count("f")) {
    auto filename = options_variables_map["f"].as<std::string>();

    if(!boost::filesystem::exists(filename)) {
      std::cout << "The specified file could not be found!" << nl;
      return 1;
    }

    auto mol = IO::read(filename);

    DistanceGeometry::MoleculeSpatialModel spatialModel {mol};

    DistanceGeometry::ImplicitGraph shortestPathsGraph {
      mol,
      spatialModel.makeBounds()
    };

    /* Prep */
    using Vertex = typename boost::graph_traits<DistanceGeometry::ImplicitGraph>::vertex_descriptor;

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

      auto splat = StdlibTypeAlgorithms::split(iter->path().filename().string(), '-');
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
