/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereopermutatorList.h"

#include "Utils/ElementTypes.h"

#include "boost/process/child.hpp"
#include "boost/process/io.hpp"
#include "boost/process/search_path.hpp"

bool graphvizInPath() {
  return !boost::process::search_path("dot").empty();
}

std::string pipeSVG(const Scine::molassembler::Molecule& molecule) {
  std::string callString = "dot -Tsvg";
  std::stringstream os;

  // Construct pipe streams for redirection
  boost::process::opstream ips;
  boost::process::pstream ps;
  boost::process::pstream err;

  // Start the child process
  boost::process::child childProcess(callString, boost::process::std_in<ips, boost::process::std_out> ps,
                                     boost::process::std_err > err);

  // Feed our graphviz into the process
  ips << molecule.dumpGraphviz();
  ips.flush();
  ips.pipe().close();

  // Wait for the child process to exit
  childProcess.wait();

  std::stringstream stderrStream;
#if BOOST_VERSION >= 107000
  /* NOTE: This implementation of buffer transfers in boost process has a bug
   * that isn't fixed before Boost 1.70.
   */
  os << ps.rdbuf();
  stderrStream << err.rdbuf();
#else
  // Workaround: cast to a parent class implementing rdbuf() correctly.
  using BasicIOSReference = std::basic_ios<char, std::char_traits<char>>&;
  // Feed the results into our ostream
  os << static_cast<BasicIOSReference>(ps).rdbuf();
  stderrStream << static_cast<BasicIOSReference>(err).rdbuf();
#endif

  return os.str();
}

bool modularCompare(
  const Scine::molassembler::Molecule& molecule,
  const Scine::molassembler::Molecule& other,
  const bool compareElementTypes,
  const bool compareBondOrders,
  const bool compareSymmetries,
  const bool compareStereopermutations
) {
  using Components = Scine::molassembler::AtomEnvironmentComponents;
  Components bitmask;

  if(compareElementTypes) {
    bitmask = bitmask | Components::ElementTypes;
  }

  if(compareBondOrders) {
    bitmask = bitmask | Components::BondOrders;
  }

  if(compareSymmetries) {
    bitmask = bitmask | Components::Symmetries;
  }

  if(compareStereopermutations) {
    bitmask = bitmask | Components::Stereopermutations;
  }

  return molecule.modularCompare(other, bitmask);
}

void canonicalizeMolecule(Scine::molassembler::Molecule& molecule) {
  molecule.canonicalize();
}

void init_molecule(pybind11::module& m) {
  using namespace Scine::molassembler;
  using namespace Scine::Utils;

  pybind11::class_<Molecule> molecule(
    m,
    "Molecule",
    "Models a molecule as a graph and a list of stereopermutators"
  );

  /* Constructors */
  molecule.def(
    pybind11::init<>(),
    "Initialize a hydrogen molecule"
  );

  molecule.def(
    pybind11::init<ElementType, ElementType, BondType>(),
    "Initialize a molecule from two element types and a mutual bond type"
  );

  molecule.def(
    pybind11::init<OuterGraph>(),
    "Initialize a molecule from connectivity alone, inferring symmetries and "
    "stereopermutators from the graph"
  );

  /* Modifiers */
  molecule.def(
    "add_atom",
    &Molecule::addAtom,
    "Add an atom to the molecule."
  );

  molecule.def(
    "add_bond",
    &Molecule::addBond,
    "Adds a bond to the molecule."
  );

  molecule.def(
    "assign_stereopermutator",
    pybind11::overload_cast<AtomIndex, const boost::optional<unsigned>&>(
      &Molecule::assignStereopermutator
    ),
    "Sets the stereopermutator at a particular atom"
  );

  molecule.def(
    "assign_stereopermutator",
    pybind11::overload_cast<const BondIndex&, const boost::optional<unsigned>&>(
      &Molecule::assignStereopermutator
    ),
    "Sets the stereopermutator at a particular bond"
  );

  molecule.def(
    "assign_stereopermutator_randomly",
    pybind11::overload_cast<AtomIndex>(
      &Molecule::assignStereopermutatorRandomly
    ),
    "Assigns an atom stereopermutator at random"
  );

  molecule.def(
    "assign_stereopermutator_randomly",
    pybind11::overload_cast<const BondIndex&>(
      &Molecule::assignStereopermutatorRandomly
    ),
    "Assigns a bond stereopermutator at random"
  );

  // TODO this is a temporary fix since being able to canonicalize some components will be necessary for the isomerism checks
  // also, access to the parts of a molecule that are canonical is needed
  molecule.def(
    "canonicalize",
    &::canonicalizeMolecule,
    "Canonicalizes the molecular graph"
  );

  molecule.def(
    "remove_atom",
    &Molecule::removeAtom,
    "Remove an atom from the graph, including bonds to it"
  );

  molecule.def(
    "remove_bond",
    &Molecule::removeBond,
    "Remove a bond from the graph"
  );

  molecule.def(
    "set_bond_type",
    &Molecule::setBondType,
    "Change the bond type of two atoms."
  );

  molecule.def(
    "set_element_type",
    &Molecule::setElementType,
    "Change the element type of an atom"
  );

  molecule.def(
    "set_geometry_at_atom",
    &Molecule::setGeometryAtAtom,
    "Change the local geometry at an atom"
  );

  /* Information */
  molecule.def(
    "dump_graphviz",
    &Molecule::dumpGraphviz,
    "Returns a graphviz string representation of the molecule"
  );

  molecule.def_property_readonly(
    "graph",
    &Molecule::graph,
    "Fetches read only access to the graph representation"
  );

  molecule.def_property_readonly(
    "stereopermutators",
    &Molecule::stereopermutators,
    "Fetches read only access to the list of stereopermutators"
  );

  molecule.def(
    "partial_compare",
    &::modularCompare,
    pybind11::arg("other"),
    pybind11::arg("compare_element_types") = true,
    pybind11::arg("compare_bond_orders") = true,
    pybind11::arg("compare_symmetries") = true,
    pybind11::arg("compare_stereopermutations") = true
  );

  molecule.def(pybind11::self == pybind11::self);
  molecule.def(pybind11::self != pybind11::self);

  /* Integration with IPython / Jupyter */
  molecule.def(
    "_repr_svg_",
    &::pipeSVG,
    "Generates an SVG representation of the molecule"
  );
}
