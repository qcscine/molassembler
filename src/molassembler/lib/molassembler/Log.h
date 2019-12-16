/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Basic Logging functionality for debugging
 */

#ifndef INCLUDE_MOLASSEMBLER_LOG_H
#define INCLUDE_MOLASSEMBLER_LOG_H

#include <set>
#include <iostream>

namespace Scine {

namespace molassembler {

namespace Log {

namespace detail {
  class NullBuffer : public std::streambuf {
  public:
    int overflow(int c);
  };

  // Some objects we need
  extern NullBuffer nullBuffer;
  extern std::ostream nullStream;
}

//! Level of logging
enum class Level : unsigned {
  Trace,
  Debug,
  Info,
  Warning,
  Error,
  Fatal,
  None
};

//! Particular cases of special logging items that may or may not be desired
enum class Particulars {
  /*! In Molecule.cpp, where a fit of AtomStereopermutators against positions is
   * performed when a Molecule is read in, you can have numerical details of
   * the fit logged. Corresponding analysis scripts also exist.
   */
  AtomStereopermutatorFit,
  /*! AtomStereopermutator's addSubstituent, removeSubstituent and propagate* functions
   */
  AtomStereopermutatorStatePropagation,
  /*! In generateConformation.cpp, when chiral constraint prototypes are
   * fully determined into chiral constraints, emit some debug information
   */
  PrototypePropagatorDebugInfo,
  //! In DgRefinementProblem, chiral constraint numerical debug information
  DgRefinementChiralityNumericalDebugInfo,
  /*! In DgRefinementProblem, the callback function can reveal some information
   * on the current status of the optimization
   */
  DgRefinementProgress,
  //! In ConformerGeneration, explain final contributions to the error function
  DgFinalErrorContributions,
  //! Explain why a structure was not accepted
  DgStructureAcceptanceFailures,
  //! In generateConformation, show the Trees generated from the molecules
  gatherDGInformationTrees,
  //! in debugDistanceGeometry, progress information
  DgDebugInfo,
  //! Ranking debug information
  RankingTreeDebugInfo
};


// Log variables
extern Level level;
extern std::set<Particulars> particulars;

// Logging calls
std::ostream& log(const Level& decisionLevel);
std::ostream& log(const Particulars& particular);
bool isSet(Particulars particular);

} // namespace Log

} // namespace molassembler

} // namespace Scine
#endif
