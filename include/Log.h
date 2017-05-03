#ifndef INCLUDE_LOG_H
#define INCLUDE_LOG_H

#include <set>
#include <iostream>

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

// Level of logging
enum class Level {
  Trace,
  Debug,
  Info,
  Warning,
  Error,
  Fatal,
  None
};

// Particular cases of special things that may or may not be desired
enum class Particulars {
  /* In AdjacencyList.cpp, where a StereocenterFit is performed for all 
   * potential Stereocenters when a Molecule is read in, you can have numerical
   * details of the fit logged. Corresponding analysis scripts also exist.
   */
  StereocenterFitAnalysisInfo,
  /* In generateConformation.cpp, when chirality constraint prototypes are
   * fully determined into chirality constraints, emit some debug information
   */
  PrototypePropagatorDebugInfo,
  // In DGRefinementProblem, chirality constraint numerical debug information
  DGRefinementChiralityNumericalDebugInfo,
  /* In DGRefinementProblem, the callback function can reveal some information
   * on the current status of the optimization
   */
  DGRefinementProgress,
  // In generateConformation, show the Trees generated from the molecules
  gatherDGInformationTrees,
  // In BFSConstraintCollector, show the node index operator() is called on
  BFSConstraintCollectorVisitCall,
  // in debugDistanceGeometry, progress information
  DGDebugInfo,
};


// Log variables
extern Level level;
extern std::set<Particulars> particulars;

// Logging calls
std::ostream& log(const Level& decisionLevel);
std::ostream& log(const Particulars& particular);

} // namespace Log

#endif
