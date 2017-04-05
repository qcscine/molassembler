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
  StereocenterFitAnalysisInfo 
};


// Log variables
extern Level level;
extern std::set<Particulars> particulars;

// Logging calls
std::ostream& log(const Level& decisionLevel);
std::ostream& log(const Particulars& particular);

}

#endif
