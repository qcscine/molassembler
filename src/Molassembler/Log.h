/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Basic Logging functionality for debugging
 */

#ifndef INCLUDE_MOLASSEMBLER_LOG_H
#define INCLUDE_MOLASSEMBLER_LOG_H

#include "Molassembler/Export.h"
#include <unordered_set>
#include <iostream>

namespace Scine {
namespace Molassembler {
namespace Log {

namespace Detail {
class NullBuffer : public std::streambuf {
public:
  int overflow(int c) override;
};

// Some objects we need
extern NullBuffer nullBuffer;
extern std::ostream nullStream;
} // namespace Detail

//! Level of logging
enum class MASM_EXPORT Level : unsigned {
  Trace,
  Debug,
  Info,
  Warning,
  Error,
  Fatal,
  None
};

//! Particular cases of special logging items that may or may not be desired
enum class MASM_EXPORT Particulars {
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


//! Library logging level
MASM_EXPORT extern Level level;
//! Library logging particulars
MASM_EXPORT extern std::unordered_set<Particulars> particulars;

/**
 * @brief Fetch a log handle with a logging level
 * @param decisionLevel logging level of a message to write to a stream
 * @return std::cout if level is greater or equal to the library logging level,
 *   a null-stream otherwise
 */
std::ostream& log(const Level& decisionLevel);
/**
 * @brief Fetch a log handle with a particular
 * @param particular The particular to which the message pertains
 * @return std::cout if the particular is part of the current library
 *   particulars, a null-stream otherwise
 */
std::ostream& log(const Particulars& particular);
//! Checks whether a particular is part of the current library particulars
bool isSet(Particulars particular);

} // namespace Log
} // namespace Molassembler
} // namespace Scine
#endif
