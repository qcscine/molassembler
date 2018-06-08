#define BOOST_TEST_MODULE EZStereocenterTestModule
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EZStereocenter.h"

BOOST_AUTO_TEST_CASE(stateConsistency) {
  using namespace molassembler;
  using namespace molassembler::Stereocenters;

  RankingInformation left, right;
  left.sortedSubstituents = {{0}, {1}};
  left.ligands = {{0}, {1}};
  left.ligandsRanking = {{0}, {1}};
  right.sortedSubstituents = {{4}, {5}};
  right.ligands = {{4}, {5}};
  right.ligandsRanking = {{0}, {1}};

  EZStereocenter trialStereocenter {
    2,
    left,
    3,
    right
  };

  trialStereocenter.assign(0u);

  // Update, where right gets flipped
  RankingInformation rightFlipped;
  rightFlipped.sortedSubstituents = {{5}, {4}};
  rightFlipped.ligands = {{4}, {5}};
  rightFlipped.ligandsRanking = {{1}, {0}};
  trialStereocenter.propagateGraphChange(left, rightFlipped);

  BOOST_CHECK_MESSAGE(
    trialStereocenter.assigned() == 1u,
    "Graph update with flipped priorities does not yield correct stereoconfiguration!"
  );

  // Notify of removal, then notify of substituent removal
  trialStereocenter.propagateVertexRemoval(0);
  trialStereocenter.removeSubstituent(
    1,
    Stereocenter::removalPlaceholder
  );

  BOOST_CHECK_MESSAGE(
    trialStereocenter.numAssignments() == 2
    && trialStereocenter.assigned() == 1u,
    "Vertex removal does not preserve chiral information!"
    << " numAssignments: " << trialStereocenter.numAssignments()
  );
}
