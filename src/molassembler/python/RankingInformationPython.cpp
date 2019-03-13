/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

#include "molassembler/RankingInformation.h"

void init_ranking_information(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<LinkInformation> linkInformation(
    m,
    "LinkInformation",
    "Information on links (graph paths) between substituents of a central atom"
  );

  linkInformation.def_readonly(
    "index_pair",
    &LinkInformation::indexPair,
    "An ordered pair of the ligand indices that are linked. See the "
    "corresponding RankingInformation's ligands member"
  );

  linkInformation.def_readonly(
    "cycle_sequence",
    &LinkInformation::cycleSequence,
    "The in-order atom sequence of the linked ligand indices"
  );

  linkInformation.def(pybind11::self == pybind11::self);
  linkInformation.def(pybind11::self != pybind11::self);
  linkInformation.def(pybind11::self < pybind11::self);


  pybind11::class_<RankingInformation> rankingInformation(
    m,
    "RankingInformation",
    "Ranking data of substituents around a central vertex"
  );

  rankingInformation.def_readonly(
    "ranked_substituents",
    &RankingInformation::substituentRanking,
    "Sorted substituents grouped by ascending priority"
  );

  rankingInformation.def_readonly(
    "sites",
    &RankingInformation::sites,
    "An unordered nested list of atom indices that constitute a binding site"
  );

  rankingInformation.def_readonly(
    "ranked_sites",
    &RankingInformation::siteRanking,
    "An ordered nested list of indices into the sites member"
  );

  rankingInformation.def_readonly(
    "links",
    &RankingInformation::links,
    "An ordered list of information on all links between ligand sites"
  );

  rankingInformation.def(pybind11::self == pybind11::self);
  rankingInformation.def(pybind11::self != pybind11::self);
}
