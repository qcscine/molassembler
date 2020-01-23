/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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
    R"delim(
      Information on links (graph paths) between substituents of a central atom

      This captures all cycles that the central atom whose substituents are
      being ranked and its sites are in.

      >>> # Simple example of links between substituents
      >>> import scine_utilities as utils
      >>> cyclopropane = io.experimental.from_smiles("C1CC1")
      >>> p = cyclopropane.stereopermutators.option(0)
      >>> # Sites are single-index, non-haptic
      >>> site_is_single_index = lambda s: len(s) == 1
      >>> all(map(site_is_single_index, p.ranking.sites))
      True
      >>> # There is a single link between carbon atom sites
      >>> is_carbon = lambda a: cyclopropane.graph.element_type(a) == utils.ElementType.C
      >>> site_is_carbon = lambda s: len(s) == 1 and is_carbon(s[0])
      >>> len(p.ranking.links) == 1
      True
      >>> single_link = p.ranking.links[0]
      >>> site_index_is_carbon = lambda s: site_is_carbon(p.ranking.sites[s])
      >>> all(map(site_index_is_carbon, single_link.index_pair))
      True
      >>> single_link.cycle_sequence # Atom indices of cycle members
      [0, 1, 2]
      >>> all(map(is_carbon, single_link.cycle_sequence)) # All carbons
      True
    )delim"
  );

  linkInformation.def_readonly(
    "index_pair",
    &LinkInformation::indexPair,
    "An ordered pair of the site indices that are linked. See the "
    "corresponding :class:`RankingInformation` sites member"
  );

  linkInformation.def_readonly(
    "cycle_sequence",
    &LinkInformation::cycleSequence,
    R"delim(
      The in-order atom sequence of the cycle involving the linked sites. The
      source vertex is always placed at the front of this sequence. The
      sequence is normalized such that second atom index is less than the last.
    )delim"
  );

  linkInformation.def(pybind11::self == pybind11::self);
  linkInformation.def(pybind11::self != pybind11::self);
  linkInformation.def(pybind11::self < pybind11::self);


  pybind11::class_<RankingInformation> rankingInformation(
    m,
    "RankingInformation",
    R"delim(
      Ranking data of substituents around a central vertex

      >>> # Model compound with a haptically bonded ethene
      >>> compound_smiles = "[Co]1(C#O)(C#O)(C#O)(C#O)(C#O)C=C1"
      >>> compound = io.experimental.from_smiles(compound_smiles)
      >>> cobalt_index = 0
      >>> p = compound.stereopermutators.option(cobalt_index)
      >>> is_haptic_site = lambda s: len(s) > 1
      >>> any(map(is_haptic_site, p.ranking.sites))
      True
      >>> # There are no links for this, none of the sites are interconnected
      >>> len(p.ranking.links)
      0
      >>> # All of the sites are ranked equally save for the haptic site
      >>> p.ranking.ranked_sites
      [[0, 1, 2, 3, 4], [5]]
      >>> p.ranking.sites[5] # The constituting atom indices of the haptic site
      [11, 12]
      >>> p.ranking.site_index_of_atom(12) # Look up atom indices
      5
      >>> p.ranking.rank_index_of_site(1) # Get ranking position of a site
      0
    )delim"
  );

  rankingInformation.def_readonly(
    "ranked_substituents",
    &RankingInformation::substituentRanking,
    "Sorted substituents grouped by ascending priority"
  );

  rankingInformation.def_readonly(
    "sites",
    &RankingInformation::sites,
    "An unordered nested list of atom indices that constitute binding sites"
  );

  rankingInformation.def_readonly(
    "ranked_sites",
    &RankingInformation::siteRanking,
    "An ordered nested list of indices into the sites member"
  );

  rankingInformation.def_readonly(
    "links",
    &RankingInformation::links,
    "An ordered list of :class:`LinkInformation` on all links between binding sites"
  );

  rankingInformation.def(
    "site_index_of_atom",
    &RankingInformation::getSiteIndexOf,
    pybind11::arg("atom_index"),
    "Fetch the site index of an atom index"
  );

  rankingInformation.def(
    "rank_index_of_site",
    &RankingInformation::getRankedIndexOfSite,
    pybind11::arg("site_index"),
    "Fetch the position of a site within the site ranking"
  );

  rankingInformation.def(pybind11::self == pybind11::self);
  rankingInformation.def(pybind11::self != pybind11::self);
}
