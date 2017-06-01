Overall general points
----------------------

- Using strict SSSR is disadvised, RDKit uses a custom "symmetrized" SSSR which
  isn't stricly smallest, but still small
- Newer concepts such as URF or Vismara's algorithm seem (judging using criteria
  from math paper) better for a new implementation
- URF exists as well-tested library and authors claim implementation is
  difficult. Haven't read URF paper properly yet though
- Generally, BFS to figure out size of smallest cycle selected atom is in should
  be easy and quick to implement, no need of advanced ring perception algorithms
- Nevertheless, having a reliable and reasonably future-proof method implemented
  would probably ease things like aromaticity detection (which I would use just
  for enforcing coplanarity)
- Overall, not too hot on imitating RDKit. Much prefer to use URF
- URF needs Vismara's algorithm


RDKIT Figueras' algorithm
-------------------------

(typically, removing / trimming means removing all bonds to and from the vertex,
not deleting the vertex itself, because this can lead to bookkeeping issues)

overall pseudocode::

  fullSet = AdjacencyList<AtomIndexType>(mol)
  SSSR = set< set<AtomIndexType> > {}
  do
    // iteratively remove all atoms with degree 1 (typically exposes more in
    // every iteration)
    while(degreeOneAtoms = fullSet.which(degree=1); degreeOneAtoms > 0) {
      fullSet.remove(degreeOneAtoms)
    }

    [lowestDegree, lowestDegreeVertices] = fullSet.lowestDegreeVertices()
    if(lowestDegree == 2) {
      // group vertices by connectivity
      chains = fullSet.makeChains(lowestDegreeVertices)

      for(chain : chains) {
        // Pick just one D2 for this chain, BFS for all D2s in chain leads to
        // same cycle
        i = chain.front

        ringSet = fullSet.getRing(i)
        if(ringSet.size() > 0) {
          if(SSSR.count(ringSet) == 0) {
            SSSR.insert(ringSet)
          }
        }
      }

      // Remove one D2 from every chain, all others degree is reduced and
      // cleaned out by removal algorithm at front of loop
      for(chain : chains) {
        fullSet.remove(chain.front)
      }
    } else if(lowestDegree == 3) {
      // Pick a random D3, fetch cycle

      i = pickRandom(lowestDegreeVertices)
      ringSets = fullSet.getThreeRings(i)
      for(ringSet : ringSets) {
        if(ringSet.size() > 0) {
          SSSR.insert(ringSet)
        }
      }
    }
  while (fullSet.size() > 0)

  // For SSSR, need to remove excess cycles from D3 nodes if present,
  // alternatively, we would now have a small (not smallest) set of smallest
  // rings
  if(SSSR.size() > fullSet.nullity()) {
    // remove largest cycles until the size fits
  }


checkEdges::
  // Selects an optimum edge for elimination in structures without D2s
  function checkEdges(n, ringSet, fullSet) {
    for(edge : ringSet.edges) {
      trialSet = fullSet
      trialSet.removeEdge(edge)
    }
  }

