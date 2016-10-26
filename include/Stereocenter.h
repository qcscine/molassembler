#ifndef INCLUDE_GRAPH_FEATURES_H
#define INCLUDE_GRAPH_FEATURES_H

#include <vector>
#include <set>
#include <algorithm>
#include <memory>

#include "common_typedefs.h"

// Detection algorithm headers
#include "AdjacencyList.h"
#include "EdgeList.h"
#include "Types/ElementTypeCollection.h" // Delib

/* TODO
 * change the abstract base class to include self-detection algorithms in 
 * Molecules. In order to directly access the private members of Molecule, all
 * GraphFeature derived classes must be friends of the Molecule class.
 */

namespace MoleculeManip {

namespace Stereocenters {

// TODO temp names
using Assignment = unsigned;

class Stereocenter {
public:
/* Public member functions */
  /* Modification */
  /*!
   * Assign this feature
   */
  virtual void assign(const Assignment& assignment) = 0;

  /* Information */
  /*!
   * Return a string specifying the type of stereocenter
   */
  virtual std::string type() const = 0;
  /*!
   * Return a set of involved atom indices
   */
  virtual std::set<AtomIndexType> involvedAtoms() const = 0;
  /*!
   * Return a list of distance constraints
   */
  virtual std::vector<DistanceConstraint> distanceConstraints() const = 0;
  /*!
   * Return a list of chirality constraints
   */
  virtual std::vector<ChiralityConstraint> chiralityConstraints() const = 0;
  /*!
   * Return the list of possible assignments at this feature
   */
  virtual std::vector<Assignment> assignments() const = 0;
  /*!
   * Return whether this feature has been assigned or not
   */
  virtual bool assigned() const = 0;
};

} // eo namespace GraphFeatures

} // eo namespace

#endif
