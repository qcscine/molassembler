#include <vector>
#include <set>
#include <algorithm>
#include <memory>

#include "common_typedefs.h"

namespace MoleculeManip {
// TODO temp names
using DistanceConstraint = std::tuple<
  AtomIndexType, // i
  AtomIndexType, // j
  double, // lower
  double // uper
>;
using ChiralityConstraint = std::tuple<
  AtomIndexType, // i
  AtomIndexType, // j
  AtomIndexType, // k
  AtomIndexType, // l
  double // target
>;
using Assignment = unsigned;

class GraphFeature {
public:
/* Public member functions */
  /* Modification */
  /*!
   * Assign this feature
   */
  virtual void assign(const Assignment& assignment) = 0;
  /* Information */
  /*!
   * Return a string specifying the type of feature 
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
  /*
   * Return the number of possible isomers stemming from this feature
   */
  //virtual unsigned numIsomers() const = 0;
  /*!
   * Return the list of possible assignments at this feature
   */
  virtual std::vector<Assignment> assignments() const = 0;
  /*!
   * Return whether this GraphFeature has Assignments
   */
  virtual bool hasAssignments() const = 0;
  /*!
   * Return whether this feature has been assigned or not
   */
  virtual bool assigned() const = 0;
};


} // eo namespace
