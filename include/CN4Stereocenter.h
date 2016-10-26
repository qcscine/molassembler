#ifndef INCLUDE_CN4_STEREOCENTER_H
#define INCLUDE_CN4_STEREOCENTER_H

#include <math.h>
#include <experimental/optional>
#include "ElementTypes.h"

#include "Stereocenter.h"
#include "CommonTrig.h"
#include "Molecule.h"

namespace MoleculeManip {

namespace Stereocenters {

class CN4Stereocenter : public Stereocenter {
private:
  std::experimental::optional<Assignment> _assignment;

public:
/* Public member functions */
  /* Modification */
  /*!
   * Assign this feature
   */
  virtual void assign(const Assignment& assignment) override final {

  }

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
  virtual std::vector<Assignment> assignments() const {

  }
  /*!
   * Return whether this feature has been assigned or not
   */
  virtual bool assigned() const {
    // Fundamentals TS
    return (bool) _assignment;
    // C++17
    // return _assignment.has_value();
  }
};

}

}

#endif
