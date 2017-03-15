#ifndef INCLUDE_STEREOCENTERS_H
#define INCLUDE_STEREOCENTERS_H

#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include <boost/optional.hpp>

#include "common_typedefs.h"
#include "StdlibTypeAlgorithms.h"
#include "AdjacencyList.h"
#include "DistanceGeometry/DistanceGeometry.h"

#include "Types/ElementTypeCollection.h" // Delib

/* TODO
 * change the abstract base class to include self-detection algorithms in 
 * Molecules. In order to directly access the private members of Molecule, all
 * GraphFeature derived classes must be friends of the Molecule class.
 */

namespace MoleculeManip {

// Predeclaration
namespace Stereocenters {
class Stereocenter;
}

std::basic_ostream<char>& operator << (
  std::basic_ostream<char>& os,
  const std::shared_ptr<Stereocenters::Stereocenter>& stereocenterPtr
);

namespace Stereocenters {

class Stereocenter {
public:
/* Public member functions */
  /* Modification */
  /*!
   * Assign this feature
   */
  virtual void assign(const unsigned& assignment) = 0;

  /* Information */
  /*!
   * Return a string specifying the type of stereocenter
   */
  virtual std::string type() const = 0;

  /*!
   * Return a set of involved atom indices
   */
  virtual std::set<AtomIndexType> involvedAtoms() const = 0;

  //!  Return a list of distance constraints
  virtual std::vector<
    DistanceGeometry::DistanceConstraint
  > distanceConstraints() const = 0;

  //!  Return a list of chirality constraints
  virtual std::vector<
    DistanceGeometry::ChiralityConstraint
  > chiralityConstraints() const = 0;

  /*!
   * Return the list of possible assignments at this feature
   */
  virtual unsigned assignments() const = 0;

  /*!
   * Return whether this feature has been assigned or not
   */
  virtual boost::optional<unsigned> assigned() const = 0;

  friend std::basic_ostream<char>& MoleculeManip::operator << (
    std::basic_ostream<char>& os,
    const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
  );
};

} // eo namespace Stereocenters

} // eo namespace

#endif
