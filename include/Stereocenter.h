#ifndef INCLUDE_STEREOCENTERS_H
#define INCLUDE_STEREOCENTERS_H

#include <algorithm>
#include <memory>
#include <boost/optional.hpp>

#include "RankingInformation.h"

// Detection algorithm headers
#include "Delib/ElementTypeCollection.h"
#include "Delib/PositionCollection.h"

/*! @file
 *
 * Contains the abstract base class / interface for all classes that model
 * conformational isomery on the graph level.
 */

/* TODO
 * - documentation
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

//! Classes that store and manipulate steric information not intrinsic to the graph
namespace Stereocenters {

enum class ChiralityConstraintTarget {
  Positive,
  Flat,
  Negative
};

struct DihedralLimits {
  const std::array<AtomIndexType, 4> indices;
  const double lower, upper;

  DihedralLimits(
    const std::array<AtomIndexType, 4>& indices,
    const std::pair<double, double>& limits
  ) : indices(indices), 
      lower(limits.first),
      upper(limits.second) 
  {
    assert(lower < upper);
  }
};

struct ChiralityConstraintPrototype {
  const std::array<AtomIndexType, 4> indices;
  const ChiralityConstraintTarget target; 

  ChiralityConstraintPrototype(
    const std::array<AtomIndexType, 4>& indices,
    const ChiralityConstraintTarget& target
  ) : indices(indices),
      target(target)
  {}
};

enum class Type {
  CNStereocenter,
  EZStereocenter
};


class Stereocenter {
public:
/* Modification */
  //!  Assign this feature
  virtual void assign(const boost::optional<unsigned>& assignment) = 0;

  //! Update vertex descriptors on vertex removal
  virtual void propagateVertexRemoval(const AtomIndexType& removedIndex) = 0;

/* Information */
  /*!
   * Return the angle imposed by the underlying symmetry defined by three
   * involved atoms in degrees. It needs to be a ternary function in order for
   * the angle requested to be clearly defined in both CNStereocenters and
   * EZStereocenters, in the latter of which two central atoms (that could be
   * j) exist.
   */
  virtual double angle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const = 0;

  /*!
   * Return whether this Stereocenter has been assigned or not
   * -> This leads to different behavior in DG! If unassigned, an Assignment is 
   *    chosen at random and adhered to during coordinate generation.
   */
  virtual boost::optional<unsigned> assigned() const = 0;

  //!  Return the number of possible assignments 
  virtual unsigned numAssignments() const = 0;

  //!  Return a list of chirality constraints
  virtual std::vector<ChiralityConstraintPrototype> chiralityConstraints() const = 0;

  /*!
   * Return the dihedral angle limits imposed by the underlying symmetries and
   * the current assignment.
   */
  virtual std::vector<DihedralLimits> dihedralLimits() const = 0;

  //!  Return a string giving information about the stereocenter
  virtual std::string info() const = 0;

  /*!
   * Return the set of center atoms (atoms that angle information is available 
   * on if asked as the central atom of an angle).
   */
  virtual std::vector<AtomIndexType> involvedAtoms() const = 0;

  //! Return the Subtype of the Stereocenter
  virtual Type type() const = 0;

/* Operators */
  //!  Ostream operator for debugging
  friend std::basic_ostream<char>& MoleculeManip::operator << (
    std::basic_ostream<char>& os,
    const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
  );
};

bool compareStereocenterEqual(
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& b
);

bool compareStereocenterLessThan(
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& b
);

} // namespace Stereocenters

} // namespace MoleculeManip

#endif
