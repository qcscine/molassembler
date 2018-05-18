#ifndef INCLUDE_STEREOCENTERS_H
#define INCLUDE_STEREOCENTERS_H

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

namespace molassembler {

/* Predeclarations */
class Molecule;

namespace DistanceGeometry {
struct ChiralityConstraint;
class MoleculeSpatialModel;
} // namespace DistanceGeometry

//! Classes that store and manipulate steric information not intrinsic to the graph
namespace Stereocenters {

enum class Type {
  CNStereocenter,
  EZStereocenter
};

class Stereocenter {
public:
/* Virtual destructor */
  virtual ~Stereocenter() = default;

/* Modification */
  //!  Assign this feature
  virtual void assign(const boost::optional<unsigned>& assignment) = 0;

  //! Assign this feature at random
  virtual void assignRandom() = 0;

  //! Update vertex descriptors on vertex removal
  virtual void propagateVertexRemoval(const AtomIndexType removedIndex) = 0;

/* Information */
  /*!
   * Return whether this Stereocenter has been assigned or not
   * -> This leads to different behavior in DG! If unassigned, an Stereopermutation is
   *    chosen at random and adhered to during coordinate generation.
   */
  virtual boost::optional<unsigned> assigned() const = 0;

  //!  Return the number of possible assignments
  virtual unsigned numStereopermutations() const = 0;

  //!  Return a list of chirality constraints
  virtual std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints() const = 0;

  //!  Return a string giving information about the stereocenter
  virtual std::string info() const = 0;

  /*!
   * Return the set of center atoms (atoms that angle information is available
   * on if asked as the central atom of an angle).
   */
  virtual std::vector<AtomIndexType> involvedAtoms() const = 0;

  virtual void setModelInformation(
    DistanceGeometry::MoleculeSpatialModel& model,
    const std::function<double(const AtomIndexType)> cycleMultiplierForIndex,
    const double looseningMultiplier
  ) const = 0;

  //! Return the Subtype of the Stereocenter
  virtual Type type() const = 0;
};

bool compareStereocenterEqual(
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& b
);

bool compareStereocenterLessThan(
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& a,
  const std::shared_ptr<molassembler::Stereocenters::Stereocenter>& b
);

} // namespace Stereocenters

} // namespace molassembler

#endif
