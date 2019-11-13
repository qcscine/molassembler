/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Point group symmetry elements
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_POINT_GROUP_ELEMENTS_H
#define INCLUDE_MOLASSEMBLER_SHAPES_POINT_GROUP_ELEMENTS_H

#include "shapes/PointGroups.h"

#include <Eigen/Core>
#include "boost/optional/optional_fwd.hpp"
#include "temple/Preprocessor.h"

#include <vector>
#include <unordered_map>
#include <memory>

namespace Scine {
namespace Shapes {
namespace elements {

//! Base class for symmetry elements
struct SymmetryElement {
  using Vector = Eigen::Vector3d;
  using Matrix = Eigen::Matrix3d;

  virtual ~SymmetryElement() = default;

  //! Returns a matrix representation of the symmetry element operation
  virtual Matrix matrix() const = 0;
  //! Returns a spatial vector unaffected by the symmetry element, if that exists
  virtual boost::optional<Vector> vector() const = 0;
  //! Returns a debug string representation of the element
  virtual std::string name() const = 0;
};

//! E symmetry element
struct Identity final : public SymmetryElement {
  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;
};

//! i symmetry element
struct Inversion final : public SymmetryElement {
  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;
};

//! Abstraction of Cn and Sn symmetry elements
struct Rotation final : public SymmetryElement {
  Rotation(
    const Eigen::Vector3d& passAxis,
    const unsigned passN,
    const unsigned passPower,
    const bool passReflect
  );

  //! Construct a proper rotation symmetry element
  static Rotation Cn(const Eigen::Vector3d& axis, const unsigned n, const unsigned power = 1);
  //! Construct an improper rotation symmetry element
  static Rotation Sn(const Eigen::Vector3d& axis, const unsigned n, const unsigned power = 1);

  //! Multiply this rotation element by another
  Rotation operator * (const Rotation& rhs) const;

  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;

  //! Axis of rotation
  Eigen::Vector3d axis;
  //! Rotation order
  unsigned n;
  //! Power of the symmetry element (i.e. C3^2 is a rotation by 240°)
  unsigned power;
  //! Proper rotations reflect, improper rotations do not
  bool reflect;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

//! Reflection by a plane symmetry element
struct Reflection final : public SymmetryElement {
  Reflection(const Eigen::Vector3d& passNormal);

  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;

  //! Normal of the reflection plane
  Eigen::Vector3d normal;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

//! Heterogeneous list of symmetry elements
using ElementsList = std::vector<std::unique_ptr<SymmetryElement>>;

/** @brief Lists all symmetry elements for a point group
 *
 * @param group Point group for which to enumerate symmetry elements
 *
 * @note Yields an empty list of elements for invalid @p group values.
 *
 * @warning Cinfv and Dinfh return the elements for C8v and D8h, respectively.
 *
 * @warning Ih and I are not implemented yet. Will assert in Debug and yield a
 * single-element list in Release.
 *
 * @return a list of symmetry elements
 */
PURITY_WEAK ElementsList symmetryElements(const PointGroup group) noexcept;

unsigned order(const PointGroup group);


struct ElementGrouping {
  using ElementIndexGroups = std::vector<std::vector<unsigned>>;

  Eigen::Vector3d probePoint;
  ElementIndexGroups groups;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * Map from number of points to element groups, including the probe point
 * E.g. there is a grouping of symmetry elements for four points in Td (G = 24)
 * So the resulting map can be indexed with four to get the symmetry element
 * groups.
 */
using NPGroupingsMapType = std::unordered_map<
  unsigned,
  std::vector<ElementGrouping>,
  std::hash<unsigned>,
  std::equal_to<unsigned>,
  Eigen::aligned_allocator<std::pair<const unsigned, std::vector<ElementGrouping>>>
>;

/**
 * @brief Generate all groupings of symmetry elements for those points in space
 * for which some symmetry elements have the same effect
 *
 * @param elements Symmetry elements list
 *
 * @return Groups of symmetry elements
 */
NPGroupingsMapType npGroupings(const ElementsList& elements);

} // namespace elements
} // namespace Shapes
} // namespace Scine

#endif
