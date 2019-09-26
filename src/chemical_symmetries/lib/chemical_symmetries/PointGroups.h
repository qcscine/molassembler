/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Point group enum and symmetry elements
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_POINT_GROUPS_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_POINT_GROUPS_H

#include <Eigen/Core>
#include "boost/optional/optional_fwd.hpp"
#include <vector>
#include <memory>
#include <unordered_map>

#include "temple/Preprocessor.h"

namespace Scine {
namespace Symmetry {

/**
 * @brief Point groups
 */
enum class PointGroup : unsigned {
  C1, Ci, Cs,
  C2, C3, C4, C5, C6, C7, C8,
  C2h, C3h, C4h, C5h, C6h, C7h, C8h,
  C2v, C3v, C4v, C5v, C6v, C7v, C8v,
  S4, S6, S8,
  D2, D3, D4, D5, D6, D7, D8,
  D2h, D3h, D4h, D5h, D6h, D7h, D8h,
  D2d, D3d, D4d, D5d, D6d, D7d, D8d,
  T, Td, Th,
  O, Oh,
  I, Ih,
  Cinfv, Dinfh
};

namespace elements {

struct SymmetryElement {
  using Vector = Eigen::Vector3d;
  using Matrix = Eigen::Matrix3d;

  virtual ~SymmetryElement() = default;
  virtual Matrix matrix() const = 0;
  virtual boost::optional<Vector> vector() const = 0;
  virtual std::string name() const = 0;
};

struct Identity final : public SymmetryElement {
  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;
};

struct Inversion final : public SymmetryElement {
  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;
};

struct Rotation final : public SymmetryElement {
  Rotation(
    const Eigen::Vector3d& passAxis,
    const unsigned passN,
    const unsigned passPower,
    const bool passReflect
  );

  static Rotation Cn(const Eigen::Vector3d& axis, const unsigned n, const unsigned power = 1);
  static Rotation Sn(const Eigen::Vector3d& axis, const unsigned n, const unsigned power = 1);

  Rotation operator * (const Rotation& rhs) const;

  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;

  Eigen::Vector3d axis;
  unsigned n;
  unsigned power;
  bool reflect;
};

struct Reflection final : public SymmetryElement {
  Reflection(const Eigen::Vector3d& passNormal);

  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
  std::string name() const final;

  Eigen::Vector3d normal;
};

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


/* There can be multiple groupings of symmetry elements of equal l for
 * different points in space. For now, we are ASSUMING that a particular
 * grouping of symmetry elements leads the folded point to lie along the
 * axis defined by the probe point used here to determine the grouping.
 *
 * So we store the point suggested by the element along with its grouping
 * so we can test this theory.
 */
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
  Eigen::aligned_allocator<std::pair<unsigned, ElementGrouping>>
>;

NPGroupingsMapType npGroupings(const ElementsList& elements);

} // namespace elements


} // namespace Symmetry
} // namespace Scine

#endif
