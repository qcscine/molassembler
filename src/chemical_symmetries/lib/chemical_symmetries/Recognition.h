/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Analyze coordinates for point group symmetry
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_RECOGNITION_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_RECOGNITION_H

#include <Eigen/Core>
#include <vector>
#include <memory>
#include <map>
#include "boost/optional.hpp"

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

using PositionCollection = Eigen::Matrix<double, 3, Eigen::Dynamic>;
PointGroup analyze(const PositionCollection& positions);

// TODO TMP for testing
namespace minimization {

struct SymmetryElement {
  using Vector = Eigen::Vector3d;
  using Matrix = Eigen::Matrix3d;

  virtual ~SymmetryElement() = default;
  virtual Matrix matrix() const = 0;
  virtual boost::optional<Vector> vector() const = 0;
};

struct Identity final : public SymmetryElement {
  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
};

struct Inversion final : public SymmetryElement {
  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;
};

struct Rotation : public SymmetryElement {
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

  Eigen::Vector3d axis;
  unsigned n;
  unsigned power;
  bool reflect;
};

struct Reflection : public SymmetryElement {
  Reflection(const Eigen::Vector3d& passNormal);

  Matrix matrix() const final;
  boost::optional<Vector> vector() const final;

  Eigen::Vector3d normal;
};

std::vector<std::unique_ptr<SymmetryElement>> symmetryElements(const PointGroup group);


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
};

std::map<unsigned, ElementGrouping> npGroupings(const PointGroup group);

} // namespace minimization

} // namespace Symmetry
} // namespace Scine

#endif
