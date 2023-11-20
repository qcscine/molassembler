/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief For transformations between coordinate systems
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_COORDINATE_SYSTEM_TRANSFORMATION_H
#define INCLUDE_MOLASSEMBLER_SHAPES_COORDINATE_SYSTEM_TRANSFORMATION_H
#include <Eigen/Geometry>
#include <Eigen/Core>

#include "Molassembler/Export.h"

namespace Scine {
namespace Molassembler {
namespace Shapes {

/**
 * @brief Coordinate system axis data class
 */
struct MASM_EXPORT CoordinateSystem {
  //! Axes
  Eigen::Vector3d x, y, z;

  //! Default constructor makes default xyz system
  CoordinateSystem();
  /*! @brief Right-handed coordinate system constructor from two perpendicular vectors
   *
   * The z coordinate is initialized from the cross product of @p a and @p b
   *
   * @param a The x coordinate of the coordinate system
   * @param b The y coordinate of the coordinate system
   *
   * @pre @p a and @p b need to be orthogonal.
   * @post x, y and z are normalized, the coordinate system is right-handed
   */
  CoordinateSystem(const Eigen::Vector3d& a, const Eigen::Vector3d& b);

  //! Yields a randomly oriented right-handed coordinate system
  static CoordinateSystem random();

  //! Yields whether the coordinate system is right handed
  bool isRightHanded() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief Gives the rotation matrix tranforming points between coordinate systems
 *
 * You can transform points from @p a using the result of this function into
 * coordinate system @p b.
 *
 * @param a Source coordinate system
 * @param b Target coordinate system
 *
 * @return Rotation matrix transforming points from system @p a into @p b
 */
MASM_EXPORT Eigen::Matrix3d rotationMatrix(const CoordinateSystem& a, const CoordinateSystem& b);

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine

#endif
