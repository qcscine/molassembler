/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief For transformations between coordinate systems
 */

#include <Eigen/Geometry>
#include <Eigen/Core>

namespace Scine {
namespace shapes {

/**
 * @brief Coordinate system axis data class
 */
struct CoordinateSystem {
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
Eigen::Matrix3d rotationMatrix(const CoordinateSystem& a, const CoordinateSystem& b);

} // namespace shapes
} // namespace Scine
