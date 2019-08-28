#include "chemical_symmetries/Recognition.h"

#include "chemical_symmetries/Diophantine.h"
#include "chemical_symmetries/Partitioner.h"

#include "boost/optional.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <vector>
#include "temple/Functional.h"
#include "temple/Optimization/TrustRegion.h"
#include "temple/Optimization/Common.h"

#include "temple/Stringify.h"
#include <iostream>
#include <fstream>
#include <random>

namespace Scine {
namespace Symmetry {

namespace elements {

template<typename EnumType>
constexpr auto underlying(const EnumType e) {
  return static_cast<std::underlying_type_t<EnumType>>(e);
}

inline bool collinear(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
  return std::fabs(std::fabs(a.dot(b) / (a.norm() * b.norm())) - 1) <= 1e-8;
}

inline bool orthogonal(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
  return std::fabs(a.dot(b) / (a.norm() * b.norm())) <= 1e-8;
}

SymmetryElement::Matrix Identity::matrix() const {
  return Matrix::Identity();
}

boost::optional<SymmetryElement::Vector> Identity::vector() const {
  return boost::none;
}

SymmetryElement::Matrix Inversion::matrix() const {
  return -Matrix::Identity();
}

boost::optional<SymmetryElement::Vector> Inversion::vector() const {
  return boost::none;
}

Rotation::Rotation(
  const Eigen::Vector3d& passAxis,
  const unsigned passN,
  const unsigned passPower,
  const bool passReflect
) : axis(passAxis.normalized()),
    n(passN),
    power(passPower),
    reflect(passReflect)
{}

Rotation Rotation::Cn(const Eigen::Vector3d& axis, const unsigned n, const unsigned power) {
  return Rotation(axis, n, power, false);
}

Rotation Rotation::Sn(const Eigen::Vector3d& axis, const unsigned n, const unsigned power) {
  return Rotation(axis, n, power, true);
}

Rotation Rotation::operator * (const Rotation& rhs) const {
  if(collinear(axis, rhs.axis)) {
    if(n == rhs.n) {
      return Rotation(axis, n, power + rhs.power, reflect xor rhs.reflect);
    } else {
      throw std::logic_error("Rotation data model cannot handle collinear multiplication of axes of different order n");
    }
  } else if(orthogonal(axis, rhs.axis)) {
    // Rotate rhs' axis by *this, but keep everything else
    return Rotation(matrix() * rhs.axis, rhs.n, rhs.power, rhs.reflect);
  } else {
    throw std::logic_error("Rotation data model cannot handle non-orthogonal multiplication of rotations");
  }
}

SymmetryElement::Matrix Rotation::matrix() const {
  if(!reflect) {
    return Eigen::AngleAxisd(2 * M_PI * power / n, axis).toRotationMatrix();
  }

  const double angle = 2 * M_PI * power / n;
  const double sine = std::sin(angle);
  const double cosine = std::cos(angle);
  const double onePlusCosine = 1 + cosine;

  const double xx = cosine - axis(0) * axis(0) * onePlusCosine;
  const double yy = cosine - axis(1) * axis(1) * onePlusCosine;
  const double zz = cosine - axis(2) * axis(2) * onePlusCosine;

  const double xy = - axis(0) * axis(1) * onePlusCosine;
  const double xz = - axis(0) * axis(2) * onePlusCosine;
  const double yz = - axis(1) * axis(2) * onePlusCosine;

  const double x = axis(0) * sine;
  const double y = axis(1) * sine;
  const double z = axis(2) * sine;

  Eigen::Matrix3d rotationMatrix;

  rotationMatrix <<
        xx, xy - z, xz + y,
    xy + z,     yy, yz - x,
    xz - y, yz + x,     zz;

  return rotationMatrix;
}

boost::optional<SymmetryElement::Vector> Rotation::vector() const {
  return axis;
}

Reflection::Reflection(const Eigen::Vector3d& passNormal) : normal(passNormal.normalized()) {}

SymmetryElement::Matrix Reflection::matrix() const {
  Eigen::Matrix3d reflection;

  const double normalSquareNorm = normal.squaredNorm();

  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      reflection(i, j) = (i == j ? 1 : 0) - 2 * normal(i) * normal(j) / normalSquareNorm;
    }
  }

  return reflection;
}

boost::optional<SymmetryElement::Vector> Reflection::vector() const {
  if(orthogonal(normal, Eigen::Vector3d::UnitZ())) {
    return normal.cross(Eigen::Vector3d::UnitZ());
  }

  if(orthogonal(normal, Eigen::Vector3d::UnitX())) {
    return normal.cross(Eigen::Vector3d::UnitX());
  }

  if(orthogonal(normal, Eigen::Vector3d::UnitY())) {
    return normal.cross(Eigen::Vector3d::UnitY());
  }

  return boost::none;
}

Rotation operator * (const Rotation& rot, const Reflection& reflection) {
  if(!collinear(rot.axis, reflection.normal)) {
    throw std::logic_error("Cannot handle off-axis Rotation / Reflection combination");
  }

  return Rotation(rot.axis, rot.n, rot.power, !rot.reflect);
}

Rotation operator * (const Reflection& reflection, const Rotation& rot) {
  return rot * reflection;
}

Eigen::Matrix3d improperRotationMatrix(
  const Eigen::Vector3d& axis,
  const double angle
) {
  const double sine = std::sin(angle);
  const double cosine = std::cos(angle);
  const double onePlusCosine = 1 + cosine;

  const double xx = cosine - axis(0) * axis(0) * onePlusCosine;
  const double yy = cosine - axis(1) * axis(1) * onePlusCosine;
  const double zz = cosine - axis(2) * axis(2) * onePlusCosine;

  const double xy = - axis(0) * axis(1) * onePlusCosine;
  const double xz = - axis(0) * axis(2) * onePlusCosine;
  const double yz = - axis(1) * axis(2) * onePlusCosine;

  const double x = axis(0) * sine;
  const double y = axis(1) * sine;
  const double z = axis(2) * sine;

  Eigen::Matrix3d rotationMatrix;

  rotationMatrix <<
        xx, xy - z, xz + y,
    xy + z,     yy, yz - x,
    xz - y, yz + x,     zz;

  return rotationMatrix;
}

Eigen::Matrix3d properRotationMatrix(
  const Eigen::Vector3d& axis,
  const double angle
) {
  return Eigen::AngleAxisd(angle, axis).toRotationMatrix();
}

Eigen::Matrix3d reflectionMatrix(const Eigen::Vector3d& planeNormal) {
  Eigen::Matrix3d reflection;

  const double normalSquareNorm = planeNormal.squaredNorm();

  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      reflection(i, j) = (i == j ? 1 : 0) - 2 * planeNormal(i) * planeNormal(j) / normalSquareNorm;
    }
  }

  return reflection;
}

//! Returns all symmetry elements of a point group
std::vector<std::unique_ptr<SymmetryElement>> symmetryElements(const PointGroup group) {
  using ElementsList = std::vector<std::unique_ptr<SymmetryElement>>;

  auto make = [](auto element) {
    using Type = decltype(element);
    return std::make_unique<Type>(element);
  };

  Identity E {};
  Inversion inversion {};

  const auto e_x = Eigen::Vector3d::UnitX();
  const auto e_y = Eigen::Vector3d::UnitY();
  const auto e_z = Eigen::Vector3d::UnitZ();

  Reflection sigma_xy {e_z};
  Reflection sigma_xz {e_y};
  Reflection sigma_yz {e_x};

  const double tetrahedronAngle = 2 * std::atan(std::sqrt(2));

  auto addProperAxisElements = [&](ElementsList& list, const Eigen::Vector3d& axis, const unsigned n) {
    // C2 gives only a C2, but C3 should also give a C3², etc.
    const Rotation element = Rotation::Cn(axis, n);
    Rotation composite = element;
    for(unsigned i = n; i > 1; --i) {
      list.push_back(make(composite));
      composite = element * composite;
    }
  };

  auto addImproperAxisElements = [&](ElementsList& list, const Eigen::Vector3d& axis, const unsigned n) {
    const Rotation element = Rotation::Sn(axis, n);
    Rotation composite = element;
    for(unsigned i = n; i > 1; --i) {
      list.push_back(make(composite));
      composite = element * composite;
    }
  };

  ElementsList elements;
  elements.push_back(make(E));

  switch(group) {
    case(PointGroup::C1):
      return elements;

    case(PointGroup::Ci):
      {
        elements.push_back(make(inversion));
        return elements;
      }

    case(PointGroup::Cs):
      {
        elements.push_back(make(sigma_xy));
        return elements;
      }

    case(PointGroup::C2):
    case(PointGroup::C3):
    case(PointGroup::C4):
    case(PointGroup::C5):
    case(PointGroup::C6):
    case(PointGroup::C7):
    case(PointGroup::C8):
      {
        const unsigned n = 2 + underlying(group) - underlying(PointGroup::C2);
        addProperAxisElements(elements, e_z, n);
        assert(elements.size() == n);
        return elements;
      }

    case(PointGroup::C2h):
    case(PointGroup::C3h):
    case(PointGroup::C4h):
    case(PointGroup::C5h):
    case(PointGroup::C6h):
    case(PointGroup::C7h):
    case(PointGroup::C8h):
      {
        elements.push_back(make(sigma_xy));
        const unsigned n = 2 + underlying(group) - underlying(PointGroup::C2h);
        std::vector<Rotation> rotations;
        const Rotation element = Rotation::Cn(e_z, n);
        Rotation composite = element;
        for(unsigned i = n; i > 1; --i) {
          rotations.push_back(composite);
          composite = element * composite;
        }
        const unsigned S = rotations.size();
        // Add sigma_xy modified Cn axes
        for(unsigned i = 0; i < S; ++i) {
          rotations.push_back(sigma_xy * rotations.at(i));
        }
        for(auto& rotation : rotations) {
          elements.emplace_back(
            make(std::move(rotation))
          );
        }
        assert(elements.size() == 2 * n);
        return elements;
      }

    case(PointGroup::C2v):
    case(PointGroup::C3v):
    case(PointGroup::C4v):
    case(PointGroup::C5v):
    case(PointGroup::C6v):
    case(PointGroup::C7v):
    case(PointGroup::C8v):
      {
        const unsigned n = 2 + underlying(group) - underlying(PointGroup::C2v);
        addProperAxisElements(elements, e_z, n);
        // Reflection planes include z and increment by pi/n along z
        const auto rotation = Rotation::Cn(e_z, 2 * n);
        Eigen::Vector3d planeNormal = e_y;
        for(unsigned i = 0; i < n; ++i) {
          elements.push_back(
            make(Reflection(planeNormal))
          );
          planeNormal = rotation.matrix() * planeNormal;
        }
        assert(elements.size() == 2 * n);
        return elements;
      }

    case(PointGroup::S4):
    case(PointGroup::S6):
    case(PointGroup::S8):
      {
        const unsigned n = 4 + 2 * (underlying(group) - underlying(PointGroup::S4));
        addImproperAxisElements(elements, e_z, n);
        assert(elements.size() == n);
        return elements;
      }

    case(PointGroup::D2):
    case(PointGroup::D3):
    case(PointGroup::D4):
    case(PointGroup::D5):
    case(PointGroup::D6):
    case(PointGroup::D7):
    case(PointGroup::D8):
      {
        const unsigned n = 2 + underlying(group) - underlying(PointGroup::D2);
        addProperAxisElements(elements, e_z, n);
        // Dn groups have C2 axes along pi/n increments in the xy plane
        const auto rotation = Rotation::Cn(e_z, 2 * n);
        Eigen::Vector3d c2axis = e_x;
        for(unsigned i = 0; i < n; ++i) {
          elements.push_back(make(Rotation::Cn(c2axis, 2)));
          c2axis = rotation.matrix() * c2axis;
        }
        assert(elements.size() == 2 * n);
        return elements;
      }

    case(PointGroup::D2h):
    case(PointGroup::D3h):
    case(PointGroup::D4h):
    case(PointGroup::D5h):
    case(PointGroup::D6h):
    case(PointGroup::D7h):
    case(PointGroup::D8h):
      {
        elements.push_back(make(sigma_xy));
        const unsigned n = 2 + underlying(group) - underlying(PointGroup::D2h);
        elements.reserve(4 * n);
        std::vector<Rotation> rotations;
        const Rotation element = Rotation::Cn(e_z, n);
        Rotation composite = element;
        for(unsigned i = n; i > 1; --i) {
          rotations.push_back(composite);
          composite = element * composite;
        }
        const unsigned S = rotations.size();
        // Generate the S_n axes from the sigma_xy * C_n
        for(unsigned i = 0; i < S; ++i) {
          rotations.push_back(sigma_xy * rotations.at(i));
        }
        for(auto& rotation : rotations) {
          elements.emplace_back(
            make(std::move(rotation))
          );
        }

        /* Dnh groups have C2 axes along pi/n increments in the xy plane
         * and sigma_v planes perpendicular to those C2 axes
         */
        const auto rotation = Rotation::Cn(e_z, 2 * n);
        Eigen::Vector3d c2axis = e_x;
        for(unsigned i = 0; i < n; ++i) {
          elements.push_back(make(Rotation::Cn(c2axis, 2)));
          elements.push_back(
            make(Reflection(
              e_z.cross(c2axis)
            ))
          );
          c2axis = rotation.matrix() * c2axis;
        }

        assert(elements.size() == 4 * n);
        return elements;
      }

    case(PointGroup::D2d):
    case(PointGroup::D3d):
    case(PointGroup::D4d):
    case(PointGroup::D5d):
    case(PointGroup::D6d):
    case(PointGroup::D7d):
    case(PointGroup::D8d):
      {
        const unsigned n = 2 + underlying(group) - underlying(PointGroup::D2d);
        addImproperAxisElements(elements, e_z, 2 * n);
        /* C2 axes */
        const auto rotationMatrix = Rotation::Cn(e_z, 2 * n).matrix();
        Eigen::Vector3d c2axis = e_x;
        for(unsigned i = 0; i < n; ++i) {
          elements.push_back(make(Rotation::Cn(c2axis, 2)));
          c2axis = rotationMatrix * c2axis;
        }
        /* sigma_ds */
        Eigen::Vector3d planeNormal = (e_x + rotationMatrix * e_x).normalized().cross(e_z);
        for(unsigned i = 0; i < n; ++i) {
          elements.push_back(make(Reflection(planeNormal)));
          planeNormal = rotationMatrix * planeNormal;
        }
        assert(elements.size() == 4 * n);
        return elements;
      }

    /* Cubic groups */
    case(PointGroup::T):
      {
        elements.reserve(12);
        const Eigen::Vector3d axis_2 = properRotationMatrix(e_y, tetrahedronAngle) * e_z;
        const Eigen::Vector3d axis_3 = Rotation::Cn(e_z, 3).matrix() * axis_2;
        const Eigen::Vector3d axis_4 = Rotation::Cn(e_z, 3).matrix() * axis_3;
        addProperAxisElements(elements, e_z, 3);
        addProperAxisElements(elements, axis_2, 3);
        addProperAxisElements(elements, axis_3, 3);
        addProperAxisElements(elements, axis_4, 3);
        elements.push_back(make(Rotation::Cn((e_z + axis_2).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_z + axis_3).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_z + axis_4).normalized(), 2)));
        assert(elements.size() == 12);
        return elements;
      }

    case(PointGroup::Td):
      {
        elements.reserve(24);
        Eigen::Matrix<double, 3, Eigen::Dynamic> positions(3, 4);
        positions.col(0) = e_z;
        positions.col(1) = properRotationMatrix(e_y, tetrahedronAngle) * e_z;
        positions.col(2) = Rotation::Cn(e_z, 3).matrix() * positions.col(1);
        positions.col(3) = Rotation::Cn(e_z, 3).matrix() * positions.col(2);
        // C3 axes
        addProperAxisElements(elements, e_z, 3);
        addProperAxisElements(elements, positions.col(1), 3);
        addProperAxisElements(elements, positions.col(2), 3);
        addProperAxisElements(elements, positions.col(3), 3);
        const Eigen::Vector3d axis_12 = (e_z + positions.col(1)).normalized();
        const Eigen::Vector3d axis_13 = (e_z + positions.col(2)).normalized();
        const Eigen::Vector3d axis_14 = (e_z + positions.col(3)).normalized();
        // S4, C2, S4^3
        addImproperAxisElements(elements, axis_12, 4);
        addImproperAxisElements(elements, axis_13, 4);
        addImproperAxisElements(elements, axis_14, 4);
        // Sigma d
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = i + 1; j < 4; ++j) {
            elements.push_back(
              make(Reflection(
                positions.col(i).cross(positions.col(j))
              ))
            );
          }
        }
        assert(elements.size() == 24);
        return elements;
      }
    case(PointGroup::Th):
      {
        // TODO
        assert(elements.size() == 24);
        return elements;
      }

    case(PointGroup::O):
      {
        // TODO
        assert(elements.size() == 24);
        return elements;
      }
    case(PointGroup::Oh):
      {
        elements.push_back(make(inversion));
        elements.reserve(48);
        /* 8 C3 and 8 S6 share the linear combinations of three axes */
        { // +++ <-> ---
          const Eigen::Vector3d axis_ppp = (  e_x + e_y + e_z).normalized();
          addProperAxisElements(elements, axis_ppp, 3);
          elements.push_back(make(Rotation::Sn(axis_ppp, 6)));
          elements.push_back(make(Rotation::Sn(-axis_ppp, 6)));
        }
        { // ++- <-> --+
          const Eigen::Vector3d axis_ppm = (  e_x + e_y - e_z).normalized();
          addProperAxisElements(elements, axis_ppm, 3);
          elements.push_back(make(Rotation::Sn(axis_ppm, 6)));
          elements.push_back(make(Rotation::Sn(-axis_ppm, 6)));
        }
        { // +-+ <-> -+-
          const Eigen::Vector3d axis_pmp = (  e_x - e_y + e_z).normalized();
          addProperAxisElements(elements, axis_pmp, 3);
          elements.push_back(make(Rotation::Sn(axis_pmp, 6)));
          elements.push_back(make(Rotation::Sn(-axis_pmp, 6)));
        }
        { // -++ <-> +--
          const Eigen::Vector3d axis_mpp = (- e_x + e_y + e_z).normalized();
          addProperAxisElements(elements, axis_mpp, 3);
          elements.push_back(make(Rotation::Sn(axis_mpp, 6)));
          elements.push_back(make(Rotation::Sn(-axis_mpp, 6)));
        }
        /* 6 C2 along linear combinations of two axes */
        elements.push_back(make(Rotation::Cn((e_x + e_y).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_x - e_y).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_x + e_z).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_x - e_z).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_y + e_z).normalized(), 2)));
        elements.push_back(make(Rotation::Cn((e_y - e_z).normalized(), 2)));
        /* 6 C4 and 3 C2 (C4^2) along axes */
        addProperAxisElements(elements, e_x, 4);
        addProperAxisElements(elements, e_y, 4);
        addProperAxisElements(elements, e_z, 4);
        /* 6 S4 along axes */
        elements.push_back(make(Rotation::Sn(  e_x, 4)));
        elements.push_back(make(Rotation::Sn(- e_x, 4)));
        elements.push_back(make(Rotation::Sn(  e_y, 4)));
        elements.push_back(make(Rotation::Sn(- e_y, 4)));
        elements.push_back(make(Rotation::Sn(  e_z, 4)));
        elements.push_back(make(Rotation::Sn(- e_z, 4)));
        /* 3 sigma h along combinations of two axes */
        elements.push_back(make(sigma_xy));
        elements.push_back(make(sigma_xz));
        elements.push_back(make(sigma_yz));
        /* 6 sigma d along linear combinations of three axes */
        elements.push_back(make(Reflection((e_x + e_y).cross(e_z))));
        elements.push_back(make(Reflection((e_x - e_y).cross(e_z))));
        elements.push_back(make(Reflection((e_x + e_z).cross(e_y))));
        elements.push_back(make(Reflection((e_x - e_z).cross(e_y))));
        elements.push_back(make(Reflection((e_y + e_z).cross(e_x))));
        elements.push_back(make(Reflection((e_y - e_z).cross(e_x))));
        assert(elements.size() == 48);
        return elements;
      }

    /* Icosahedral groups */
    case(PointGroup::I):
      {
        // TODO
        assert(elements.size() == 60);
        return elements;
      }
    case(PointGroup::Ih):
      {
        // TODO
        assert(elements.size() == 120);
        return elements;
      }

    case(PointGroup::Cinfv):
    case(PointGroup::Dinfh):
      throw std::logic_error("Do not use symmetry elements to analyze linear point groups!");

    default:
      throw std::logic_error("Used invalid PointGroup value");
  }
}

std::map<unsigned, ElementGrouping> npGroupings(
  const std::vector<std::unique_ptr<SymmetryElement>>& elements
) {
  assert(elements.front()->matrix() == elements::Identity().matrix());
  const unsigned E = elements.size();

  /* There can be multiple groupings of symmetry elements of equal l for
   * different points in space. For now, we are ASSUMING that a particular
   * grouping of symmetry elements leads the folded point to lie along the
   * axis defined by the probe point used here to determine the grouping.
   *
   * So we store the point suggested by the element along with its grouping
   * so we can test this theory.
   *
   * TODO mark this resolved after it's been confirmed
   */

  std::map<unsigned, ElementGrouping> npGroupings;

  for(const auto& elementPtr : elements) {
    if(auto axisOption = elementPtr->vector()) {
      PositionCollection mappedPoints(3, E);
      mappedPoints.col(0) = *axisOption;
      unsigned np = 1;
      std::vector<
        std::vector<unsigned>
      > groups {
        {0}
      };

      for(unsigned i = 1; i < E; ++i) {
        Eigen::Vector3d mapped = elements.at(i)->matrix() * mappedPoints.col(0);
        bool found = false;
        for(unsigned j = 0; j < np; ++j) {
          if(mappedPoints.col(j).isApprox(mapped, 1e-8)) {
            found = true;
            groups.at(j).push_back(i);
            break;
          }
        }

        if(!found) {
          mappedPoints.col(np) = mapped;
          ++np;
          groups.push_back(std::vector<unsigned> {i});
        }
      }

      /* If all symmetry operations leave the point in place, this grouping
       * is very uninteresting because we cannot unfold the point at all.
       */
      if(np == 1) {
        continue;
      }

      if(npGroupings.count(np) == 0) {
        ElementGrouping grouping;
        grouping.probePoint = mappedPoints.col(0);
        grouping.groups = std::move(groups);

        npGroupings.emplace(
          np,
          std::move(grouping)
        );
      } else {
        throw std::logic_error("Multiple l mappings found!");
      }
    }
  }

  return npGroupings;
}

} // namespace elements

namespace detail {

Eigen::MatrixXd dropColumn(const Eigen::MatrixXd& matrix, unsigned colToRemove) {
  const unsigned c = matrix.cols();
  const unsigned leftColCount = colToRemove - 1;
  const unsigned rightColCount = c - colToRemove;

  Eigen::MatrixXd copy(matrix.rows(), c - 1);
  copy.leftCols(leftColCount) = matrix.leftCols(leftColCount);
  copy.rightCols(rightColCount) = matrix.rightCols(rightColCount);

  return copy;
}

/*! @brief Normalize positions for continuous symmetry measure analysis
 *
 * Reframes to center of mass frame (although no masses exist) and rescales
 * vectors so that the maximum distance is 1)
 */
PositionCollection normalize(const PositionCollection& positions) {
  const unsigned N = positions.cols();

  // Translate the origin to the average position
  const Eigen::Vector3d center = positions.rowwise().sum() / positions.cols();
  PositionCollection transformed = positions.colwise() - center;

  // Rescale all distances so that the longest is a unit vector
  const double longestDistance = std::sqrt(positions.colwise().squaredNorm().maxCoeff());
  for(unsigned i = 0; i < N; ++i) {
    transformed.col(i) /= longestDistance;
  }

  // Drop any vectors that are very close to the center of mass
  for(int i = 0; i < transformed.cols(); ++i) {
    if(transformed.col(i).norm() < 1e-3) {
      const int r = transformed.rows();
      const int c = transformed.cols() - 1;
      if(i < c) {
        transformed.block(0, i, r, c - i) = transformed.rightCols(c - i);
      }
      transformed.conservativeResize(r, c);
      --i;
    }
  }

  assert(transformed.cols() >= 2);

  return transformed;
}

//! Determine degeneracy of intertial moments
unsigned degeneracy(const Eigen::Vector3d& inertialMoments) {
  constexpr double degeneracyEpsilon = 0.1;
  unsigned mdeg = 0;
  if(
    std::fabs(
      (inertialMoments(2) - inertialMoments(1)) / inertialMoments(2)
    ) <= degeneracyEpsilon
  ) {
    mdeg += 1;
  }
  if(
    std::fabs(
      (inertialMoments(1) - inertialMoments(0)) / inertialMoments(2)
    ) <= degeneracyEpsilon
  ) {
    mdeg += 2;
  }

  return 1 + (mdeg + 1) / 2;
}

//! Determine point group of an approximately linear system
PointGroup linear(const PositionCollection& normalizedPositions) {
  const unsigned P = normalizedPositions.cols();
  for(unsigned i = 0; i < P / 2; ++i) {
    const unsigned opposite = P - i - 1;
    if(!normalizedPositions.col(i).isApprox(-normalizedPositions.col(opposite), 1e-4)) {
      /* Positions do not have inversion symmetry */
      return PointGroup::Cinfv;
    }
  }

  return PointGroup::Dinfh;
}

} // namespace detail

namespace csm {

double all_symmetry_elements(
  const PositionCollection& normalizedPositions,
  const PointGroup pointGroup,
  std::vector<unsigned> particleIndices
) {
}

double grouped_symmetry_elements(
  const PositionCollection& normalizedPositions,
  const PointGroup pointGroup,
  std::vector<unsigned> particleIndices,
  const std::vector<std::unique_ptr<elements::SymmetryElement>>& elements,
  const elements::ElementGrouping& elementGrouping
) {
  double value = 100;
}

double point_group(
  const PositionCollection& normalizedPositions,
  const PointGroup pointGroup
) {
  const auto elements = elements::symmetryElements(pointGroup);
  const auto npGroups = npGroupings(elements);

  const unsigned G = elements.size();
  const unsigned P = normalizedPositions.cols();

  /* Rough algorithm:
   * - Select an orientation of the point group
   * - For each possible subdivision of points:
   *   - Minimize CSM for each set according to either G = P or G = l * P
   *     over all permutations
   * - Simplex minimize over point group orientation
   */

  // Eigen::Vector3d pointGroupOrientation = Eigen::Vector3d::UnitZ();

  /* Subdivide P */
  // The sizes of groups that we can subdivide into must contain the whole set
  std::vector<unsigned> subdivisionGroupSizes {G};
  // Add npGroupings' group sizes to subdivisionGroupSizes
  for(const auto& groupMapPair : npGroups) {
    subdivisionGroupSizes.push_back(groupMapPair.first);
  }

  /* If there is one more point than ony of the possible groupings, then one
   * point can reasonably be symmetrized to the origin.
   */
  if(
    temple::any_of(
      subdivisionGroupSizes,
      [P](const unsigned size) { return size == P + 1; }
    )
  ) {
    subdivisionGroupSizes.push_back(1);
  }

  std::sort(
    std::begin(subdivisionGroupSizes),
    std::end(subdivisionGroupSizes),
    std::greater<>()
  );

  /* First check if we can solve the diophantine at all, for that the greatest
   * common divisor has to divide P
   */
  if(!diophantine::has_solution(subdivisionGroupSizes, P)) {
    std::cout << "Diophantine specified by " << temple::stringify(subdivisionGroupSizes) << " = " << P << " has no solutions\n";
    return 100.0;
  }

  /* A number of different cases can arise here, and we have to consider them
   * VERY CAREFULLY before proceeding.
   *
   * - If there are two or more elements in subdivisionGroupSizes, we can
   *   enumerate all solutions of the diophantine
   *   BUT skip those where the coefficient of group size 1 (if present) is greater than 1
   *   (i.e. it is NOT reasonable to symmetrize more than one point to the
   *   origin)
   * - If there is only one element in subdivisionGroupSizes, then P is G and
   *   we do not have to subdivide at all -> directly match points to symmetry
   *   elements
   */

  if(subdivisionGroupSizes.size() == 1) {
    // call p = g minimization
  } else {
    std::vector<unsigned> subdivisionMultipliers;
    std::cout << temple::stringify(subdivisionGroupSizes) << " subdivisions\n";

    if(diophantine::first_solution(subdivisionMultipliers, subdivisionGroupSizes, P)) {
      do {
        std::cout << temple::stringify(subdivisionMultipliers) << "\n";

        /* We have a composition of groups to subdivide our points:
         *
         * E.g.: P = 12
         * subdivisionGroupSizes {4, 3, 2}
         * one solution: subdivisionMultipliers {1, 2, 1}
         *
         * Next we make every possible partition of our points into those
         * respective group sizes. We populate a flat map of point index to
         * group size index according to the group sizes and multipliers:
         *
         * p0 = 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2
         *
         * Of this, there are 12! / (4! 6! 2!) = 13860 permutations.
         */
        std::vector<unsigned> flatGroupMap;
        flatGroupMap.reserve(P);
        for(unsigned i = 0; i < subdivisionMultipliers.size(); ++i) {
          const unsigned groupMultiplier = subdivisionMultipliers.at(i);

          if(groupMultiplier == 0) {
            continue;
          }

          const unsigned groupSize = subdivisionGroupSizes.at(i);
          flatGroupMap.resize(flatGroupMap.size() + groupSize * groupMultiplier, i);
        }
        assert(flatGroupMap.size() == P);

        do {
          /* For each of these permutations, we sub-partition each group using
           * Partitioner. This is necessary to treat multipliers > 1 correctly.
           *
           * In the previous case, where we have a two groups of size three, we
           * have 10 sub-partitions to consider. All other groups with
           * multiplier == 1 have only a single sub-partition.
           */

          /* Collect the indices mapped to groups of equal size */
          const unsigned numSizeGroups = subdivisionGroupSizes.size();
          std::vector<
            std::vector<unsigned>
          > sameSizeIndexGroups(numSizeGroups);
          for(unsigned i = 0; i < P; ++i) {
            sameSizeIndexGroups.at(flatGroupMap.at(i)).push_back(i);
          }

          for(unsigned i = 0; i < numSizeGroups; ++i) {
            const unsigned groupSize = subdivisionGroupSizes.at(i);
            const unsigned multiplier = subdivisionMultipliers.at(i);
            const auto& sameSizeParticleIndices = sameSizeIndexGroups.at(i);

            Partitioner partitioner {multiplier, groupSize};

            do {
              // Group the particle indices
              auto particleGroups = temple::map(
                partitioner.partitions(),
                [&](const auto& group) -> std::vector<unsigned> {
                  return temple::map(
                    group,
                    [&](const unsigned indexOfParticle) -> unsigned {
                      return sameSizeParticleIndices.at(indexOfParticle);
                    }
                  );
                }
              );


            } while(partitioner.next_partition());
          }
        } while(std::next_permutation(std::begin(flatGroupMap), std::end(flatGroupMap)));
      } while(diophantine::next_solution(subdivisionMultipliers, subdivisionGroupSizes, P));
    }
  }

  return 100.0;
}

namespace cn {

PositionCollection symmetrize(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis,
  const unsigned n,
  const std::vector<
    std::vector<unsigned>
  > groups
) {
  PositionCollection symmetrizedPositions = normalizedPositions;

  for(const auto& group : groups) {
    const double angleRadians = 2 * M_PI / n;
    /* Fold */
    // The first position is unchanged, hence we start with i = 1
    for(unsigned i = 1; i < n; ++i) {
      symmetrizedPositions.col(group.front()) += Eigen::AngleAxisd(i * angleRadians, axis) * normalizedPositions.col(group.at(i));
    }

    // Average the folded positions
    symmetrizedPositions.col(group.front()) /= n;

    /* Unfold */
    for(unsigned i = 1; i < n; ++i) {
      symmetrizedPositions.col(group.at(i)) = Eigen::AngleAxisd(-i * angleRadians, axis) * symmetrizedPositions.col(group.front());
    }
  }

  return symmetrizedPositions;
}

double evaluatePermutation(
  const PositionCollection& normalizedPositions,
  const Eigen::Vector3d& axis,
  const unsigned n,
  const std::vector<unsigned>& groupPermutation
) {
  assert(n >= 2);

  const double angleRadians = 2 * M_PI / n;
  Eigen::Vector3d averagePoint = normalizedPositions.col(groupPermutation.front());

  /* Fold and average */
  averagePoint = normalizedPositions.col(groupPermutation.front());
  for(unsigned i = 1; i < n; ++i) {
    averagePoint.noalias() += Eigen::AngleAxisd(i * angleRadians, axis) * normalizedPositions.col(groupPermutation.at(i));
  }
  averagePoint /= n;

  /* Calculate CSM while unfolding */
  double csm = (
    normalizedPositions.col(groupPermutation.front())
    - averagePoint
  ).squaredNorm();

  for(unsigned i = 1; i < n; ++i) {
    csm += (
      normalizedPositions.col(groupPermutation.at(i))
      - Eigen::AngleAxisd(i * angleRadians, -axis) * averagePoint
    ).squaredNorm();
  }

  return csm;
}

/**
 * @brief Functor for minimizing csm as a function of cn axis
 */
struct AxisMinimizationFunctor {
  using GroupsType = std::vector<
    std::vector<unsigned>
  >;
  std::reference_wrapper<const PositionCollection> normalizedPositionsRef;
  std::reference_wrapper<const GroupsType> groupsRef;
  unsigned axis_n;

  AxisMinimizationFunctor(
    const PositionCollection& positions,
    const GroupsType& groups,
    const unsigned axisOrder
  ) : normalizedPositionsRef(positions),
      groupsRef(groups),
      axis_n(axisOrder)
  {
    assert(
      temple::all_of(groups,
        [axisOrder](const auto& g) { return g.size() == axisOrder; }
      )
    );
  }

  double evaluate(const Eigen::VectorXd& parameters) {
    const Eigen::Vector3d axis {
      std::sin(parameters(0)) * std::cos(parameters(1)),
      std::sin(parameters(0)) * std::sin(parameters(1)),
      std::cos(parameters(0))
    };

    double csm = 0;
    for(const auto& group : groupsRef.get()) {
      csm += evaluatePermutation(
        normalizedPositionsRef.get(),
        axis,
        axis_n,
        group
      );
    }

    return csm;
  }

  void numericalEvaluation(
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    value = evaluate(parameters);
    gradient = temple::optimization::numericalGradient(
      [this](const Eigen::VectorXd& p) -> double {
        return evaluate(p);
      },
      parameters
    );
    hessian = temple::optimization::numericalHessian(
      [this](const Eigen::VectorXd& p) -> double {
        return evaluate(p);
      },
      parameters
    );
  }

  void explicitEvaluation(
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    assert(parameters.size() == 2);
    assert(gradient.size() == 2);
    assert(hessian.rows() == 2 && hessian.cols() == 2);

    const double& theta = parameters(0);
    const double& phi = parameters(1);
    const Eigen::Vector3d n {
      std::sin(theta) * std::cos(phi),
      std::sin(theta) * std::sin(phi),
      std::cos(theta)
    };

    assert(std::fabs(n.squaredNorm() - 1) <= 1e-10);

    value = 0;

    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    Eigen::Vector3d h = Eigen::Vector3d::Zero();
    const double alpha = 2.0 * M_PI / axis_n;

    const auto& groups = groupsRef.get();
    const auto& normalizedPositions = normalizedPositionsRef.get();

    for(const auto& group : groups) {
      Eigen::Vector3d a, b, c;
      for(unsigned i = 0; i < axis_n; ++i) {
        a.setZero();
        b.setZero();
        c.setZero();

        /* Calculate a, b and c */
        // Terms of j < i
        for(unsigned j = 0; j < i; ++j) {
          const double compoundAngle = (axis_n + j - i) * alpha;
          const double cosine = std::cos(compoundAngle);
          const auto& position = normalizedPositions.col(group.at(j));

          a.noalias() += cosine * position;
          b.noalias() += (1 - cosine) * position;
          c.noalias() += std::sin(compoundAngle) * position;
        }
        // Terms for j == i
        a.noalias() += (1 - static_cast<int>(axis_n)) * normalizedPositions.col(group.at(i));
        // Terms for j > i
        for(unsigned j = i + 1; j < axis_n; ++j) {
          const double compoundAngle = (j - i) * alpha;
          const double cosine = std::cos(compoundAngle);
          const auto& position = normalizedPositions.col(group.at(j));

          a.noalias() += cosine * position;
          b.noalias() += (1 - cosine) * position;
          c.noalias() += std::sin(compoundAngle) * position;
        }

        // Average out all vectors
        a /= axis_n;
        b /= axis_n;
        c /= axis_n;

        // Add value contributions
        value += (a + n.dot(b) * n + n.cross(c)).squaredNorm();

        /* Add contributions to M and h */
        const Eigen::Matrix3d intermediate = a * b.transpose();
        M.noalias() += (
          a * a.transpose()
          - c * c.transpose()
          + intermediate
          + intermediate.transpose()
        );
        h.noalias() += c.cross(a);
      }
    }

    // Calculate gradient components
    const Eigen::Vector3d Mnph = M * n + h;
    const Eigen::Vector3d nPartialTheta {
      std::cos(theta) * std::cos(phi),
      std::cos(theta) * std::sin(phi),
      - std::sin(theta)
    };
    const Eigen::Vector3d nPartialPhi {
      - std::sin(theta) * std::sin(phi),
      std::sin(theta) * std::cos(phi),
      0.0
    };

    gradient(0) = 2 * Mnph.dot(nPartialTheta);
    gradient(1) = 2 * Mnph.dot(nPartialPhi);

    // Calculate hessian components
    const Eigen::Vector3d nPartialThetaPartialPhi {
      - std::cos(theta) * std::sin(phi),
      std::cos(theta) * std::cos(phi),
      0
    };

    // Resolve a 1x1 matrix to its single entry
    auto resolve = [](const Eigen::MatrixXd& matr) -> double {
      assert(matr.rows() == 1 && matr.cols() == 1);
      return matr(0);
    };

    hessian(0, 0) = 2 * (
      resolve(nPartialTheta.transpose() * M * nPartialTheta)
      + Mnph.dot(-n)
    );
    hessian(1, 1) = 2 * (
      resolve(nPartialPhi.transpose() * M * nPartialPhi)
      + Mnph.dot(Eigen::Vector3d {-n(0), -n(1), 0})
    );
    hessian(0, 1) = 2 * (
      resolve(nPartialTheta.transpose() * M * nPartialPhi)
      + Mnph.dot(nPartialThetaPartialPhi)
    );
    hessian(1, 0) = 2 * (
      resolve(nPartialPhi.transpose() * M * nPartialTheta)
      + Mnph.dot(nPartialThetaPartialPhi)
    );
  }

  void operator() (
    const Eigen::VectorXd& parameters,
    double& value,
    Eigen::Ref<Eigen::VectorXd> gradient,
    Eigen::Ref<Eigen::MatrixXd> hessian
  ) {
    numericalEvaluation(parameters, value, gradient, hessian);
    //explicitEvaluation(parameters, value, gradient, hessian);
  }
};

/**
 * @brief Checker for LBFGS minimization
 */
struct AxisMinimizationChecker {
  using FloatType = double;
  using VectorType = Eigen::VectorXd;

  bool shouldContinue(
    const unsigned iteration,
    const FloatType /* value */,
    const VectorType& gradient
  ) {
    return (
      iteration <= 1000
      && gradient.squaredNorm() > 1e-3
    );
  }
};

//! Data struct for axis minimization
struct axis_optimization_t {
  Eigen::Vector3d axis;
  double csm;
};

axis_optimization_t optimize_axis_two_parameters(
  const PositionCollection& normalizedPositions,
  const unsigned n,
  const std::vector<
    std::vector<unsigned>
  >& groups,
  const Eigen::Vector3d& initialAxis
) {
  assert(std::fabs(initialAxis.squaredNorm() - 1) <= 1e-5);

  Eigen::VectorXd parameters (2);
  // Theta from z = cos(theta)
  parameters(0) = std::acos(initialAxis(2));
  // Phi from atan(y / x)
  parameters(1) = std::atan2(initialAxis(1), initialAxis(0));

  AxisMinimizationFunctor functor {
    normalizedPositions,
    groups,
    n
  };

  auto optimizationResult = temple::TrustRegionOptimizer<>::minimize(
    parameters,
    functor,
    AxisMinimizationChecker {}
  );

  if(optimizationResult.iterations >= 1000) {
    throw std::logic_error("Could not minimize axis! Maximum iterations reached.");
  }

  const double theta = parameters(0);
  const double phi = parameters(1);
  Eigen::Vector3d axis {
    std::sin(theta) * std::cos(phi),
    std::sin(theta) * std::sin(phi),
    std::cos(theta)
  };

  return {
    axis,
    100 * optimizationResult.value / (groups.size() * n)
  };
}

struct cn_csm_t {
  double csm;
  std::vector<
    std::vector<unsigned>
  > groups;
};

//! @brief Greedy csm minimization at fixed axis. Not fully permutational.
cn_csm_t fixedAxisGreedy(
  const PositionCollection& normalizedPositions,
  const unsigned n,
  const Eigen::Vector3d& axis
) {
  std::mt19937_64 urbg;
  std::random_device rd;
  urbg.seed(rd());

  assert(n >= 2);
  const unsigned points = normalizedPositions.cols();
  assert(n <= points);

  // Axis vector must be normalized!
  assert(std::fabs(axis.norm() - 1) < 1e-10);

  /* On-axis points should be ignored */
  Eigen::ParametrizedLine<double, 3> axisLine(
    Eigen::Vector3d::Zero(),
    axis
  );

  std::vector<unsigned> offAxisPointIndices;
  offAxisPointIndices.reserve(points);
  for(unsigned i = 0; i < points; ++i) {
    if(axisLine.distance(normalizedPositions.col(i)) > 0.1) {
      offAxisPointIndices.push_back(i);
    }
  }

  const unsigned validPoints = offAxisPointIndices.size();
  /* If there are fewer off-axis points than the rotation order being tested,
   * the element is clearly not present
   */
  if(validPoints < n) {
    cn_csm_t result;
    result.csm = 100;
    return result;
  }

  /* Precalculate the required rotation matrices */
  const double angleRadians = 2 * M_PI / n;
  std::vector<Eigen::Matrix3d> foldMatrices(n - 1);
  for(unsigned i = 0; i < n - 1; ++i) {
    foldMatrices.at(i) = Eigen::AngleAxisd((i + 1) * angleRadians, axis).toRotationMatrix();
  }

  std::vector<Eigen::Matrix3d> unfoldMatrices(n - 1);
  for(unsigned i = 0; i < n - 1; ++i) {
    unfoldMatrices.at(i) = Eigen::AngleAxisd((i + 1) * angleRadians, -axis).toRotationMatrix();
  }

  /* Divide off-axis particle indices into groups of size n each */
  const unsigned groups = validPoints / n;
  assert(groups != 0);

  cn_csm_t result;
  result.csm = 100;
  result.groups.resize(groups);

  std::vector<unsigned> groupIndices(validPoints);
  for(unsigned i = 0; i < validPoints; ++i) {
    groupIndices[i] = i / n;
  }

  do {
    /* Groups are subdivided, but unordered. Calculate csm for each subgroup
     * and minimize csm over all permutations
     */
    double groupCollectiveCSM = 0;
    for(unsigned g = 0; g < groups; ++g) {
      // Collect all indices of the current group number
      std::vector<unsigned> group;
      group.reserve(n);
      for(unsigned i = 0; i < validPoints; ++i) {
        if(groupIndices.at(i) == g) {
          group.push_back(offAxisPointIndices.at(i));
        }
      }

      auto calculateCSM = [&](const std::vector<unsigned>& permutation) -> double {
        Eigen::Vector3d averagePoint = normalizedPositions.col(permutation.front());

        /* Fold and average in-place */
        averagePoint.noalias() = normalizedPositions.col(permutation.front());
        for(unsigned i = 1; i < n; ++i) {
          averagePoint.noalias() += foldMatrices.at(i - 1) * normalizedPositions.col(permutation.at(i));
        }
        averagePoint /= n;

        /* Calculate CSM while unfolding */
        double csm = (
          normalizedPositions.col(permutation.front())
          - averagePoint
        ).squaredNorm();

        for(unsigned i = 1; i < n; ++i) {
          csm += (
            normalizedPositions.col(permutation.at(i))
            - unfoldMatrices.at(i - 1) * averagePoint
          ).squaredNorm();
        }

        return csm;
      };

      assert(group.size() == n);

      std::shuffle(std::begin(group), std::end(group), urbg);

      bool foundBetterAdjacentPermutation = false;
      double currentCSM = calculateCSM(group);
      do {
        foundBetterAdjacentPermutation = false;
        for(unsigned i = 0; i < n && !foundBetterAdjacentPermutation; ++i) {
          for(unsigned j = i + 1; j < n; ++j) {
            std::swap(group.at(i), group.at(j));

            double adjacentCSM = calculateCSM(group);
            if(adjacentCSM < currentCSM) {
              currentCSM = adjacentCSM;
              foundBetterAdjacentPermutation = true;
              break;
            }

            std::swap(group.at(i), group.at(j));
          }
        }
      } while(foundBetterAdjacentPermutation);
      result.groups.at(g) = std::move(group);
      groupCollectiveCSM += currentCSM;
    }

    result.csm = std::min(result.csm, groupCollectiveCSM);
  } while(std::next_permutation(std::begin(groupIndices), std::end(groupIndices)));

  result.csm *= 100.0 / (groups * n);
  return result;
}

//! Minimizes CSM at fixed axis. Fully permutational
double fixedAxis(
  const PositionCollection& normalizedPositions,
  const unsigned n,
  const Eigen::Vector3d& axis
) {
  assert(n >= 2);
  const unsigned points = normalizedPositions.cols();
  assert(n <= points);

  // Axis vector must be normalized!
  assert(std::fabs(axis.norm() - 1) < 1e-10);

  /* On-axis points should be ignored */
  Eigen::ParametrizedLine<double, 3> axisLine(
    Eigen::Vector3d::Zero(),
    axis
  );

  std::vector<unsigned> offAxisPointIndices;
  offAxisPointIndices.reserve(points);
  for(unsigned i = 0; i < points; ++i) {
    if(axisLine.distance(normalizedPositions.col(i)) > 0.1) {
      offAxisPointIndices.push_back(i);
    }
  }

  const unsigned validPoints = offAxisPointIndices.size();
  /* If there are fewer off-axis points than the rotation order being tested,
   * the element is clearly not present
   */
  if(validPoints < n) {
    return 100;
  }

  std::vector<unsigned> groupIndices(validPoints);
  for(unsigned i = 0; i < validPoints; ++i) {
    groupIndices[i] = i / n;
  }

  /* Precalculate the required rotation matrices */
  const double angleRadians = 2 * M_PI / n;
  std::vector<Eigen::Matrix3d> foldMatrices(n - 1);
  for(unsigned i = 0; i < n - 1; ++i) {
    foldMatrices.at(i) = Eigen::AngleAxisd((i + 1) * angleRadians, axis).toRotationMatrix();
  }

  std::vector<Eigen::Matrix3d> unfoldMatrices(n - 1);
  for(unsigned i = 0; i < n - 1; ++i) {
    unfoldMatrices.at(i) = Eigen::AngleAxisd((i + 1) * angleRadians, -axis).toRotationMatrix();
  }

  /* Divide off-axis particle indices into groups of size n each */
  const unsigned groups = validPoints / n;
  double overallLowestCSM = 100;
  do {
    /* Groups are subdivided, but unordered. Calculate csm for each subgroup
     * and minimize csm there over all permutations
     */
    double groupCollectiveCSM = 0;
    for(unsigned g = 0; g < groups; ++g) {
      // Collect all indices of the current group number
      std::vector<unsigned> group;
      group.reserve(n);
      for(unsigned i = 0; i < validPoints; ++i) {
        if(groupIndices.at(i) == g) {
          group.push_back(offAxisPointIndices.at(i));
        }
      }

      assert(group.size() == n);

      /* All point indices of the current group are ordered ascending right now,
       * so we can iterate through their permutations:
       */
      double lowestCSM = 100;
      do {
        /* Fold */
        Eigen::MatrixXd foldedPositions(3, n);
        foldedPositions.col(0) = normalizedPositions.col(group.front());

        // The first position is unchanged, hence we start with i = 1
        for(unsigned i = 1; i < n; ++i) {
          foldedPositions.col(i) = foldMatrices.at(i - 1) * normalizedPositions.col(group.at(i));
        }

        /* Average */
        Eigen::MatrixXd unfoldedPositions(3, n);
        unfoldedPositions.col(0) = foldedPositions.rowwise().sum() / n;

        /* Unfold */
        for(unsigned i = 1; i < n; ++i) {
          unfoldedPositions.col(i) = unfoldMatrices.at(i - 1) * unfoldedPositions.col(0);
        }

        /* Calculate CSM */
        double csm = 0;
        for(unsigned i = 0; i < n; ++i) {
          csm += (normalizedPositions.col(group.at(i)) - unfoldedPositions.col(i)).squaredNorm();
        }
        lowestCSM = std::min(lowestCSM, csm);
      } while(std::next_permutation(std::begin(group), std::end(group)));

      groupCollectiveCSM += lowestCSM;
    }

    overallLowestCSM = std::min(overallLowestCSM, groupCollectiveCSM);
  } while(std::next_permutation(std::begin(groupIndices), std::end(groupIndices)));

  return 100 * overallLowestCSM / (groups * n);
}

} // namespace cn

} // namespace csm

PointGroup analyze(const PositionCollection& positions) {
  auto transformed = detail::normalize(positions);

  const unsigned N = transformed.cols();
  assert(N > 1);

  // Construct the lower triangle of the inertial matrix in the new frame
  Eigen::Matrix3d inertialMatrix = Eigen::Matrix3d::Zero(3, 3);
  for(unsigned i = 0; i < N; ++i) {
    const auto& vec = transformed.col(i);
    inertialMatrix(0, 0) += vec(1) * vec(1) + vec(2) * vec(2); // xx
    inertialMatrix(1, 1) += vec(0) * vec(0) + vec(2) * vec(2); // yy
    inertialMatrix(2, 2) += vec(0) * vec(0) + vec(1) * vec(1); // zz
    inertialMatrix(1, 0) += vec(0) * vec(1); // xy
    inertialMatrix(2, 0) += vec(0) * vec(2); // xz
    inertialMatrix(2, 1) += vec(1) * vec(2); // yz
  }

  // Decompose the inertial matrix to get principal axes and inertial moments
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> decomposition(inertialMatrix);

  for(unsigned i = 0; i < 3; ++i) {
    std::cout << "Eigenvalue " << i
      << " has value " << decomposition.eigenvalues()(i)
      << " and principal axis " << decomposition.eigenvectors().col(i).transpose() << "\n";
  }

  const unsigned degeneracy = detail::degeneracy(decomposition.eigenvalues());

  std::cout << "Degeneracy degree: " << degeneracy << "\n";

  /* If degeneracy == 2, rotate unique axis to coincide with z, and one of the
   * degenerate axes to coincide with x. There could be a C2 on x (?)
   *
   * If degeneracy == 3, look for C3 axes parallel to vectors, then rotate
   * the coordinates so that they coincide with x, ±x, ±x. (???)
   */

  /* Key insights:
   * An element of symmetry must be parallel or perpendicular to the
   * combination of an adequate number of vectors.
   *
   * Examples:
   * - If you have C3h, and you're testing for a C2 axis perpendicular to the
   *   C3 axis, it will be parallel to a vector i + j (i != j) if n = 6
   *   or to a single i if n = 3.
   * - If you are looking for a mirror plane in C3v, it will be perpendicular
   *   to a vector i - j (i != j)
   * - A C3 axis in a cubic group will be parallel to i if n = 4 or to some
   *   i + j + k (i != j != k) if n > 4
   */

  /* Linearity can be checked really easily by looking at the principal moments
   * of inertia. They're sorted in ascending order by the solver.
   *
   * - Linear molecules: IA << IB == IC (typically IA ≃ 0)
   * - Spherical tops: IA = IB = IC => only Td, Oh
   * - Symmetric tops: IA = IB or IB = IC. Must have a threefold or higher
   *   rotation axis by definition.
   *   - Oblate (disc): IA = IB < IC
   *   - Prolate (rugby football): IA < IB = IC
   * - Asymmetric tops
   */


  /* Linear cases are best tested separately since we cannot represent infinite
   * symmetry elements
   */
  if(decomposition.eigenvalues()(0) < 0.1 && degeneracy == 2) {
    /* Do an extra test for collinearity */
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> rankDecomposition(transformed);
    rankDecomposition.setThreshold(1e-3);
    if(rankDecomposition.rank() <= 1) {
      return detail::linear(transformed);
    }
  }

  /* If the positions are already symmetrical enough to be recognized as a
   * spherical top, then we might as well make use of the knowledge
   */
  //if(degeneracy == 3) {
  //  const double tetrahedral_csm = csm(transformed, PointGroup::Td);
  //  const double octahedral_csm = csm(transformed, PointGroup::Oh);
  //  // TODO icosahedral is also a spherical top!
  //  return tetrahedral_csm < octahedral_csm ? PointGroup::Td : PointGroup::Oh;
  //}

  /* Find a main axis and rotate it to z */
  /* Rotate secondary axis to x */
  for(unsigned axisIndex = 0; axisIndex < 3; ++axisIndex) {
    std::cout << "At axis " << axisIndex << " along " << decomposition.eigenvectors().col(axisIndex).transpose() << "\n";
    for(unsigned n = 2; n <= N; ++n) {
      const double axis_csm = csm::cn::fixedAxis(
        transformed,
        n,
        decomposition.eigenvectors().col(axisIndex)
      );

      const auto csm_greedy = csm::cn::fixedAxisGreedy(
        transformed,
        n,
        decomposition.eigenvectors().col(axisIndex)
      );

      if(std::fabs(axis_csm - csm_greedy.csm) > 0.1) {
        std::cout << "C" << n << " reg and greedy differ strongly: " << axis_csm  << " and " << csm_greedy.csm << "\n";
      }

      std::cout << "S(C" << n << ") = " << csm_greedy.csm << "\n";

      if(axis_csm < 0.1) {
        // Cn symmetry csm
        const PointGroup cn_point_group = static_cast<PointGroup>(
          elements::underlying(PointGroup::C2) + n - 2
        );
        double symm_csm = csm::point_group(
          transformed,
          cn_point_group
        );
        std::cout << "C" << n << " point group symmetry csm: " << symm_csm << "\n";
      }
    }
  }

  return PointGroup::C1;
}

} // namespace Symmetry
} // namespace Scine
