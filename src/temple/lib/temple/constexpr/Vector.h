// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_VECTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_VECTOR_H

#include "temple/constexpr/Math.h"

#include <array>

/*! @file
 *
 * @brief 3D Vector class with some operations defined
 *
 * Provides a very basic \c constexpr three-dimensional vector class with some
 * limited geometric functionality such as dot- and cross-product. Directly
 * includes the calculation of angles between two vectors.
 */

namespace temple {

//! Constexpr three-dimensional vector math class.
struct Vector {
//!@name State
//!@{
  std::array<double, 3> data;
//!@}

//!@name Constructors
//!@{
  constexpr Vector() : data({{0, 0, 0}}) {}
  constexpr Vector(std::array<double, 3> positions) : data(positions) {}
  constexpr Vector(double x, double y, double z) : data({{x, y, z}}) {}
//!@}

//!@name Common vector operations
//!@{
  //! Computes the dot product with another vector
  PURITY_WEAK constexpr double dot(const Vector& other) const noexcept {
    return (
      this->data[0] * other.data[0]
      + this->data[1] * other.data[1]
      + this->data[2] * other.data[2]
    );
  }

  //! Computes the cross product with another vector
  PURITY_WEAK constexpr Vector cross(const Vector& other) const noexcept {
    return Vector {
      {{
        this -> data[1] * other.data[2] - this -> data[2] * other.data[1],
        this -> data[2] * other.data[0] - this -> data[0] * other.data[2],
        this -> data[0] * other.data[1] - this -> data[1] * other.data[0]
       }}
    };
  }

  //! Calculates the L2 norm of the vector (cartesian length)
  PURITY_WEAK constexpr double norm() const {
    return Math::sqrt(
      this -> data[0] * this -> data[0]
      + this -> data[1] * this -> data[1]
      + this -> data[2] * this -> data[2]
    );
  }

  PURITY_WEAK constexpr Vector operator + (const Vector& other) const noexcept {
    return Vector {
      {{
        this->data[0] + other.data[0],
        this->data[1] + other.data[1],
        this->data[2] + other.data[2]
      }}
    };
  }

  PURITY_WEAK constexpr Vector operator - (const Vector& other) const noexcept {
    return Vector {
      {{
        this->data[0] - other.data[0],
        this->data[1] - other.data[1],
        this->data[2] - other.data[2]
      }}
    };
  }

  PURITY_WEAK constexpr Vector operator * (const double constant) const noexcept {
    return Vector {
      {{
        constant * this->data[0],
        constant * this->data[1],
        constant * this->data[2]
      }}
    };
  }

  //! Division by double operator, throws on division by zero
  PURITY_WEAK constexpr Vector operator / (const double constant) const {
    if(constant == 0) {
      throw "Constexpr::Vector divided by zero!";
    }

    return Vector {
      {{
        this->data[0] / constant,
        this->data[1] / constant,
        this->data[2] / constant
      }}
    };
  }

  //! Unitary inversion operator
  PURITY_WEAK constexpr Vector operator - () const noexcept {
    return Vector {
      {{
        - this -> data[0],
        - this -> data[1],
        - this -> data[2]
      }}
    };
  }
//!@}
};

//! Free-standing binary angle in radians calculation
PURITY_WEAK constexpr double angle(const Vector& a, const Vector& b) {
  return Math::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

} // namespace temple

#endif
