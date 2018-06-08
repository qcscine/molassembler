#ifndef INCLUDE_CONSTEXPR_MAGIC_VECTOR_H
#define INCLUDE_CONSTEXPR_MAGIC_VECTOR_H

#include "Math.h"

#include <array>

/*! @file
 *
 * Provides a very basic \c constexpr three-dimensional vector class with some
 * limited geometric functionality such as dot- and cross-product. Directly
 * includes the calculation of angles between two vectors.
 */

namespace temple {

//! Constexpr three-dimensional vector math class.
struct Vector {
  std::array<double, 3> data;

  constexpr Vector() : data({{0, 0, 0}}) {}
  constexpr Vector(std::array<double, 3> positions) : data(positions) {}
  constexpr Vector(double x, double y, double z) : data({{x, y, z}}) {}

  //! Computes the dot product with another vector
  constexpr double dot(const Vector& other) const noexcept PURITY_WEAK {
    return (
      this->data[0] * other.data[0]
      + this->data[1] * other.data[1]
      + this->data[2] * other.data[2]
    );
  }

  //! Computes the cross product with another vector
  constexpr Vector cross(const Vector& other) const noexcept PURITY_WEAK {
    return Vector {
      {{
        this -> data[1] * other.data[2] - this -> data[2] * other.data[1],
        this -> data[2] * other.data[0] - this -> data[0] * other.data[2],
        this -> data[0] * other.data[1] - this -> data[1] * other.data[0]
       }}
    };
  }

  //! Calculates the L2 norm of the vector (cartesian length)
  constexpr double norm() const PURITY_WEAK {
    return Math::sqrt(
      this -> data[0] * this -> data[0]
      + this -> data[1] * this -> data[1]
      + this -> data[2] * this -> data[2]
    );
  }

  constexpr Vector operator + (const Vector& other) const noexcept PURITY_WEAK {
    return Vector {
      {{
        this->data[0] + other.data[0],
        this->data[1] + other.data[1],
        this->data[2] + other.data[2]
      }}
    };
  }

  constexpr Vector operator - (const Vector& other) const noexcept PURITY_WEAK {
    return Vector {
      {{
        this->data[0] - other.data[0],
        this->data[1] - other.data[1],
        this->data[2] - other.data[2]
      }}
    };
  }

  constexpr Vector operator * (const double& constant) const noexcept PURITY_WEAK {
    return Vector {
      {{
        constant * this->data[0],
        constant * this->data[1],
        constant * this->data[2]
      }}
    };
  }

  //! Division by double operator, throws on division by zero
  constexpr Vector operator / (const double& constant) const PURITY_WEAK {
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
  constexpr Vector operator - () const noexcept PURITY_WEAK {
    return Vector {
      {{
        - this -> data[0],
        - this -> data[1],
        - this -> data[2]
      }}
    };
  }
};

//! Free-standing binary angle in radians calculation
constexpr double angle(const Vector& a, const Vector& b) PURITY_WEAK;
constexpr double angle(const Vector& a, const Vector& b) {
  return Math::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

} // namespace temple

#endif
