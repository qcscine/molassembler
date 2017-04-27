#ifndef INCLUDE_CONSTEXPR_MAGIC_VECTOR_H
#define INCLUDE_CONSTEXPR_MAGIC_VECTOR_H

/* TODO
 * - remove dependency
 */
#include "static_math/cmath.h"

#include "Math.h"

#include <array>

namespace ConstexprMagic {

struct Vector {
  std::array<double, 3> data;

  constexpr Vector(std::array<double, 3> positions) : data(positions) {}

  constexpr double dot(const Vector& other) const {
    return (
      this->data[0] * other.data[0]
      + this->data[1] * other.data[1]
      + this->data[2] * other.data[2]
    );
  }

  constexpr Vector cross(const Vector& other) const {
    return Vector {
      {
        this -> data[1] * other.data[2] - this -> data[2] * other.data[1],
        this -> data[2] * other.data[0] - this -> data[0] * other.data[2],
        this -> data[0] * other.data[1] - this -> data[1] * other.data[0]
      }
    };
  }

  constexpr double norm() const {
    return smath::sqrt(
      this -> data[0] * this -> data[0]
      + this -> data[1] * this -> data[1]
      + this -> data[2] * this -> data[2]
    );
  }

  constexpr Vector operator + (const Vector& other) const {
    return Vector {
      {
        this->data[0] + other.data[0],
        this->data[1] + other.data[1],
        this->data[2] + other.data[2]
      }
    };
  }

  constexpr Vector operator - (const Vector& other) const {
    return Vector {
      {
        this->data[0] - other.data[0],
        this->data[1] - other.data[1],
        this->data[2] - other.data[2]
      }
    };
  }

  constexpr Vector operator * (const double& constant) const {
    return Vector {
      {
        constant * this->data[0],
        constant * this->data[1],
        constant * this->data[2]
      }
    };
  }

  constexpr Vector operator / (const double& constant) const {
    return Vector {
      {
        this->data[0] / constant,
        this->data[1] / constant,
        this->data[2] / constant
      }
    };
  }
};

constexpr double angle(const Vector& a, const Vector& b) {
  return Math::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

constexpr double toDegrees(const double& radians) {
  return 180 * radians / M_PI;
}

} // namespace ConstexprMagic

#endif
