#ifndef INCLUDE_MOLECULE_MANIP_DISTANCE_GEOMETRY_ERROR_H
#define INCLUDE_MOLECULE_MANIP_DISTANCE_GEOMETRY_ERROR_H

#include <system_error>

enum class DGError {
  ZeroPermutationStereocenters = 1,
  GraphImpossible = 2,
  TooManyFailures = 3
};

// Boilerplate to allow interoperability of DGerror with std::error_code
namespace std {
  template<> struct is_error_code_enum<DGError> : std::true_type {};
} // namespace std

namespace detail {
  struct DGError_category : public std::error_category {
    virtual const char* name() const noexcept override final {
      return "DistanceGeometryError";
    }

    virtual std::string message(int c) const override final {
      switch(static_cast<DGError>(c)) {
        case DGError::ZeroPermutationStereocenters:
          return "Graph contains Stereocenters with zero possible permutations.";
        case DGError::GraphImpossible:
          return "Graph cannot be modeled in three-dimensional space.";
        case DGError::TooManyFailures:
          return "Refinement yielded too many failures.";
        default:
          return "Unknown error.";
      };
    }
  };
} // namespace detail

#define THIS_MODULE_API_DECL extern inline

THIS_MODULE_API_DECL const detail::DGError_category& DGError_category() {
  static detail::DGError_category c;
  return c;
}

inline std::error_code make_error_code(DGError e) {
  return {static_cast<int>(e), DGError_category()};
}

#endif
