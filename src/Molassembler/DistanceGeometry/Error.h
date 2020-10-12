/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Defines DG result type error categories and strings for boost::outcome
 *
 * Contains the DG error_category definitions for use with boost::outcome
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_ERROR_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_ERROR_H

#include <system_error>

/**
 * @brief Errors that can occur in Distance Geometry related algorithms.
 */
enum class DgError {
  /**
   * @brief The molecule you are trying to generate conformers for has
   *   zero-assignment stereopermutators, meaning that it is not representable
   *   in three dimensional space.
   *
   * AtomStereopermutators remove those stereopermutations from the
   * user-accessible set that it deems obviously impossible. This includes
   * overlapping haptic ligand binding cones and multidentate ligand bridges
   * that are too short to span the angle needed in a stereopermutation. In most
   * cases, this will simply eliminate trans-arranged multidentate ligands with
   * too short bridges. It is conservative, however, and may not be strict
   * enough to eliminate all stereopermutators with bridges that are too short
   * to comply with the spatial modeling. In that case, you may get
   * GraphImpossible.
   *
   * If you get this error, reconsider whether your input graph is reasonable
   * and representable in three dimensions. If you believe the atom
   * stereopermutator incorrectly has zero assignments, please contact us and
   * open an issue.
   */
  ZeroAssignmentStereopermutators = 1,
  /**
   * @brief The molecule you are trying to generate conformers for is either
   *   incompatible with the applied spatial model or plain not representable
   *   in three dimensions.
   *
   * The applied spatial model is not very smart and mostly applies simple
   * geometric considerations. One one hand, it may be that centers whose
   * shapes are heavily distorted due to e.g. multiple small cycles are not
   * recognized correctly or modeled loosely enough in order for a conformer to
   * be possible. On the other hand, it is also possible to create graphs that
   * are not representable in three dimensions. In both circumstances, you will
   * get this error to indicate that the spatial model cannot deal with your
   * input graph.
   *
   * If you get this error, reconsider whether your input graph is reasonable
   * and representable in three dimensions. If you are sure it is, please
   * contact us and open an issue.
   */
  GraphImpossible = 2,
  /**
   * @brief An exception occurred during refinement
   *
   * The form of the potential during Distance Geometry refinement can be
   * exceptional due to e.g. divisions by zero. These exceptions can occur
   * completely randomly and have no bearing on validity of input.
   *
   * If you get this error, generate some more conformers. If all of your inputs
   * yield refinement exceptions, there might be a modeling problem, so please
   * contact us and open an issue.
   */
  RefinementException = 3,
  /**
   * @brief Refinement could not find a minimum in your specified maximum
   *   number of iterations
   *
   * Typically, this may mean that the Molecule you are trying to generate
   * conformers is either way too big for the number of iterations in the
   * potential minimization or that refinement got stuck in some sort of bad
   * back and forwards.
   *
   * Try adjusting your number of iterations. If the problem persists, there
   * may be a problem with the form of the refinement potential, so please
   * contact us and open an issue.
   */
  RefinementMaxIterationsReached = 4,
  /**
   * @brief The result of a refinement did not meet criteria for acceptance
   *
   * In Distance Geometry, we generate a list of atom-pairwise distance bounds
   * that indicate the minimum and maximum distance that atom pair should have
   * in a final conformation. Additionally, chiral constraints are generated
   * (similar to improper dihedrals) that indicate whether a chiral element is
   * arranged correctly. Refinement, which tries to minimize these errors, may
   * end up in a local minimum that still violates some of these bounds.
   *
   * This is purely a stochastic problem and should not reflect on your inputs.
   * If you get this error, generate some more conformers. If all of your
   * inputs yield this error, there may be a problem with the refinement
   * potential, so please contact us and open an issue.
   */
  RefinedStructureInacceptable = 5,
  /**
   * @brief Chiral constraints on the refined structure are still incorrect
   *
   * After refinement, chiral constraints are completely wrong. This is a rare
   * stochastic error and should not reflect on your inputs.
   *
   * Generate more conformers. If you get this error a lot on your inputs,
   * there may be a problem with the refinement potential, so please contact us
   * and open an issue.
   */
  RefinedChiralsWrong = 6,
  /**
   * @brief In directed conformer generation, failed to generate decision list
   */
  DecisionListMismatch = 7,
  /**
   * @brief Unknown exception
   */
  UnknownException = 8
};

// Boilerplate to allow interoperability of DgError with std::error_code
namespace std {
  template<> struct is_error_code_enum<DgError> : std::true_type {};
} // namespace std

namespace Detail {
  struct DGError_category : public std::error_category {
    virtual const char* name() const noexcept override final {
      return "DistanceGeometryError";
    }

    virtual std::string message(int c) const override final {
      switch(static_cast<DgError>(c)) {
        case DgError::ZeroAssignmentStereopermutators:
          return "Graph contains Stereopermutators with zero possible permutations.";
        case DgError::GraphImpossible:
          return "Graph cannot be modeled in three-dimensional space.";
        case DgError::RefinementException:
          return "Refinement encountered an exception during minimization.";
        case DgError::RefinementMaxIterationsReached:
          return "Refinement did not converge in maximum number of iterations.";
        case DgError::RefinedStructureInacceptable:
          return "Refined structure deemed inacceptable.";
        case DgError::RefinedChiralsWrong:
          return "Refined structure has chiral constraints with wrong sign.";
        case DgError::DecisionListMismatch:
          return "Failed to generate decision list.";
        case DgError::UnknownException:
          return "Conformer generation encountered an unexpected exception.";
        default:
          return "Unknown error.";
      };
    }
  };
} // namespace Detail

extern inline const Detail::DGError_category& DGError_category() {
  static Detail::DGError_category c;
  return c;
}

inline std::error_code make_error_code(DgError e) {
  return {static_cast<int>(e), DGError_category()};
}

#endif
