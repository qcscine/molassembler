#ifndef LIB_SYMMETRY_INFORMATION_HPP
#define LIB_SYMMETRY_INFORMATION_HPP

#ifndef UNUSED
#define UNUSED(x) (void)(x)
#endif

#include <vector>
#include <functional>
#include <cassert>

/* TODO
 * - conditional rotations?
 */

namespace PermSymmetry {

struct SymmetryInformation {
  /*static const unsigned size;

  static std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > rotations;

  static double angle(
    const uint8_t& a,
    const uint8_t& b
  );
  */

  /* TODO 
   * - more members to extract embedding and refinement constraints 
   */
};

/* 4 */
template<typename T>
struct Tetrahedral : public SymmetryInformation {
  /* Positions are enumerated as
   *
   *    2
   *    |
   *   (1) (1 is on top, ( ) signifies the central atom 
   *  /   \
   * 3     4
   *
   * Remember Newman projections? This is sort of supposed to be that.
   *
   * Alternatively:
   *
   *    1
   *    |
   *   (_)   
   *  /  \ °4
   * 2    3 
   *
   * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
   *
   */
  static constexpr unsigned size = 4;

  static const std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > rotations;

  static double angle(
    const uint8_t& a,
    const uint8_t& b
  ) {
    /* Signal unused on purpose */
    UNUSED(a);
    UNUSED(b);
    return 109.5;
  }
};

template<typename T>
struct SquarePlanar : public SymmetryInformation {
  /* Positions are enumerated as
   *
   * 4   3
   *  \_/
   *  (_) <- central atom
   *  / \
   * 1   2
   *
   */
  static constexpr unsigned size = 4;

  static const std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > rotations;

  static double angle(
    const uint8_t& a,
    const uint8_t& b
  ) {
    assert(a < size && b < size && a != b);
    if((a + b) % 2 == 1) {
      // this expression indicates cis
      return 90;
    } else {
      // are trans
      return 180;
    }
  }
};

/* 5 */
template<typename T>
struct SquarePyramidal : public SymmetryInformation {
  /* Positions are enumerated as
   *
   * 4   3
   *  \_/
   *  (5)   
   *  / \
   * 1   2
   *
   * Viewed from the top of the pyramid. The central atom is ( ), 5 is axial.
   *
   * Alternatively
   *
   *    5
   *    |
   * 4  |  3
   *  : | :    <- behind view plane
   *   (_)
   *  /   \    <- in front of view plane
   * 1     2
   *
   */
  static constexpr unsigned size = 5;

  static const std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > rotations;

  static double angle(
    const uint8_t& a,
    const uint8_t& b
  ) {
    assert(a < size && b < size && a != b && a < b);
    if(b < 4 && (a + b) % 2 == 0) {
      /* as long as b is not the pyramidal ligand, 
       * (a + b) % 2 == 0 is true for trans ligands only
       */
      return 180;
    } else { 
      // all other cases are
      return 90;
    }
  }
};

template<typename T>
struct TrigonalBiPyramidal : public SymmetryInformation {
  /* Positions are enumerated as
   *
   * Viewed from the top of the pyramid. The central atom is ( ), 4 and 5 
   * are axial.
   *
   *     4
   *     |  3
   *     | :    <- behind view plane
   * 1--(_)
   *     | \    <- in front of view plane
   *     |  2
   *     5
   */
  static constexpr unsigned size = 5;

  static const std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > rotations;

  static const std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > pseudorotations;

  static double angle(
    const uint8_t& a,
    const uint8_t& b
  ) {
    assert(a < size && b < size && a != b && a < b);
    if(b < 3) {
      // since a < b, this means either 0,1 0,2 1,2 
      return 120;
    } else if(b == 3) {
      // since a < b, this means {1,2,3}, 3 
      return 90;
    } else if(a < 3) {
      // now, b must be 4, and if a is not 3, then
      return 90;
    } else {
      // only case left: 3,4
      return 180;
    }
  }

};


/* 6 */
template<typename T>
struct Octahedral : public SymmetryInformation {
  /* Positions are enumerated as
   *
   * Viewed from the top of the pyramid. The central atom is ( ), 5 and 6 
   * are axial.
   *
   *     5
   *  4  |  3
   *   : | :
   *    (_)        
   *   / | \
   *  1  |  2
   *     6
   *
   * Where /, \ denote bonds in front of the view plane, : denotes bonds
   * behind the view plane.
   */
  static constexpr unsigned size = 6;

  static const std::vector<
    std::pair<
      std::function<
        std::vector<T>(
          const std::vector<T>&
        )
      >,
      uint8_t
    >
  > rotations;

  static double angle(
    const uint8_t& a,
    const uint8_t& b
  ) {
    assert(a < size && b < size && a != b && a < b);
    if(
      (
        b < 4 // if b < 4, then only equatorial ligands
        && (a + b) % 2 == 0 // this gives trans eq ligands
      ) || a == 4 // this indicates 4,5 (axial trans)
    ) {
      return 180;
    } else {
      return 90;
    }
  }
};

} // eo namespace PermSymmetry

#include "SymmetryInformation.hxx"

#endif
