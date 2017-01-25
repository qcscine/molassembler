#ifndef LIB_SYMMETRY_INFORMATION_H
#define LIB_SYMMETRY_INFORMATION_H

#include <functional>
#include <cassert>
#include <array>

/* TODO
 * - conditional rotations?
 */

namespace PermSymmetry {

struct SymmetryInformation {
  /*static const unsigned size;

  static const std::array<
    std::array<unsigned, _size_>,
    _number_of_rotations_
  > rotations;

  static double constexpr angle(
    const unsigned& a,
    const unsigned& b
  );
  */
};

/* 2 */
struct Linear : public SymmetryInformation {
  /* 1 – (_) – 2 */
  static constexpr unsigned size = 2;

  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{ // c++11 quirk
    {{1, 0}}
  }};

  static double constexpr angle(
    const unsigned& a __attribute__((unused)),
    const unsigned& b __attribute__((unused))
  ) {
    return 180;
  }
};

/* NOTE: There should be Bent (2), but since the angle is completely 
 * unspecifiable, this has not been added.
 */

/* 3 */
struct TrigonalPlanar : public SymmetryInformation {
  /* Positions are enumerated as
   *
   *     1
   *     |
   *    (_)
   *   /   \
   *  2     3
   *
   * This is not quite ideal since the angles are thoroughly misrepresented, 
   * but all positions including the central atom are in one plane. The angles
   * are idealized as 120°.
   */

  static constexpr unsigned size = 3;

  static constexpr std::array<
    std::array<unsigned, 3>,
    2
  > rotations {{ //c++11 quirk
    {{1, 2, 0}}, // C3
    {{0, 2, 1}} // C2
  }};

  static double constexpr angle(
    const unsigned& a __attribute__((unused)),
    const unsigned& b __attribute__((unused))
  ) {
    return 120;
  }
};

/* 4 */
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

  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > rotations {{
    {{0, 3, 1, 2}},
    {{2, 1, 3, 0}},
    {{3, 0, 2, 1}},
    {{1, 2, 0, 3}}
  }};

  static double constexpr angle(
    const unsigned& a __attribute__((unused)),
    const unsigned& b __attribute__((unused))
  ) {
    return 109.5;
  }
};

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

  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > rotations {{
    {{3, 0, 1, 2}}, // C4
    {{1, 0, 3, 2}}, // C2
    {{3, 2, 1, 0}}  // C2'
  }};

  static double constexpr angle(
    const unsigned& a,
    const unsigned& b
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

  static constexpr std::array<
    std::array<unsigned, 5>,
    1
  > rotations {{
    {{3, 0, 1, 2, 4}} // C4
  }};

  static double constexpr angle(
    const unsigned& a,
    const unsigned& b
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

  static constexpr std::array<
    std::array<unsigned, 5>,
    4
  > rotations {{
    {{2, 0, 1, 3, 4}}, // C3
    {{0, 2, 1, 4, 3}}, // C2 on 1
    {{2, 1, 0, 4, 3}}, // C2 on 2
    {{1, 0, 2, 4, 3}} // C2 on 3
  }};

  static constexpr std::array<
    std::array<unsigned, 5>,
    1
  > pseudorotations {{
    {{0, 4, 3, 2, 1}} // Pseudorotation with 1 fixed combines with C3 for all
  }};

  static double constexpr angle(
    const unsigned& a,
    const unsigned& b
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

  static constexpr std::array<
    std::array<unsigned, 6>,
    3
  > rotations {{
    {{3, 0, 1, 2, 4, 5}}, // vertical C4
    {{0, 5, 2, 4, 1, 3}}, // horizontal C4
    {{4, 1, 5, 3, 2, 0}} // horizontal C4'
  }};

  static double constexpr angle(
    const unsigned& a,
    const unsigned& b
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

#endif
