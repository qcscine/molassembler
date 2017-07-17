#include "Symmetries.h"

#include <cassert>
#include <algorithm>
#include <Eigen/Geometry>

#include "constexpr_magic/Math.h"

namespace Symmetry {

const std::map<Name, SymmetryInformation> symmetryData {
  {
    Name::Linear, 
    SymmetryInformation {
      "linear",
      2,
      RotationsList {
        {1, 0}
      },
      TetrahedronList {},
      CoordinateList {
        { 1 , 0, 0 },
        { -1, 0, 0 }
      }
    },
  },
  {
    Name::Bent,
    SymmetryInformation {
      /*
       *  1
       *   \
       *    (_) – 0
       *
       */
      "bent",
      2,
      RotationsList {
        {1, 0}
      },
      TetrahedronList {},
      CoordinateList {
        {1, 0, 0},
        Eigen::AngleAxisd(
          M_PI * 107 / 180,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d({1, 0, 0})
      }
    }
  },
  {
    Name::TrigonalPlanar,
    SymmetryInformation {
      /*
       *     0
       *     |
       *    (_)
       *   /   \
       *  1     2
       *
       * This is not quite ideal since the angles are thoroughly misrepresented, 
       * but all positions including the central atom are in one plane. The angles
       * are idealized as 120°.
       */
      "trigonal planar",
      3,
      RotationsList {
        {1, 2, 0}, // C3
        {0, 2, 1} // C2
      },
      TetrahedronList {},
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d({1, 0, 0}),
        Eigen::AngleAxisd(
          - 2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d({1, 0, 0}),
      }
    }
  }, 
  {
    Name::TrigonalPyramidal,
    SymmetryInformation {
      /*
       *    
       *   (_)   
       *  : | :
       * 0  1  2
       *
       */
      "trigonal pyramidal",
      3,
      RotationsList {
        {2, 0, 1} // C3
      },
      TetrahedronList {
        {{boost::none, 0, 1, 2}}
      },
      CoordinateList {
        Eigen::AngleAxisd(
          M_PI * 111.5 / 180,
          Eigen::Vector3d::UnitX()
        ) * Eigen::Vector3d::UnitY(),
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitY()
        ) * Eigen::AngleAxisd(
          M_PI * 111.5 / 180,
          Eigen::Vector3d::UnitX()
        ) * Eigen::Vector3d::UnitY(),
        Eigen::AngleAxisd(
          - 2 * M_PI / 3,
          Eigen::Vector3d::UnitY()
        ) * Eigen::AngleAxisd(
          M_PI * 111.5 / 180,
          Eigen::Vector3d::UnitX()
        ) * Eigen::Vector3d::UnitY(),
      }
    }
  }, 
  {
    Name::TShaped,
    SymmetryInformation {
      /*
       * 0 – (_) – 2
       *      |
       *      1
       *
       */
      "T-shaped",
      3,
      RotationsList {
        {2, 1, 0} // C2
      },
      TetrahedronList {},
      CoordinateList {
        - Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitY(),
        Eigen::Vector3d::UnitX()
      }
    }
  }, 
  {
    Name::Tetrahedral,
    SymmetryInformation {
      /* 
       *    1
       *    |
       *   (0) (0 is on top, ( ) signifies the central atom  beneath it
       *  /   \
       * 2     3
       *
       * Remember Newman projections? This is sort of supposed to be that.
       *
       * Alternatively:
       *
       *    0
       *    |
       *   (_)   
       *  /  \ °3
       * 1    2 
       *
       * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
       *
       */
      "tetrahedral",
      4,
      RotationsList {
        {0, 3, 1, 2}, // C4, 1
        {2, 1, 3, 0}, // C4, 2
        {3, 0, 2, 1}, // C4, 3
        {1, 2, 0, 3}  // C4, 4
      },
      TetrahedronList {
        {{0, 1, 2, 3}}
      },
      CoordinateList {
        Eigen::Vector3d::UnitY(),
        Eigen::AngleAxisd(
          M_PI * 109.5 / 180,
          Eigen::Vector3d::UnitX()
        ) * Eigen::Vector3d::UnitY(),
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitY()
        ) * Eigen::AngleAxisd(
          M_PI * 109.5 / 180,
          Eigen::Vector3d::UnitX()
        ) * Eigen::Vector3d::UnitY(),
        Eigen::AngleAxisd(
          - 2 * M_PI / 3,
          Eigen::Vector3d::UnitY()
        ) * Eigen::AngleAxisd(
          M_PI * 109.5 / 180,
          Eigen::Vector3d::UnitX()
        ) * Eigen::Vector3d::UnitY(),
      }
    }
  }, 
  {
    Name::SquarePlanar,
    SymmetryInformation {
      /* 
       * 3   2
       *  \_/
       *  (_) <- central atom
       *  / \
       * 0   1
       *
       */
      "square planar",
      4,
      RotationsList {
        {3, 0, 1, 2}, // C4
        {1, 0, 3, 2}, // C2
        {3, 2, 1, 0}  // C2'
      },
      TetrahedronList {},
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitY(),
        - Eigen::Vector3d::UnitX(),
        - Eigen::Vector3d::UnitY()
      }
    }
  }, 
  {
    Name::Seesaw,
    SymmetryInformation {
      /*
       * 0 – (_) – 3
       *     / :
       *    1   2
       *
       */
      "seesaw",
      4,
      RotationsList {
        {3, 2, 1, 0} // C2
      },
      TetrahedronList {
#ifdef USE_ALTERNATE_TETRAHEDRA
        // Alternate
        {{0, 1, 2, 3}}
#else
        // Regular
        {{0, boost::none, 1, 2}},
        {{boost::none, 3, 1, 2}},
#endif
      },
      CoordinateList {
        Eigen::Vector3d::UnitY(),
        Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
        - Eigen::Vector3d::UnitY()
      }
    }
  }, 
  {
    Name::SquarePyramidal,
    SymmetryInformation {
      /* 
       * 3   2
       *  \_/
       *  (4)   
       *  / \
       * 0   1
       *
       * Viewed from the top of the pyramid. The central atom is ( ), 5 is axial.
       *
       * Alternatively,
       *
       *    4
       * 3  |  2
       *  : | :    <- behind view plane
       *   (_)
       *  /   \    <- in front of view plane
       * 0     1
       *
       */
      "square pyramidal",
      5,
      RotationsList {
        {3, 0, 1, 2, 4} // C4
      },
      TetrahedronList {
#ifdef USE_ALTERNATE_TETRAHEDRA
        // Alternate
        {{0, 1, 4, 2}}, 
        {{0, 3, 2, 4}}
#else 
        // Regular
        {{0, 1, 4, boost::none}},
        {{1, 2, 4, boost::none}},
        {{2, 3, 4, boost::none}},
        {{3, 0, 4, boost::none}}
#endif
      },
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitY(),
        - Eigen::Vector3d::UnitX(),
        - Eigen::Vector3d::UnitY(),
        Eigen::Vector3d::UnitZ()
      }
    }
  }, 
  {
    Name::TrigonalBiPyramidal,
    SymmetryInformation {
      /* Viewed from the top of the pyramid. The central atom is ( ), 3 and 4 
       * are axial.
       *
       *     3
       *     |  2
       *     | :    <- behind view plane
       * 0--(_)
       *     | \    <- in front of view plane
       *     |  1
       *     4
       */
      "trigonal bipyramidal",
      5,
      RotationsList {
        {2, 0, 1, 3, 4}, // C3
        {0, 2, 1, 4, 3}, // C2 on 0
        {2, 1, 0, 4, 3}, // C2 on 1
        {1, 0, 2, 4, 3} // C2 on 2
      },
      TetrahedronList {
        {{0, 1, 3, 2}},
        {{0, 1, 2, 4}}
      },
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          - 2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitZ(),
        - Eigen::Vector3d::UnitZ()
      }
    }
  }, 
  {
    Name::PentagonalPlanar,
    SymmetryInformation {
      /* 
       * All in plane:
       *
       *      0
       *  1.  |  .4
       *    °(_)°
       *    /   \
       *   2     3
       *
       */
      "pentagonal planar",
      5,
      RotationsList {
        {4, 0, 1, 2, 3}, // C5
        {0, 4, 3, 2, 1} // C2
      },
      TetrahedronList {},
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          3 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          4 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX()
      }
    }
  }, 
  {
    Name::Octahedral,
    SymmetryInformation {
    /* The central atom is ( ), 4 and 5 are axial, the rest equatorial.
     *
     *     4
     *  3  |  2
     *   : | :
     *    (_)        
     *   / | \
     *  0  |  1
     *     5
     *
     * Where /, \ denote bonds in front of the view plane, : denotes bonds
     * behind the view plane.
     */
      "octahedral",
      6,
      RotationsList {
        {3, 0, 1, 2, 4, 5}, // vertical C4
        {0, 5, 2, 4, 1, 3}, // horizontal C4
        {4, 1, 5, 3, 2, 0} // horizontal C4'
      },
      TetrahedronList {
#ifdef USE_ALTERNATE_TETRAHEDRA
        // Alternate
        {{3, 0, 4, 5}},
        {{0, 1, 4, 5}},
        {{1, 2, 4, 5}},
        {{2, 3, 4, 5}}
#else
        // Regular
        {{3, 0, 4, boost::none}},
        {{0, 1, 4, boost::none}},
        {{1, 2, 4, boost::none}},
        {{2, 3, 4, boost::none}},
        {{3, 0, boost::none, 5}},
        {{0, 1, boost::none, 5}},
        {{1, 2, boost::none, 5}},
        {{2, 3, boost::none, 5}}
#endif
      },
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitY(),
        - Eigen::Vector3d::UnitX(),
        - Eigen::Vector3d::UnitY(),
        Eigen::Vector3d::UnitZ(),
        - Eigen::Vector3d::UnitZ()
      }
    }
  },
  {
    Name::TrigonalPrismatic,
    SymmetryInformation {
      /* 
       *  3  4  5
       *   : | :
       *    (_)        
       *   : | :
       *  0  1  2
       *
       * Where /, \ denote bonds in front of the view plane, : denotes bonds
       * behind the view plane.
       *
       * Angles 
       *  0-1, 0-2 -> 86°
       *  0-3 -> 76°
       *  0-4, 0-5 -> 134°
       * 
       * From [W(CH3)6], Haaland, Hammel, Rypdal, Volden, J. Am. Chem. Soc. 1990 
       */
      "trigonal prismatic",
      6,
      RotationsList {
        {2, 0, 1, 5, 3, 4}, // C3 axial
        {5, 4, 3, 2, 1, 0} // C2 betw. 1, 4
      },
      TetrahedronList {
        // TODO dubious if this captures all relevant information, too limited
        {{boost::none, 0, 1, 2}},
        {{3, boost::none, 4, 5}}
      },
      CoordinateList {
        // 0, lower X by 76/2° into -Z
        Eigen::AngleAxisd(
          M_PI * 38 / 180,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
        // 1, take 0 and rotate 2 * pi / 3 around Z
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::AngleAxisd(
          M_PI * 38 / 180,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
        // 2, take 0 and rotate - 2 * pi / 3 around Z
        Eigen::AngleAxisd(
          - 2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::AngleAxisd(
          M_PI * 38 / 180,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
        // 3, raise X by 76/2° into Z
        Eigen::AngleAxisd(
          - M_PI * 38 / 180,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
        // 4, take 3 and rotate 2 * pi / 3 around Z
        Eigen::AngleAxisd(
          2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::AngleAxisd(
          - M_PI * 38 / 180,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
        // 5, take 3 and rotate - 2 * pi / 3 around Z
        Eigen::AngleAxisd(
          - 2 * M_PI / 3,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::AngleAxisd(
          - M_PI * 38 / 180,
          Eigen::Vector3d::UnitY()
        ) * Eigen::Vector3d::UnitX(),
      }
    }
  },
  {
    Name::PentagonalPyramidal,
    SymmetryInformation {
      /* 
       *      0
       *  1.  |  .4
       *    °(5)°
       *    :   :
       *   2     3
       *
       * 0-4 are in plane,
       * 5 is above plane,
       * ( ) signifies the central atom beneath it
       */
      "pentagonal pyramidal",
      6,
      RotationsList {
        {4, 0, 1, 2, 3, 5} // C5 axial
      },
      TetrahedronList {
#ifdef USE_ALTERNATE_TETRAHEDRA
        // Alternate
        {{0, 1, 5, 2}},
        {{2, 3, 5, 4}},
        {{4, 5, boost::none, 0}}
#else
        // Regular
        {{0, boost::none, 1, 5}},
        {{1, boost::none, 2, 5}},
        {{2, boost::none, 3, 5}},
        {{3, boost::none, 4, 5}},
        {{4, boost::none, 0, 5}}
#endif
      },
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          3 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          4 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitZ()
      }
    }
  },
  {
    Name::PentagonalBiPyramidal,
    SymmetryInformation {
      /* 
       * 3, 5, (_) and 6 in plane, 1 and 2 in front, 0 and 4 in back
       *
       *      5
       *  0_  | .4
       *   _:(_) – 3
       *  1   |°·2
       *      6
       *
       * 0-4 are equatorial, 
       * 5,6 are axial
       */
      "pentagonal bipyramidal",
      7,
      RotationsList {
        {4, 0, 1, 2, 3, 5, 6}, // C5 axial
        {1, 0, 4, 3, 2, 6, 5} // C2 equatorial on 3
      },
      TetrahedronList {
#ifdef USE_ALTERNATE_TETRAHEDRA
        // Alternate
        {{0, 1, 5, 6}},
        {{1, 2, 5, 6}},
        {{2, 3, 5, 6}},
        {{3, 4, 5, 6}},
        {{4, 0, 5, 6}}
#else
        // Regular
        {{0, 1, 5, boost::none}},
        {{1, 2, 5, boost::none}},
        {{2, 3, 5, boost::none}},
        {{3, 4, 5, boost::none}},
        {{4, 0, 5, boost::none}},
        {{0, 1, boost::none, 6}},
        {{1, 2, boost::none, 6}},
        {{2, 3, boost::none, 6}},
        {{3, 4, boost::none, 6}},
        {{4, 0, boost::none, 6}}
#endif
      },
      CoordinateList {
        Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          2 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          3 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::AngleAxisd(
          4 * 2 * M_PI / 5,
          Eigen::Vector3d::UnitZ()
        ) * Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::UnitZ(),
        - Eigen::Vector3d::UnitZ()
      }
    }
  },
  {
    Name::SquareAntiPrismatic,
    SymmetryInformation {
      /*  Two representations, one oblique, the other more helpful.
       *  The first is a side-on view. 0 is mostly hidden by 1, 3 mostly hidden 
       *  by 2. 4 and 6 are in the viewing plane while 5 juts out above plane and
       *  7 dips behind plane. 
       *
       *  4   7 5 6
       *    : ·/ :
       *     (__)
       *    ·/  \·
       *   01    23 
       *
       * Below is a top-down view. Strong lines indicate above-plane bonds, dots
       * indicate below-plane bonds.
       *
       *   0   7   3
       *     · | ·
       *   4 –( )– 6
       *     · | ·
       *   1   5   2
       *
       * Angles:
       *
       *   in-plane cis (4, 5) -> 55°
       *   in-plane trans (4, 6) -> 148°
       *   short between planes (0, 4) -> 51°
       *   long between planes (0, 5) -> 175°
       *
       * from [ReF8]2-, Koz'min, P.A., Zhurnal Strukturnoi Khimii 1964 
       * HINT: use the ICSD for angle calculations and preview
       *
       */
      "square antiprismatic",
      8,
      RotationsList {
        {3, 0, 1, 2, 7, 4, 5, 6}, // C4 axial
        /* 180° on equatorial axis in plane with 4, 6 
         *
         *   1   5   2
         *     \ · /
         * – 4 ·( )· 6 –––––– equatorial axis
         *     / · \
         *   0   7   3
         *
         * and 45° anticlockwise on axial axis
         *           
         *   5   2   6
         *     · | ·
         *   1 –( )– 3
         *     · | ·
         *   4   0   7
         *
         */
        {5, 4, 7, 6, 1, 0, 3, 2}, 
      },
      TetrahedronList {
#ifdef USE_ALTERNATE_TETRAHEDRA
        // Alternate
        {{0, 1, 4, 6}},
        {{1, 2, 5, 7}},
        {{2, 3, 6, 4}},
        {{3, 0, 7, 5}},
#else 
        // Regular
        {{7, 0, 4, boost::none}},
        {{0, 4, boost::none, 1}},
        {{4, 1, 5, boost::none}},
        {{1, 5, boost::none, 2}},
        {{5, 2, 6, boost::none}},
        {{2, 6, boost::none, 3}},
        {{6, 3, 7, boost::none}},
        {{3, 7, boost::none, 0}}
#endif
      },
      // TODO not quite ideal, not oriented along any axis...
      CoordinateList { // generated w/ reference/minimal.py
        {-0.00928803, 0.61156848, 0.79113698},
        {0.79562737, 0.60564101, -0.01326839},
        {0.79562737, -0.60564101, -0.01326839},
        {-0.00928803, -0.61156848, 0.79113698},
        {-0.3961716, 0.85216935, -0.34184129},
        {0.29375817, 0., -0.95587977},
        {-0.3961716, -0.85216935, -0.34184129},
        {-0.98308669, 0., 0.18314084}
      }
    }
  }
};

constexpr unsigned replaceMe = std::numeric_limits<unsigned>::max();

namespace data {

struct Linear {
  static constexpr Symmetry::Name name = Symmetry::Name::Linear;
  static constexpr unsigned size = 2;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return M_PI;
  }
  static constexpr std::array<ConstexprMagic::Vector, 2> coordinates {{
    { 1 , 0, 0 },
    { -1, 0, 0 }
  }};
  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{
    {1, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

struct Bent {
  /*
   *  1
   *   \
   *    (_) – 0
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Bent;
  static constexpr unsigned size = 2;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    /* subject to a lot of variation, between 90 and 109 degrees pursuant to 
     * english wikipedia, using experimental data here to improve instances
     * of this geometry on e.g. O center would a big improvement to DG runs
     */
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(107); 
  }
  static constexpr std::array<ConstexprMagic::Vector, 2> coordinates {{
    {1., 0., 0.},
    {-0.292372, 0.956305, 0.}
  }};
  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{
    {1, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

struct TrigonalPlanar {
  /*
   *     0
   *     |
   *    (_)
   *   /   \
   *  1     2
   *
   * This is not quite ideal since the angles are thoroughly misrepresented, 
   * but all positions including the central atom are in one plane. The angles
   * are idealized as 120°.
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPlanar;
  static constexpr unsigned size = 3;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(120);
  }
  static constexpr std::array<ConstexprMagic::Vector, 3> coordinates {{
      {1, 0, 0},
      {-0.5, 0.866025, 0},
      {-0.5, -0.866025, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    2
  > rotations {{
    {1, 2, 0}, // C3
    {0, 2, 1} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

struct TrigonalPyramidal {
  /*
   *     0
   *     |
   *    (_)
   *   /   \
   *  1     2
   *
   * This is not quite ideal since the angles are thoroughly misrepresented, 
   * but all positions including the central atom are in one plane. The angles
   * are idealized as 120°.
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPyramidal;
  static constexpr unsigned size = 3;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(107.5);
  }
  static constexpr std::array<ConstexprMagic::Vector, 3> coordinates {{
    {0, -0.366501, 0.930418},
    {0.805765, -0.366501, -0.465209},
    {-0.805765, -0.366501, -0.465209}
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    1
  > rotations {{
    {2, 0, 1} // C3
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{replaceMe, 0, 1, 2}}
  }};
};

struct TShaped {
  /*
   * 0 – (_) – 2
   *      |
   *      1
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TShaped;
  static constexpr unsigned size = 3;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if((a + b) % 2 == 1) {
      return M_PI / 2;
    } 

    return M_PI;
  }
  static constexpr std::array<ConstexprMagic::Vector, 3> coordinates {{
    {-1, -0, -0},
    {0, 1, 0},
    {1, 0, 0},
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    1
  > rotations {{
    {2, 1, 0} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

struct Tetrahedral {
  /* 
   *    1
   *    |
   *   (0) (0 is on top, ( ) signifies the central atom  beneath it
   *  /   \
   * 2     3
   *
   * Remember Newman projections? This is sort of supposed to be that.
   *
   * Alternatively:
   *
   *    0
   *    |
   *   (_)   
   *  /  \ °3
   * 1    2 
   *
   * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Tetrahedral;
  static constexpr unsigned size = 4;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(109.5);
  }
  static constexpr std::array<ConstexprMagic::Vector, 4> coordinates {{
    {0, 1, 0},
    {0, -0.333807, 0.942641},
    {0.816351, -0.333807, -0.471321},
    {-0.816351, -0.333807, -0.471321}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > rotations {{
    {0, 3, 1, 2}, // C4, 1
    {2, 1, 3, 0}, // C4, 2
    {3, 0, 2, 1}, // C4, 3
    {1, 2, 0, 3}  // C4, 4
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{0, 1, 2, 3}}
  }};
};

struct SquarePlanar {
  /* 
   * 3   2
   *  \_/
   *  (_) <- central atom
   *  / \
   * 0   1
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::SquarePlanar;
  static constexpr unsigned size = 4;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if((a + b) % 2 == 1) {
      // this expression indicates cis
      return M_PI / 2;
    } 

    // leftover case is trans
    return M_PI;
  }
  static constexpr std::array<ConstexprMagic::Vector, 4> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > rotations {{
    {3, 0, 1, 2}, // C4
    {1, 0, 3, 2}, // C2
    {3, 2, 1, 0}  // C2'
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

struct Seesaw {
  /*
   * 0 – (_) – 3
   *     / :
   *    1   2
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Seesaw;
  static constexpr unsigned size = 4;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    const auto& smaller = std::min(a, b);
    const auto& larger = std::max(a, b);
    if(smaller == 0 && larger == 3) {
      return M_PI;
    }

    if(smaller == 1 && larger == 2) {
      return ConstexprMagic::Math::toRadians<double>(120);
    }

    return M_PI / 2;
  }
  static constexpr std::array<ConstexprMagic::Vector, 4> coordinates {{
    {0, 1, 0},
    {1, 0, 0},
    {-0.5, 0, -0.866025},
    {-0, -1, -0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > rotations {{
    {3, 2, 1, 0} // C2
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{0, 1, 2, 3}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, replaceMe, 1, 2}},
    {{replaceMe, 3, 1, 2}},
  }};
#endif
};

struct SquarePyramidal {
  /* 
   * 3   2
   *  \_/
   *  (4)   
   *  / \
   * 0   1
   *
   * Viewed from the top of the pyramid. The central atom is ( ), 5 is axial.
   *
   * Alternatively,
   *
   *    4
   * 3  |  2
   *  : | :    <- behind view plane
   *   (_)
   *  /   \    <- in front of view plane
   * 0     1
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::SquarePyramidal;
  static constexpr unsigned size = 5;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if(a == 4 || b == 4) { // all bonds to axial ligand are 90°
      return M_PI / 2; 
    }

    if((a + b) % 2 == 0) { // 0 + 2 or 1 + 3 are trans
      return M_PI;
    }

    // rest are cis
    return M_PI / 2;
  }
  static constexpr std::array<ConstexprMagic::Vector, 5> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0},
    {0, 0, 1}
  }};
  static constexpr std::array<
    std::array<unsigned, 5>,
    1
  > rotations {{
    {3, 0, 1, 2, 4} // C4
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, 1, 4, 2}}, 
    {{0, 3, 2, 4}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
    {{0, 1, 4, replaceMe}},
    {{1, 2, 4, replaceMe}},
    {{2, 3, 4, replaceMe}},
    {{3, 0, 4, replaceMe}}
  }};
#endif
};
struct TrigonalBiPyramidal {
  /* Viewed from the top of the pyramid. The central atom is ( ), 3 and 4 
   * are axial.
   *
   *     3
   *     |  2
   *     | :    <- behind view plane
   * 0--(_)
   *     | \    <- in front of view plane
   *     |  1
   *     4
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalBiPyramidal;
  static constexpr unsigned size = 5;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    unsigned smaller = std::min(a, b), larger = std::max(a, b);
    if(larger < 3) {
      // -> smaller < 2, this means either 0,1 0,2 1,2 axial
      return ConstexprMagic::Math::toRadians<double>(120);
    } else if(larger == 3) {
      // -> smaller < 3, this means {1,2,3}, 3 
      return M_PI / 2;
    } else if(smaller < 3) {
      // now, larger must be 4 (process of elimination), so if a is not 3:
      return M_PI / 2;
    } else {
      // only case left: 3,4
      return M_PI;
    }
  }
  static constexpr std::array<ConstexprMagic::Vector, 5> coordinates {{
    {1, 0, 0},
    {-0.5, 0.866025, 0},
    {-0.5, -0.866025, 0},
    {0, 0, 1},
    {-0, -0, -1}
  }};
  static constexpr std::array<
    std::array<unsigned, 5>,
    4
  > rotations {{
    {2, 0, 1, 3, 4}, // C3
    {0, 2, 1, 4, 3}, // C2 on 0
    {2, 1, 0, 4, 3}, // C2 on 1
    {1, 0, 2, 4, 3} // C2 on 2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
      {{0, 1, 3, 2}},
      {{0, 1, 2, 4}}
  }};
};

struct PentagonalPlanar {
  /* 
   * All in plane:
   *
   *      0
   *  1.  |  .4
   *    °(_)°
   *    /   \
   *   2     3
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPlanar;
  static constexpr unsigned size = 5;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * ConstexprMagic::Math::toRadians<double>(72);
  }
  static constexpr std::array<ConstexprMagic::Vector, 5> coordinates {{
    {1, 0, 0},
    {0.309017, 0.951057, 0},
    {-0.809017, 0.587785, 0},
    {-0.809017, -0.587785, 0},
    {0.309017, -0.951057, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 5>,
    2
  > rotations {{
    {4, 0, 1, 2, 3}, // C5
    {0, 4, 3, 2, 1} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

struct Octahedral {
  /* The central atom is ( ), 4 and 5 are axial, the rest equatorial.
   *
   *     4
   *  3  |  2
   *   : | :
   *    (_)        
   *   / | \
   *  0  |  1
   *     5
   *
   * Where /, \ denote bonds in front of the view plane, : denotes bonds
   * behind the view plane.
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Octahedral;
  static constexpr unsigned size = 6;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }
    
    if(
      (
        std::max(a, b) < 4 // if the largest is < 4, then equatorial 
        && (a + b) % 2 == 0 // this gives trans eq ligands
      ) || std::min(a, b) == 4 // this indicates 4,5 (axial trans)
    ) {
      return M_PI;
    } 

    return M_PI / 2;
  }
  static constexpr std::array<ConstexprMagic::Vector, 6> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0},
    {0, 0, 1},
    {-0, -0, -1}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    3
  > rotations {{
    {3, 0, 1, 2, 4, 5}, // vertical C4
    {0, 5, 2, 4, 1, 3}, // horizontal C4
    {4, 1, 5, 3, 2, 0} // horizontal C4'
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
    {{3, 0, 4, 5}},
    {{0, 1, 4, 5}},
    {{1, 2, 4, 5}},
    {{2, 3, 4, 5}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    8
  > tetrahedra {{
    {{3, 0, 4, replaceMe}},
    {{0, 1, 4, replaceMe}},
    {{1, 2, 4, replaceMe}},
    {{2, 3, 4, replaceMe}},
    {{3, 0, replaceMe, 5}},
    {{0, 1, replaceMe, 5}},
    {{1, 2, replaceMe, 5}},
    {{2, 3, replaceMe, 5}}
  }};
#endif
};

struct TrigonalPrismatic {
  /* 
   *  3  4  5
   *   : | :
   *    (_)        
   *   : | :
   *  0  1  2
   *
   * Where /, \ denote bonds in front of the view plane, : denotes bonds
   * behind the view plane.
   *
   * Angles 
   *  0-1, 0-2 -> 86°
   *  0-3 -> 76°
   *  0-4, 0-5 -> 134°
   * 
   * From [W(CH3)6], Haaland, Hammel, Rypdal, Volden, J. Am. Chem. Soc. 1990 
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPrismatic;
  static constexpr unsigned size = 6;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    } 
    
    // Between plane symmetric
    if(std::min(a - b, b - a) == 3) {
      return ConstexprMagic::Math::toRadians<double>(76);
    } 

    // In plane triangle
    if(
      (a < 3 && b < 3)
      || (a >= 3 && b >= 3)
    ) {
      return ConstexprMagic::Math::toRadians<double>(86);
    } 

    // Between plane asymmetric
    return ConstexprMagic::Math::toRadians<double>(134);
  }
  static constexpr std::array<ConstexprMagic::Vector, 6> coordinates {{
    {0.788011, 0, -0.615661},
    {-0.394005, 0.682437, -0.615661},
    {-0.394005, -0.682437, -0.615661},
    {0.788011, 0, 0.615661},
    {-0.394005, 0.682437, 0.615661},
    {-0.394005, -0.682437, 0.615661}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    2
  > rotations {{
    {2, 0, 1, 5, 3, 4}, // C3 axial
    {5, 4, 3, 2, 1, 0} // C2 betw. 1, 4
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    // TODO dubious if this captures all relevant information, too limited
    {{replaceMe, 0, 1, 2}},
    {{3, replaceMe, 4, 5}}
  }};
};

struct PentagonalPyramidal {
  /* 
   *      0
   *  1.  |  .4
   *    °(5)°
   *    :   :
   *   2     3
   *
   * 0-4 are in plane,
   * 5 is above plane,
   * ( ) signifies the central atom beneath it
   */
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPyramidal;
  static constexpr unsigned size = 6;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if(a == 5 || b == 5) {
      return M_PI / 2;
    }  
    
    // remainder are identical to PentagonalPlanar
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * ConstexprMagic::Math::toRadians<double>(72);
  }
  static constexpr std::array<ConstexprMagic::Vector, 6> coordinates {{
    {1, 0, 0},
    {0.309017, 0.951057, 0},
    {-0.809017, 0.587785, 0},
    {-0.809017, -0.587785, 0},
    {0.309017, -0.951057, 0},
    {0, 0, 1}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    1
  > rotations {{
    {4, 0, 1, 2, 3, 5} // C5 axial
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {{0, 1, 5, 2}},
    {{2, 3, 5, 4}},
    {{4, 5, replaceMe, 0}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    5
  > tetrahedra {{
    {{0, replaceMe, 1, 5}},
    {{1, replaceMe, 2, 5}},
    {{2, replaceMe, 3, 5}},
    {{3, replaceMe, 4, 5}},
    {{4, replaceMe, 0, 5}}
  }};
#endif
};

struct PentagonalBiPyramidal {
    /* 
     * 3, 5, (_) and 6 in plane, 1 and 2 in front, 0 and 4 in back
     *
     *      5
     *  0_  | .4
     *   _:(_) – 3
     *  1   |°·2
     *      6
     *
     * 0-4 are equatorial, 
     * 5,6 are axial
     */
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPyramidal;
  static constexpr unsigned size = 7;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if(a + b == 11) {
      return M_PI; // trans 5,6
    }

    if(ConstexprMagic::Math::XOR(a > 4, b > 4)) {
      return M_PI / 2; // any angle to axial index
    }

    // remainder are equatorial angles, like PentagonalPlanar
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * ConstexprMagic::Math::toRadians<double>(72);
  }
  static constexpr std::array<ConstexprMagic::Vector, 7> coordinates {{
    {1, 0, 0},
    {0.309017, 0.951057, 0},
    {-0.809017, 0.587785, 0},
    {-0.809017, -0.587785, 0},
    {0.309017, -0.951057, 0},
    {0, 0, 1},
    {-0, -0, -1}
  }};
  static constexpr std::array<
    std::array<unsigned, 7>,
    2
  > rotations {{
    {4, 0, 1, 2, 3, 5, 6}, // C5 axial
    {1, 0, 4, 3, 2, 6, 5} // C2 equatorial on 3
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    5
  > tetrahedra {{
      {{0, 1, 5, 6}},
      {{1, 2, 5, 6}},
      {{2, 3, 5, 6}},
      {{3, 4, 5, 6}},
      {{4, 0, 5, 6}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    10
  > tetrahedra {{
    {{0, 1, 5, replaceMe}},
    {{1, 2, 5, replaceMe}},
    {{2, 3, 5, replaceMe}},
    {{3, 4, 5, replaceMe}},
    {{4, 0, 5, replaceMe}},
    {{0, 1, replaceMe, 6}},
    {{1, 2, replaceMe, 6}},
    {{2, 3, replaceMe, 6}},
    {{3, 4, replaceMe, 6}},
    {{4, 0, replaceMe, 6}}
  }};
#endif
};

struct SquareAntiPrismatic {
  /*  Two representations, one oblique, the other more helpful.
   *  The first is a side-on view. 0 is mostly hidden by 1, 3 mostly hidden 
   *  by 2. 4 and 6 are in the viewing plane while 5 juts out above plane and
   *  7 dips behind plane. 
   *
   *  4   7 5 6
   *    : ·/ :
   *     (__)
   *    ·/  \·
   *   01    23 
   *
   * Below is a top-down view. Strong lines indicate above-plane bonds, dots
   * indicate below-plane bonds.
   *
   *   0   7   3
   *     · | ·
   *   4 –( )– 6
   *     · | ·
   *   1   5   2
   *
   * Angles:
   *
   *   in-plane cis (4, 5) -> 55°
   *   in-plane trans (4, 6) -> 148°
   *   short between planes (0, 4) -> 51°
   *   long between planes (0, 5) -> 175°
   *
   * from [ReF8]2-, Koz'min, P.A., Zhurnal Strukturnoi Khimii 1964 
   * HINT: use the ICSD for angle calculations and preview
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::SquareAntiPrismatic;
  static constexpr unsigned size = 8;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
    if(a == b) {
      return 0;
    }

    return squareAntiprismaticAngles.at(
      std::min(a, b),
      std::max(a, b)
    );
#else
    if(a == b) {
      return 0;
    } 
    
    if(
      (a < 4 && b < 4)
      || (a >= 4 && b >= 4)
    ) { // in plane
      if((a + b) % 2 == 1) { // cis
        return ConstexprMagic::Math::toRadians<double>(72.9875); 
      } 
      
      // otherwise trans
      return ConstexprMagic::Math::toRadians<double>(114.475);
    }
    
    // remaining cases are between planes
    unsigned minDiff = std::min(a - b, b - a);
    if(minDiff == 3 || minDiff == 4 || minDiff == 7) { // short
      return ConstexprMagic::Math::toRadians<double>(78.05);
    } 

    // last case is long between planes
    return ConstexprMagic::Math::toRadians<double>(142.275);
#endif
  }
  static constexpr std::array<ConstexprMagic::Vector, 8> coordinates {{
    {-0.00928803, 0.611568, 0.791137},
    {0.795627, 0.605641, -0.0132684},
    {0.795627, -0.605641, -0.0132684},
    {-0.00928803, -0.611568, 0.791137},
    {-0.396172, 0.852169, -0.341841},
    {0.293758, 0, -0.95588},
    {-0.396172, -0.852169, -0.341841},
    {-0.983087, 0, 0.183141}
  }};
  static constexpr std::array<
    std::array<unsigned, 8>,
    2
  > rotations {{
    {3, 0, 1, 2, 7, 4, 5, 6}, // C4 axial
    /* 180° on equatorial axis in plane with 4, 6 
     *
     *   1   5   2
     *     \ · /
     * – 4 ·( )· 6 –––––– equatorial axis
     *     / · \
     *   0   7   3
     *
     * and 45° anticlockwise on axial axis
     *           
     *   5   2   6
     *     · | ·
     *   1 –( )– 3
     *     · | ·
     *   4   0   7
     *
     */
    {5, 4, 7, 6, 1, 0, 3, 2}, 
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
      {{0, 1, 4, 6}},
      {{1, 2, 5, 7}},
      {{2, 3, 6, 4}},
      {{3, 0, 7, 5}},
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    8
  > tetrahedra {{
    {{7, 0, 4, replaceMe}},
    {{0, 4, replaceMe, 1}},
    {{4, 1, 5, replaceMe}},
    {{1, 5, replaceMe, 2}},
    {{5, 2, 6, replaceMe}},
    {{2, 6, replaceMe, 3}},
    {{6, 3, 7, replaceMe}},
    {{3, 7, replaceMe, 0}}
  }};
#endif
};

const std::string Linear::stringName {"linear"};
const std::string Bent::stringName {"bent"};
const std::string TrigonalPlanar::stringName {"trigonal planar"};
const std::string TrigonalPyramidal::stringName {"trigonal pyramidal"};
const std::string TShaped::stringName {"T-shaped"};
const std::string Tetrahedral::stringName {"tetrahedral"};
const std::string SquarePlanar::stringName {"square planar"};
const std::string Seesaw::stringName {"seesaw"};
const std::string TrigonalBiPyramidal::stringName {"trigonal bipyramidal"};
const std::string SquarePyramidal::stringName {"square pyramidal"};
const std::string PentagonalPlanar::stringName {"pentagonal planar"};
const std::string Octahedral::stringName {"octahedral"};
const std::string TrigonalPrismatic::stringName {"trigonal prismatic"};
const std::string PentagonalPyramidal::stringName {"pentagonal pyramidal"};
const std::string PentagonalBiPyramidal::stringName {"pentagonal bipyramidal"};
const std::string SquareAntiPrismatic::stringName {"square antiprismatic"};

} // namespace data

template<typename ... SymmetryClasses>
struct maxSizeCalculator {
  static constexpr unsigned op() {
    const std::array<unsigned, sizeof...(SymmetryClasses)> sizes = {{
      SymmetryClasses::size...
    }};

    // C++17 max_element (is constexpr then)
    unsigned largestSize = 0;

    for(unsigned i = 0; i < sizes.size(); ++i) {
      if(sizes.at(i) > largestSize) {
        largestSize = sizes.at(i);
      }
    }

    return largestSize;
  }
};

template<typename ... SymmetryClasses>
struct maxTetrahedraCalculator {
  static constexpr unsigned op() {
    const std::array<unsigned, sizeof...(SymmetryClasses)> sizes = {{
      SymmetryClasses::tetrahedra.size()...
    }};

    // C++17 max_element (is constexpr then)
    unsigned largestSize = 0;

    for(unsigned i = 0; i < sizes.size(); ++i) {
      if(sizes.at(i) > largestSize) {
        largestSize = sizes.at(i);
      }
    }

    return largestSize;
  }
};

template<typename ... SymmetryClasses>
struct maxRotationsCalculator {
  static constexpr unsigned op() {
    const std::array<unsigned, sizeof...(SymmetryClasses)> sizes = {{
      SymmetryClasses::rotations.size()...
    }};

    // C++17 max_element (is constexpr then)
    unsigned largestSize = 0;

    for(unsigned i = 0; i < sizes.size(); ++i) {
      if(sizes.at(i) > largestSize) {
        largestSize = sizes.at(i);
      }
    }

    return largestSize;
  }
};

using allSymmetryDataTypes = std::tuple<
  data::Linear, // 2
  data::Bent,
  data::TrigonalPlanar, // 3
  data::TrigonalPyramidal,
  data::TShaped,
  data::Tetrahedral, // 4
  data::SquarePlanar,
  data::Seesaw,
  data::SquarePyramidal, // 5
  data::TrigonalBiPyramidal,
  data::PentagonalPlanar,
  data::Octahedral, // 6
  data::TrigonalPrismatic,
  data::PentagonalPyramidal,
  data::PentagonalBiPyramidal, // 7
  data::SquareAntiPrismatic // 8
>;

template<
  typename Tuple,
  template<typename ...> class TemplateFunction,
  std::size_t... I
> constexpr unsigned unpackHelper(std::index_sequence<I...>) {
  return TemplateFunction<
    std::tuple_element_t<
      I,
      Tuple
    >...
  >::op();
}

template<
  typename Tuple,
  template<typename ...> class TemplateFunction
> constexpr unsigned unpackTupleToTemplateFunctor() {
  return unpackHelper<Tuple, TemplateFunction>(
    std::make_index_sequence<
      std::tuple_size<Tuple>::value
    >()
  );
}

constexpr unsigned maxSize = unpackTupleToTemplateFunctor<
  allSymmetryDataTypes,
  maxSizeCalculator
>();

constexpr unsigned maxTetrahedra = unpackTupleToTemplateFunctor<
  allSymmetryDataTypes,
  maxTetrahedraCalculator
>();

constexpr unsigned maxRotations = unpackTupleToTemplateFunctor<
  allSymmetryDataTypes,
  maxRotationsCalculator
>();

struct ConstexprSymmetryInfo {
  using AngleFunctionPtr = double(*)(const unsigned&, const unsigned&);
  using CoordinatesType = std::array<ConstexprMagic::Vector, maxSize>;
  using TetrahedraType = std::array<
    std::array<unsigned, 4>,
    maxTetrahedra
  >;
  using RotationsType = std::array<
    std::array<unsigned, maxSize>,
    maxRotations
  >;

  const Symmetry::Name name;
  const unsigned size;
  const AngleFunctionPtr angleFunction;
  const CoordinatesType coordinates;
  const RotationsType rotations;
  const TetrahedraType tetrahedra;

  constexpr ConstexprSymmetryInfo(
    Symmetry::Name&& name,
    unsigned&& size,
    AngleFunctionPtr angleFunction,
    CoordinatesType&& coordinates,
    RotationsType&& rotations,
    TetrahedraType&& tetrahedra
  ) : name(name),
      size(size),
      angleFunction(angleFunction),
      coordinates(coordinates),
      rotations(rotations),
      tetrahedra(tetrahedra)
  {}
};

/* Maybe everything below is entirely unneeded if I can just iterate through
 * the data types in the tuple instead of instances of Names.
 */

/*template<unsigned long baseSize>
constexpr ConstexprSymmetryInfo::CoordinatesType makeWastefulCoordinates(
  const std::array<ConstexprMagic::Vector, baseSize>& coordinates
) {
  return {
    coordinates.begin(),
    coordinates.end()
  };
}

constexpr auto test = makeWastefulCoordinates(data::Linear::coordinates);*/

/*template<typename SymmetryClass>
constexpr ConstexprSymmetryInfo makeConstexprSymmetryInfo() {
  return {
    SymmetryClass::name,
    SymmetryClass::size,
    SymmetryClass::angleFunction,
    ConstexprSymmetryInfo::CoordinatesType {SymmetryClass::coordinates},
    ConstexprSymmetryInfo::RotationsType {SymmetryClass::rotations},
    ConstexprSymmetryInfo::TetrahedraType {SymmetryClass::tetrahedra}
  };
}

constexpr auto cLinear = makeConstexprSymmetryInfo<data::Linear>();*/

constexpr std::array<ConstexprSymmetryInfo, nSymmetries> constSymmetryData {{
  {
    Symmetry::Name::Linear,
    2,
    &AngleFunctions::linear,
    ConstexprSymmetryInfo::CoordinatesType {{
      { 1 , 0, 0 },
      { -1, 0, 0 }
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {1, 0}
    }},
    ConstexprSymmetryInfo::TetrahedraType {{}}
  },
  {
    Symmetry::Name::Bent,
    2,
    /*
     *  1
     *   \
     *    (_) – 0
     *
     */
    &AngleFunctions::bent,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1., 0., 0.},
      {-0.292372, 0.956305, 0.}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {1, 0}
    }},
    ConstexprSymmetryInfo::TetrahedraType {{}}
  },
  {
    Symmetry::Name::TrigonalPlanar,
    3,
    /*
     *     0
     *     |
     *    (_)
     *   /   \
     *  1     2
     *
     * This is not quite ideal since the angles are thoroughly misrepresented, 
     * but all positions including the central atom are in one plane. The angles
     * are idealized as 120°.
     */
    &AngleFunctions::trigonalPlanar,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {-0.5, 0.866025, 0},
      {-0.5, -0.866025, 0}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {1, 2, 0}, // C3
      {0, 2, 1} // C2
    }},
    ConstexprSymmetryInfo::TetrahedraType {{}}
  },
  {
    Symmetry::Name::TrigonalPyramidal,
    3,
    /*
     *    
     *   (_)   
     *  : | :
     * 0  1  2
     *
     */
    &AngleFunctions::trigonalPyramidal,
    ConstexprSymmetryInfo::CoordinatesType {{
      {0, -0.366501, 0.930418},
      {0.805765, -0.366501, -0.465209},
      {-0.805765, -0.366501, -0.465209}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {2, 0, 1} // C3
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
      {{replaceMe, 0, 1, 2}}
    }}
  },
  {
    Symmetry::Name::TShaped,
    3,
    /*
     * 0 – (_) – 2
     *      |
     *      1
     *
     */
    &AngleFunctions::tShaped,
    ConstexprSymmetryInfo::CoordinatesType {{
      {-1, -0, -0},
      {0, 1, 0},
      {1, 0, 0},
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {2, 1, 0} // C2
    }},
    ConstexprSymmetryInfo::TetrahedraType {{}}
  },
  {
    Symmetry::Name::Tetrahedral,
    4,
    /* 
     *    1
     *    |
     *   (0) (0 is on top, ( ) signifies the central atom  beneath it
     *  /   \
     * 2     3
     *
     * Remember Newman projections? This is sort of supposed to be that.
     *
     * Alternatively:
     *
     *    0
     *    |
     *   (_)   
     *  /  \ °3
     * 1    2 
     *
     * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
     *
     */
    &AngleFunctions::tetrahedral,
    ConstexprSymmetryInfo::CoordinatesType {{
      {0, 1, 0},
      {0, -0.333807, 0.942641},
      {0.816351, -0.333807, -0.471321},
      {-0.816351, -0.333807, -0.471321}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {0, 3, 1, 2}, // C4, 1
      {2, 1, 3, 0}, // C4, 2
      {3, 0, 2, 1}, // C4, 3
      {1, 2, 0, 3}  // C4, 4
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
      {{0, 1, 2, 3}}
    }}
  },
  {
    Symmetry::Name::SquarePlanar,
    4,
    /* 
     * 3   2
     *  \_/
     *  (_) <- central atom
     *  / \
     * 0   1
     *
     */
    &AngleFunctions::squarePlanar,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {0, 1, 0},
      {-1, -0, -0},
      {-0, -1, -0}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {3, 0, 1, 2}, // C4
      {1, 0, 3, 2}, // C2
      {3, 2, 1, 0}  // C2'
    }},
    ConstexprSymmetryInfo::TetrahedraType {{}}
  },
  {
    Symmetry::Name::Seesaw,
    4,
    /*
     * 0 – (_) – 3
     *     / :
     *    1   2
     *
     */
    &AngleFunctions::seesaw,
    ConstexprSymmetryInfo::CoordinatesType {{
      {0, 1, 0},
      {1, 0, 0},
      {-0.5, 0, -0.866025},
      {-0, -1, -0}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {3, 2, 1, 0} // C2
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
#ifdef USE_ALTERNATE_TETRAHEDRA
      // Alternate
      {{0, 1, 2, 3}}
#else
      // Regular
      {{0, replaceMe, 1, 2}},
      {{replaceMe, 3, 1, 2}},
#endif
    }}
  },
  {
    Symmetry::Name::SquarePyramidal,
    5,
    /* 
     * 3   2
     *  \_/
     *  (4)   
     *  / \
     * 0   1
     *
     * Viewed from the top of the pyramid. The central atom is ( ), 5 is axial.
     *
     * Alternatively,
     *
     *    4
     * 3  |  2
     *  : | :    <- behind view plane
     *   (_)
     *  /   \    <- in front of view plane
     * 0     1
     *
     */
    &AngleFunctions::squarePyramidal,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {0, 1, 0},
      {-1, -0, -0},
      {-0, -1, -0},
      {0, 0, 1}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {3, 0, 1, 2, 4} // C4
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
#ifdef USE_ALTERNATE_TETRAHEDRA
      // Alternate
      {{0, 1, 4, 2}}, 
      {{0, 3, 2, 4}}
#else 
      // Regular
      {{0, 1, 4, replaceMe}},
      {{1, 2, 4, replaceMe}},
      {{2, 3, 4, replaceMe}},
      {{3, 0, 4, replaceMe}}
#endif
    }}
  },
  {
    Symmetry::Name::TrigonalBiPyramidal,
    5,
    /* Viewed from the top of the pyramid. The central atom is ( ), 3 and 4 
     * are axial.
     *
     *     3
     *     |  2
     *     | :    <- behind view plane
     * 0--(_)
     *     | \    <- in front of view plane
     *     |  1
     *     4
     */
    &AngleFunctions::trigonalBiPyramidal,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {-0.5, 0.866025, 0},
      {-0.5, -0.866025, 0},
      {0, 0, 1},
      {-0, -0, -1}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {2, 0, 1, 3, 4}, // C3
      {0, 2, 1, 4, 3}, // C2 on 0
      {2, 1, 0, 4, 3}, // C2 on 1
      {1, 0, 2, 4, 3} // C2 on 2
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
      {{0, 1, 3, 2}},
      {{0, 1, 2, 4}}
    }}
  },
  {
    Symmetry::Name::PentagonalPlanar,
    5,
    /* 
     * All in plane:
     *
     *      0
     *  1.  |  .4
     *    °(_)°
     *    /   \
     *   2     3
     *
     */
    &AngleFunctions::pentagonalPlanar,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {0.309017, 0.951057, 0},
      {-0.809017, 0.587785, 0},
      {-0.809017, -0.587785, 0},
      {0.309017, -0.951057, 0}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {4, 0, 1, 2, 3}, // C5
      {0, 4, 3, 2, 1} // C2
    }},
    ConstexprSymmetryInfo::TetrahedraType {{}}
  },
  {
    Symmetry::Name::Octahedral,
    6,
      /* The central atom is ( ), 4 and 5 are axial, the rest equatorial.
       *
       *     4
       *  3  |  2
       *   : | :
       *    (_)        
       *   / | \
       *  0  |  1
       *     5
       *
       * Where /, \ denote bonds in front of the view plane, : denotes bonds
       * behind the view plane.
       */
    &AngleFunctions::octahedral,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {0, 1, 0},
      {-1, -0, -0},
      {-0, -1, -0},
      {0, 0, 1},
      {-0, -0, -1}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {3, 0, 1, 2, 4, 5}, // vertical C4
      {0, 5, 2, 4, 1, 3}, // horizontal C4
      {4, 1, 5, 3, 2, 0} // horizontal C4'
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
#ifdef USE_ALTERNATE_TETRAHEDRA
      // Alternate
      {{3, 0, 4, 5}},
      {{0, 1, 4, 5}},
      {{1, 2, 4, 5}},
      {{2, 3, 4, 5}}
#else
      // Regular
      {{3, 0, 4, replaceMe}},
      {{0, 1, 4, replaceMe}},
      {{1, 2, 4, replaceMe}},
      {{2, 3, 4, replaceMe}},
      {{3, 0, replaceMe, 5}},
      {{0, 1, replaceMe, 5}},
      {{1, 2, replaceMe, 5}},
      {{2, 3, replaceMe, 5}}
#endif
    }}
  },
  {
    Symmetry::Name::TrigonalPrismatic,
    6,
    /* 
     *  3  4  5
     *   : | :
     *    (_)        
     *   : | :
     *  0  1  2
     *
     * Where /, \ denote bonds in front of the view plane, : denotes bonds
     * behind the view plane.
     *
     * Angles 
     *  0-1, 0-2 -> 86°
     *  0-3 -> 76°
     *  0-4, 0-5 -> 134°
     * 
     * From [W(CH3)6], Haaland, Hammel, Rypdal, Volden, J. Am. Chem. Soc. 1990 
     */
    &AngleFunctions::trigonalPrismatic,
    ConstexprSymmetryInfo::CoordinatesType {{
      {0.788011, 0, -0.615661},
      {-0.394005, 0.682437, -0.615661},
      {-0.394005, -0.682437, -0.615661},
      {0.788011, 0, 0.615661},
      {-0.394005, 0.682437, 0.615661},
      {-0.394005, -0.682437, 0.615661}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {2, 0, 1, 5, 3, 4}, // C3 axial
      {5, 4, 3, 2, 1, 0} // C2 betw. 1, 4
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
      // TODO dubious if this captures all relevant information, too limited
      {{replaceMe, 0, 1, 2}},
      {{3, replaceMe, 4, 5}}
    }}
  },
  {
    Symmetry::Name::PentagonalPyramidal,
    6,
    /* 
     *      0
     *  1.  |  .4
     *    °(5)°
     *    :   :
     *   2     3
     *
     * 0-4 are in plane,
     * 5 is above plane,
     * ( ) signifies the central atom beneath it
     */
    &AngleFunctions::pentagonalPyramidal,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {0.309017, 0.951057, 0},
      {-0.809017, 0.587785, 0},
      {-0.809017, -0.587785, 0},
      {0.309017, -0.951057, 0},
      {0, 0, 1}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {4, 0, 1, 2, 3, 5} // C5 axial
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
#ifdef USE_ALTERNATE_TETRAHEDRA
      // Alternate
      {{0, 1, 5, 2}},
      {{2, 3, 5, 4}},
      {{4, 5, replaceMe, 0}}
#else
      // Regular
      {{0, replaceMe, 1, 5}},
      {{1, replaceMe, 2, 5}},
      {{2, replaceMe, 3, 5}},
      {{3, replaceMe, 4, 5}},
      {{4, replaceMe, 0, 5}}
#endif
    }}
  },
  {
    Symmetry::Name::PentagonalBiPyramidal,
    7,
    /* 
     * 3, 5, (_) and 6 in plane, 1 and 2 in front, 0 and 4 in back
     *
     *      5
     *  0_  | .4
     *   _:(_) – 3
     *  1   |°·2
     *      6
     *
     * 0-4 are equatorial, 
     * 5,6 are axial
     */
    &AngleFunctions::pentagonalBiPyramidal,
    ConstexprSymmetryInfo::CoordinatesType {{
      {1, 0, 0},
      {0.309017, 0.951057, 0},
      {-0.809017, 0.587785, 0},
      {-0.809017, -0.587785, 0},
      {0.309017, -0.951057, 0},
      {0, 0, 1},
      {-0, -0, -1}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {4, 0, 1, 2, 3, 5, 6}, // C5 axial
      {1, 0, 4, 3, 2, 6, 5} // C2 equatorial on 3
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
#ifdef USE_ALTERNATE_TETRAHEDRA
      // Alternate
      {{0, 1, 5, 6}},
      {{1, 2, 5, 6}},
      {{2, 3, 5, 6}},
      {{3, 4, 5, 6}},
      {{4, 0, 5, 6}}
#else
      // Regular
      {{0, 1, 5, replaceMe}},
      {{1, 2, 5, replaceMe}},
      {{2, 3, 5, replaceMe}},
      {{3, 4, 5, replaceMe}},
      {{4, 0, 5, replaceMe}},
      {{0, 1, replaceMe, 6}},
      {{1, 2, replaceMe, 6}},
      {{2, 3, replaceMe, 6}},
      {{3, 4, replaceMe, 6}},
      {{4, 0, replaceMe, 6}}
#endif
    }}
  },
  {
    Symmetry::Name::SquareAntiPrismatic,
    8,
    /*  Two representations, one oblique, the other more helpful.
     *  The first is a side-on view. 0 is mostly hidden by 1, 3 mostly hidden 
     *  by 2. 4 and 6 are in the viewing plane while 5 juts out above plane and
     *  7 dips behind plane. 
     *
     *  4   7 5 6
     *    : ·/ :
     *     (__)
     *    ·/  \·
     *   01    23 
     *
     * Below is a top-down view. Strong lines indicate above-plane bonds, dots
     * indicate below-plane bonds.
     *
     *   0   7   3
     *     · | ·
     *   4 –( )– 6
     *     · | ·
     *   1   5   2
     *
     * Angles:
     *
     *   in-plane cis (4, 5) -> 55°
     *   in-plane trans (4, 6) -> 148°
     *   short between planes (0, 4) -> 51°
     *   long between planes (0, 5) -> 175°
     *
     * from [ReF8]2-, Koz'min, P.A., Zhurnal Strukturnoi Khimii 1964 
     * HINT: use the ICSD for angle calculations and preview
     *
     */
    &AngleFunctions::squareAntiprismatic,
    ConstexprSymmetryInfo::CoordinatesType {{
      {-0.00928803, 0.611568, 0.791137},
      {0.795627, 0.605641, -0.0132684},
      {0.795627, -0.605641, -0.0132684},
      {-0.00928803, -0.611568, 0.791137},
      {-0.396172, 0.852169, -0.341841},
      {0.293758, 0, -0.95588},
      {-0.396172, -0.852169, -0.341841},
      {-0.983087, 0, 0.183141}
    }},
    ConstexprSymmetryInfo::RotationsType {{
      {3, 0, 1, 2, 7, 4, 5, 6}, // C4 axial
      /* 180° on equatorial axis in plane with 4, 6 
       *
       *   1   5   2
       *     \ · /
       * – 4 ·( )· 6 –––––– equatorial axis
       *     / · \
       *   0   7   3
       *
       * and 45° anticlockwise on axial axis
       *           
       *   5   2   6
       *     · | ·
       *   1 –( )– 3
       *     · | ·
       *   4   0   7
       *
       */
      {5, 4, 7, 6, 1, 0, 3, 2}, 
    }},
    ConstexprSymmetryInfo::TetrahedraType {{
#ifdef USE_ALTERNATE_TETRAHEDRA
      // Alternate
      {{0, 1, 4, 6}},
      {{1, 2, 5, 7}},
      {{2, 3, 6, 4}},
      {{3, 0, 7, 5}},
#else 
      // Regular
      {{7, 0, 4, replaceMe}},
      {{0, 4, replaceMe, 1}},
      {{4, 1, 5, replaceMe}},
      {{1, 5, replaceMe, 2}},
      {{5, 2, 6, replaceMe}},
      {{2, 6, replaceMe, 3}},
      {{6, 3, 7, replaceMe}},
      {{3, 7, replaceMe, 0}}
#endif
    }}
  }
}};

} // namespace Symmetry
