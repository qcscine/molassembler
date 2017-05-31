#include "Symmetries.h"

#include <cassert>
#include <algorithm>
#include <Eigen/Geometry>

#include "constexpr_magic/Math.h"
#include "template_magic/TemplateMagic.h"

namespace Symmetry {

const std::vector<Name> allNames {
  Name::Linear, // 2
  Name::Bent,
  Name::TrigonalPlanar, // 3
  Name::TrigonalPyramidal,
  Name::TShaped,
  Name::Tetrahedral, // 4
  Name::SquarePlanar,
  Name::Seesaw,
  Name::SquarePyramidal, // 5
  Name::TrigonalBiPyramidal,
  Name::PentagonalPlanar,
  Name::Octahedral, // 6
  Name::TrigonalPrismatic,
  Name::PentagonalPyramidal,
  Name::PentagonalBiPyramidal, // 7
  Name::SquareAntiPrismatic // 8
};

const std::map<Name, SymmetryInformation> symmetryData {
  {
    Name::Linear, 
    SymmetryInformation {
      "linear",
      2,
      RotationsList {
        {1, 0}
      },
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        return M_PI;
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double { 
        /* subject to a lot of variation, between 90 and 109 degrees pursuant to 
         * english wikipedia, using experimental data here to improve instances
         * of this geometry on e.g. O center would a big improvement to DG runs
         */
        if(a == b) {
          return 0;
        }

        return ConstexprMagic::Math::toRadians(107); 
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        return ConstexprMagic::Math::toRadians(120);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        return ConstexprMagic::Math::toRadians(107.5);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        if((a + b) % 2 == 1) {
          return M_PI / 2;
        } 

        return M_PI;
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        return ConstexprMagic::Math::toRadians(109.5);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        if((a + b) % 2 == 1) {
          // this expression indicates cis
          return M_PI / 2;
        } 

        // leftover case is trans
        return M_PI;
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        const auto& smaller = std::min(a, b);
        const auto& larger = std::max(a, b);
        if(smaller == 0 && larger == 3) {
          return M_PI;
        }

        if(smaller == 1 && larger == 2) {
          return ConstexprMagic::Math::toRadians(120);
        }

        return M_PI / 2;
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        unsigned smaller = std::min(a, b), larger = std::max(a, b);
        if(larger < 3) {
          // -> smaller < 2, this means either 0,1 0,2 1,2 axial
          return ConstexprMagic::Math::toRadians(120);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        unsigned absDiff = std::min(a - b, b - a);
        return std::min(
          absDiff,
          std::min(absDiff - 5, 5 - absDiff)
        ) * ConstexprMagic::Math::toRadians(72);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        } 
        
        // Between plane symmetric
        if(std::min(a - b, b - a) == 3) {
          return ConstexprMagic::Math::toRadians(76);
        } 

        // In plane triangle
        if(
          (a < 3 && b < 3)
          || (a >= 3 && b >= 3)
        ) {
          return ConstexprMagic::Math::toRadians(86);
        } 

        // Between plane asymmetric
        return ConstexprMagic::Math::toRadians(134);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
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
        ) * ConstexprMagic::Math::toRadians(72);
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
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        if(a + b == 11) {
          return M_PI; // trans 5,6
        }

        if(TemplateMagic::XOR(a > 4, b > 4)) {
          return M_PI / 2; // any angle to axial index
        }

        // remainder are equatorial angles, like PentagonalPlanar
        unsigned absDiff = std::min(a - b, b - a);
        return std::min(
          absDiff,
          std::min(absDiff - 5, 5 - absDiff)
        ) * ConstexprMagic::Math::toRadians(72);
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
#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        }

        return squareAntiprismaticAngles.at(
          std::min(a, b),
          std::max(a, b)
        );
      },
#else
      [](
        const unsigned& a,
        const unsigned& b
      ) -> double {
        if(a == b) {
          return 0;
        } 
        
        if(
          (a < 4 && b < 4)
          || (a >= 4 && b >= 4)
        ) { // in plane
          if((a + b) % 2 == 1) { // cis
            return ConstexprMagic::Math::toRadians(72.9875); 
          } 
          
          // otherwise trans
          return ConstexprMagic::Math::toRadians(114.475);
        }
        
        // remaining cases are between planes
        unsigned minDiff = std::min(a - b, b - a);
        if(minDiff == 3 || minDiff == 4 || minDiff == 7) { // short
          return ConstexprMagic::Math::toRadians(78.05);
        } 

        // lsat case is long between planes
        return ConstexprMagic::Math::toRadians(142.275);
      },
#endif
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

} // namespace Symmetry
