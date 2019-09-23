/* @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher PointGroup.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "chemical_symmetries/CoordinateSystemTransformation.h"
#include "chemical_symmetries/Names.h"
#include "chemical_symmetries/Symmetries.h"
#include "chemical_symmetries/Recognition.h"

#include "temple/Functional.h"
#include "temple/Stringify.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Scine;
using namespace Symmetry;

template<typename EnumType>
constexpr inline auto underlying(const EnumType e) {
  return static_cast<std::underlying_type_t<EnumType>>(e);
}

PositionCollection addOrigin(const PositionCollection& vs) {
  const unsigned N = vs.cols();
  PositionCollection positions(3, N + 1);
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) = vs.col(i);
  }

  // Add origin point explicitly to consideration
  positions.col(N) = Eigen::Vector3d::Zero(3);
  return positions;
}

void distort(Eigen::Ref<PositionCollection> positions, const double distortionNorm = 0.01) {
  const unsigned N = positions.cols();
  for(unsigned i = 0; i < N; ++i) {
    positions.col(i) += distortionNorm * Eigen::Vector3d::Random().normalized();
  }
}

const std::vector<std::string>& pointGroupStrings() {
  static const std::vector<std::string> strings {
    "C1", "Ci", "Cs",
    "C2", "C3", "C4", "C5", "C6", "C7", "C8",
    "C2h", "C3h", "C4h", "C5h", "C6h", "C7h", "C8h",
    "C2v", "C3v", "C4v", "C5v", "C6v", "C7v", "C8v",
    "S4", "S6", "S8",
    "D2", "D3", "D4", "D5", "D6", "D7", "D8",
    "D2h", "D3h", "D4h", "D5h", "D6h", "D7h", "D8h",
    "D2d", "D3d", "D4d", "D5d", "D6d", "D7d", "D8d",
    "T", "Td", "Th",
    "O", "Oh",
    "I", "Ih",
    "Cinfv", "Dinfh"
  };

  return strings;
}

const std::string& pointGroupString(PointGroup group) {
  return pointGroupStrings().at(underlying(group));
}

const std::vector<std::string>& topNames() {
  static const std::vector<std::string> strings {
    "Linear",
    "Asymmetric",
    "Prolate",
    "Oblate",
    "Spherical"
  };

  return strings;
}

const std::string& topName(Top top) {
  return topNames().at(underlying(top));
}

BOOST_AUTO_TEST_CASE(Recognition) {
  const std::map<Name, PointGroup> expected {
    {Name::Bent, PointGroup::C2v}
  };

  // const std::map<Name, PointGroup> expected {
  //   {Name::Linear, PointGroup::Dinfh},
  //   {Name::Bent, PointGroup::C2v},
  //   {Name::TrigonalPlanar, PointGroup::D3h},
  //   {Name::CutTetrahedral, PointGroup::C3v},
  //   {Name::TShaped, PointGroup::C2v},
  //   {Name::Tetrahedral, PointGroup::Td},
  //   {Name::SquarePlanar, PointGroup::D4h},
  //   {Name::Seesaw, PointGroup::C2v},
  //   {Name::TrigonalPyramidal, PointGroup::C3v},
  //   {Name::SquarePyramidal, PointGroup::C4v},
  //   {Name::TrigonalBiPyramidal, PointGroup::D3h},
  //   {Name::PentagonalPlanar, PointGroup::D5h},
  //   {Name::Octahedral, PointGroup::Oh},
  //   {Name::TrigonalPrismatic, PointGroup::D3h},
  //   {Name::PentagonalPyramidal, PointGroup::C5v},
  //   {Name::PentagonalBiPyramidal, PointGroup::D5h},
  //   {Name::SquareAntiPrismatic, PointGroup::D4d}
  // };

  for(const auto& nameGroupPair : expected) {
    auto positions = addOrigin(symmetryData().at(nameGroupPair.first).coordinates);
    // distort(positions);
    auto normalized = detail::normalize(positions);

    // Add a random coordinate transformation
    const Eigen::Matrix3d rot = rotationMatrix(CoordinateSystem {}, CoordinateSystem::random());
    for(unsigned i = 0; i < normalized.cols(); ++i) {
      normalized.col(i) = rot * normalized.col(i);
    }

    // Standardize the top
    Top top = standardizeTop(normalized);

    const double pgCSM = csm::pointGroup(
      normalized,
      nameGroupPair.second
    ).value_or(1000);
    BOOST_REQUIRE_MESSAGE(
      pgCSM != 1000,
      "Could not calculate "
      << pointGroupString(nameGroupPair.second)
      << " CSM for " << Symmetry::name(nameGroupPair.first)
    );

    BOOST_CHECK_MESSAGE(
      pgCSM < 1,
      "Expected CSM(" << pointGroupString(nameGroupPair.second)
      << ") < 1 for " << Symmetry::name(nameGroupPair.first)
      << ", got " << pgCSM << " (top is " << topName(top) << ")"
    );
  }
}

// TODO reframe as a test, not just as output testing
BOOST_AUTO_TEST_CASE(PointGroupElements) {
  auto writeXYZ = [](const std::string& filename, const PositionCollection& positions) {
    std::ofstream outfile(filename);
    const unsigned N = positions.cols();
    outfile << N << "\n\n";
    outfile << std::fixed << std::setprecision(10);
    for(unsigned i = 0; i < N; ++i) {
      outfile << std::left << std::setw(3) << "H";
      outfile << std::right
        << std::setw(16) << positions.col(i).x()
        << std::setw(16) << positions.col(i).y()
        << std::setw(16) << positions.col(i).z()
        << "\n";
    }
    outfile.close();
  };

  /* For each point group, create a point at unit x and z and apply all
   * transformations to each point
   */
  auto writePointGroup = [&](const PointGroup group) {
    const auto elements = elements::symmetryElements(group);
    const auto groupings = elements::npGroupings(elements);
    std::cout << pointGroupStrings().at(underlying(group)) << ": {";
    for(const auto& element : elements) {
      std::cout << element -> name() << ", ";
    }
    std::cout << "}\n";
    for(const auto& iterPair : groupings) {
      std::cout << "  np = " << iterPair.first << " along " << iterPair.second.probePoint.transpose() << " -> " << temple::stringify(iterPair.second.groups) << "\n";
      std::cout << "  ";
      std::cout << temple::stringifyContainer(iterPair.second.groups,
        [&](const auto& grp) -> std::string {
          return temple::stringifyContainer(grp,
            [&elements](const unsigned elementIdx) -> std::string {
              return elements.at(elementIdx)->name();
            }
          );
        }
      ) << "\n";
    }
  };

  // const PointGroup limit = PointGroup::Td;
  // for(unsigned g = 0; g < underlying(limit); ++g) {
  //   writePointGroup(static_cast<PointGroup>(g));
  // }
  writePointGroup(PointGroup::C2v);
}

BOOST_AUTO_TEST_CASE(PaperCSMExamples) {
  /* From Continuous Symmetry Measures. 2. Symmetry Groups and the Tetrahedron,
   * Zabrodsky, Peleg, Avnir. J. Am. Chem. Soc. 1993, 115, 8278-8289
   *
   * Example on p. 8283 from their ref 12
   */
  PositionCollection tetrahedralPositions(3, 4);
  tetrahedralPositions.col(0) = Eigen::Vector3d {0.0, 0.0, 1.645};
  tetrahedralPositions.col(1) = Eigen::Vector3d {0.0, 1.518860, -0.347028};
  tetrahedralPositions.col(2) = Eigen::Vector3d {-1.286385, -0.700083, -0.391603};
  tetrahedralPositions.col(3) = Eigen::Vector3d {1.179085, -0.755461, -0.372341};

  const double expectedTdCSM = 0.17;
  const double calculatedTd = csm::pointGroup(
    detail::normalize(tetrahedralPositions),
    PointGroup::Td
  ).value_or(1000);

  std::cout << "Calculated Td CSM: " << calculatedTd << " (expect " << expectedTdCSM << ")\n";
  BOOST_CHECK_CLOSE(expectedTdCSM, calculatedTd, 1);

  /* Example from page 8284 */
  PositionCollection c3vPositions(3, 4);
  c3vPositions.col(0) = Eigen::Vector3d {0.0, 0.0, 1.645};
  c3vPositions.col(1) = Eigen::Vector3d {0.0, 0.87461971, -0.48460962};
  c3vPositions.col(2) = Eigen::Vector3d {0.75128605, -0.39338892, -0.52991926};

  /* NOTE: Sign inversion on x is speculative, paper reprints P3 = P4 positions
   * although this does not make sense to me.
   */
  c3vPositions.col(3) = Eigen::Vector3d {-0.75128605, -0.39338892, -0.52991926};

  const double expectedC3vCSM = 1.16;
  const double calculatedC3v = csm::pointGroup(
    detail::normalize(c3vPositions),
    PointGroup::C3v
  ).value_or(1000);

  std::cout << "Calculated C3v CSM: " << calculatedC3v << " (expect " << expectedC3vCSM << ")\n";
  BOOST_CHECK_CLOSE(expectedC3vCSM, calculatedC3v, 1);
}

BOOST_AUTO_TEST_CASE(SquarePlanarC4D4PointGroups) {
  const double pointGroupCSM = csm::pointGroup(
    detail::normalize(symmetryData().at(Name::SquarePlanar).coordinates),
    PointGroup::C4
  ).value_or(1000);
  BOOST_CHECK_MESSAGE(
    std::fabs(pointGroupCSM) < 1e-10,
    "C4 point group CSM on square planar coordinates is not zero, but " << pointGroupCSM
  );

  const double D4CSM = csm::pointGroup(
    detail::normalize(symmetryData().at(Name::SquarePlanar).coordinates),
    PointGroup::D4
  ).value_or(1000);

  BOOST_CHECK_MESSAGE(
    0 < D4CSM && D4CSM < 1e-10,
    "D4 CSM on square planar coordinates is not zero, but " << D4CSM
  );
}

BOOST_AUTO_TEST_CASE(InertialStandardization) {
  const std::vector<
    std::pair<Name, Top>
  > tops {
    {Name::Linear, Top::Linear},
    {Name::Bent, Top::Asymmetric},
    {Name::Seesaw, Top::Asymmetric},
    {Name::TrigonalBiPyramidal, Top::Prolate},
    {Name::PentagonalBiPyramidal, Top::Oblate},
    {Name::SquarePlanar, Top::Oblate},
    {Name::Octahedral, Top::Spherical},
    {Name::Tetrahedral, Top::Spherical}
  };

  for(const auto nameTopPair : tops) {
    auto positions = addOrigin(symmetryData().at(nameTopPair.first).coordinates);

    // Apply a random coordinate transformation
    const Eigen::Matrix3d R = rotationMatrix(CoordinateSystem {}, CoordinateSystem::random());
    for(unsigned i = 0; i < positions.cols(); ++i) {
      positions.col(i) = R * positions.col(i);
    }

    // Analyze it
    auto normalizedPositions = detail::normalize(positions);
    Top top = standardizeTop(normalizedPositions);
    BOOST_CHECK_MESSAGE(
      top == nameTopPair.second,
      "Top standardization failed. Expected "
      << topNames().at(static_cast<unsigned>(nameTopPair.second))
      << " for symmetry " << Symmetry::name(nameTopPair.first)
      << ", got " << topNames().at(static_cast<unsigned>(top)) << " instead."
    );
  }
}

BOOST_AUTO_TEST_CASE(FixedCnAxis) {
  const std::vector<std::pair<Name, unsigned>> highestOrderAxis {
    {Name::Bent, 2},
    {Name::TrigonalPlanar, 3},
    {Name::CutTetrahedral, 3},
    {Name::TShaped, 2},
    {Name::Tetrahedral, 3},
    {Name::SquarePlanar, 4},
    {Name::Seesaw, 2},
    {Name::TrigonalPyramidal, 3},
    {Name::SquarePyramidal, 4},
    {Name::TrigonalBiPyramidal, 3},
    {Name::PentagonalPlanar, 5},
    {Name::Octahedral, 4},
    {Name::TrigonalPrismatic, 3},
    {Name::PentagonalPyramidal, 5},
    {Name::PentagonalBiPyramidal, 5},
    {Name::SquareAntiPrismatic, 4}
  };

  const Eigen::Matrix3d axes = Eigen::Matrix3d::Identity();

  constexpr double acceptanceThreshold = 0.1;

  for(const auto& nameOrderPair : highestOrderAxis) {
    std::cout << "For symmetry " << Symmetry::name(nameOrderPair.first) << "\n";
    auto positions = addOrigin(symmetryData().at(nameOrderPair.first).coordinates);
    distort(positions);

    auto normalizedPositions = detail::normalize(positions);
    standardizeTop(normalizedPositions);

    boost::optional<unsigned> highestFoundOrder;
    for(unsigned i = 0; i < 3; ++i) {
      Eigen::Vector3d axis = axes.col(i);
      char axisChar = "xyz"[i];
      for(unsigned n = 2; n < 6; ++n) {
        const double cn = csm::element(normalizedPositions, elements::Rotation::Cn(axis, n));

        if(cn < acceptanceThreshold) {
          std::cout << "Found ";
          highestFoundOrder = n;
        }

        std::cout << "C" << n << " = " << cn << " along " << axisChar << "-axis\n";
      }
    }

    BOOST_CHECK_MESSAGE(
      highestFoundOrder,
      "No Cn axis found along any principal moment axis for "
      << Symmetry::name(nameOrderPair.first)
    );
    if(highestFoundOrder) {
      BOOST_CHECK_MESSAGE(
        highestFoundOrder.value() == nameOrderPair.second,
        "Expected to find Cn of order " << nameOrderPair.second << ", but found " << highestFoundOrder.value() << " instead."
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(AsymmetricTopStandardization) {
  auto coordinates = addOrigin(symmetryData().at(Symmetry::Name::Bent).coordinates);
  auto normalizedPositions = detail::normalize(coordinates);
  Top top = standardizeTop(normalizedPositions);
  BOOST_REQUIRE(top == Top::Asymmetric);

  unsigned highestAxisOrder = reorientAsymmetricTop(normalizedPositions);
  BOOST_CHECK(highestAxisOrder == 2);

  // Ensure rotation of highest order axis to z worked
  BOOST_CHECK(csm::element(normalizedPositions, elements::Rotation::Cn(Eigen::Vector3d::UnitZ(), 2)) < 1e-10);
}
