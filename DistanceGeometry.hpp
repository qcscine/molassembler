#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

#include <vector>
#include <Eigen/Core>

#include "AtomSet.h"
#include "Types/PositionCollection.h"
#include "Molecule.hpp"

namespace DistanceGeometry {

    template<typename T>
    using ContainerType = std::vector<T>;
    
    /* Enum types */
    enum MetrizationOption {
        off,
        partial,
        full
    };

    enum RefinementOption {
        threeDimensional,
        fourDimensional
    };

    /* Interfacing functions */
    Delib::AtomSet embed(
        const MoleculeManip::Molecule& mol,
        const MetrizationOption metrizationOption,
        const RefinementOption refinementOption
    );


    /* Internal functions */
    // Embedding implementation
    Delib::AtomSet  _embed(
        const MoleculeManip::Molecule& mol
    );

    Eigen::MatrixXd _generateDistanceBounds(
        const MoleculeManip::Molecule& mol
    );

    Eigen::MatrixXd _triangleInequalitySmooth(
        const Eigen::MatrixXd& matrix
    );

    Eigen::MatrixXd _generateDistancesMatrix(
        const Eigen::MatrixXd& distanceBoundsMatrix,
        const enum MetrizationOption = MetrizationOption::off
    );

    ContainerType<Eigen::VectorXd> _refine(
        const MoleculeManip::Molecule& mol,
        const ContainerType<Eigen::VectorXd>& coordinates,
        const enum RefinementOption = RefinementOption::fourDimensional
    );

}

// Implementation
/* Internal functions */

/* distance geometry embedding steps: 
 *
 * - generate a distance bounds matrix using empirical information
 * - smooth it using triangle bounds
 * - generate a random distances matrix, ideally with partial metrization
 *   (meaning you store an additional 4N coordinates, after each random
 *   distance generation update the bounds using the inequalities)
 * - convert the distances matrix to a metric matrix
 * - calculate the eigenvalues and eigenvectors
 * - project the top four eigenvalues into four dimensional space
 * - conjugate gradient minimization of 4D space, modification of the error
 *   function when all chiral centers have correct stereochemistry
 */

#endif
