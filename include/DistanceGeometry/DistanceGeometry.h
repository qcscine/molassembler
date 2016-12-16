#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

namespace MoleculeManip {

namespace DistanceGeometry {

/* Enum types */
enum class MetrizationOption {
  off,
  partial,
  full
};

enum class EmbeddingOption : unsigned {
  threeDimensional = 3,
  fourDimensional = 4
};

} // eo namespace DistanceGeometry

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

} // eo namespace MoleculeManip

#endif
