#ifndef INCLUDE_MOLASSEMBLER_CONFORMER_GENERATION_H
#define INCLUDE_MOLASSEMBLER_CONFORMER_GENERATION_H

#include "DistanceGeometry/generateConformation.h"

namespace molassembler {

/*! Generate a conformational ensemble of a Molecule
 *
 * Returns a result type which may or may not contain a vector of
 * PositionCollections. The result type is much like an optional, except that
 * in the error case it carries data about the error in order to help diagnose
 * possible mistakes made in the molecular graph specification.
 */
outcome::result<
  std::vector<Delib::PositionCollection>
> generateEnsemble(
  const Molecule& molecule,
  const unsigned numStructures
);

/*! Generate a 3D structure of a Molecule
 *
 * Returns a result type which may or may not contain a PositionCollection. The
 * result type is much like an optional, except that in the error case it
 * carries data about the error in order to help diagnose possible mistakes
 * made in the molecular graph specification.
 */
outcome::result<Delib::PositionCollection> generateConformation(const Molecule& molecule);

} // namespace molassembler

#endif
