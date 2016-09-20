#include <iostream>
#include "AtomSet.h"

namespace MoleculeManip {
    struct Molecule {
    /* Data members */
        /* A Molecule conceptually contains a graph:
         * - Atoms are vertices (and thus have values)
         * - Bonds are edges (and thus weighted)
         * - The ensuing graph is
         *   - connected: a path from any node to any other exists
         *   - sparse: few edges present compared to the number of possible 
         *     edges -> use adjacency list instead of adjacency matrix
         */
        Delib::AtomSet atoms;
        // + Custom implementation of an AdjacencyList
        // + Custom implementation of a StereocenterList

    /* Operators */
        /* An efficient implementation of the following two is imperative.
         */
        bool operator == (const Molecule& b) const;
        bool operator != (const Molecule& b) const;

        /* Output stream operator for easier debugging. */
        friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);
    };
}
