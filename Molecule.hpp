#include <iostream>
#include "AtomSet.h"
#include "ConnectivityManager/ConnectivityManager.hpp"

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

        // The set of QC data on the atoms
        Delib::AtomSet atoms;
        // The information on interconnectedness of the atoms
        FastConnectivityManager connectivity;
        // + Custom implementation of a StereocenterList

    /* Operators */
        /* An efficient implementation of the following two is imperative.
         * Some ideas for fast differentiation can probably be found from the 
         * wikipedia category "Graph invariants"
         */
        bool operator == (const Molecule& b) const;
        bool operator != (const Molecule& b) const;

        /* Output stream operator for easier debugging. */
        friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);
    };
}
