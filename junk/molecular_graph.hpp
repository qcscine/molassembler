#ifndef LIB_INCLUDE_MOLECULAR_GRAPH_HPP
#define LIB_INCLUDE_MOLECULAR_GRAPH_HPP

#include <vector>
#include "Atom.h"


namespace molecular_graph {

    template<typename T>
    using Container = std::vector<T>;

    template<typename T>
    using vec2d = std::vector< std::vector<T> >;

    enum BondType { 
        Single,
        Double,
        Triple,
        Quadruple,
        Quintuple,
        Sextuple 
    };

    class BondAdjacencyList {
    public:
        BondAdjacencyList();
        uint8_t connected_components();

    private:
        std::vector< std::vector< BondType > > _adjacencies;

        uint8_t _connected_components(std::vector<uint8_t>& component_count);

    };


    class MolecularGraph {
    public:
        /* Constructors
         * - Having one from const Delib::AtomSet& would be very beneficial
         */
        MolecularGraph();

        // A Simple property calls
        // - Get the number of atoms
        auto num_atoms();

        // - Get the number of connected groups
        auto num_groups();

    private:
        Container<Delib::Atom> _atoms;
        BondAdjacencyList _adjacencies;
    };

}

#endif
