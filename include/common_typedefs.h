#ifndef INCLUDE_COMMON_TYPEDEFS_H
#define INCLUDE_COMMON_TYPEDEFS_H

#include <tuple>

namespace MoleculeManip {

    /* Common typedefs */
    enum BondType : uint8_t {
        Single,
        Double,
        Triple,
        Quadruple,
        Quintuple,
        Sextuple,
        Aromatic
    };

    // TODO how to handle subtypes?
    enum StereocenterType : uint8_t {
        Tetrahedral,
        Octahedral,
        CisTrans
    };

    /*
    using unsigned_type = unsigned;
    using EdgeIndexType = uint32_t;
    using AtomIndexType = uint16_t; // 65k max is sufficient
    */
    // TODO temp for testing
    using unsigned_type = long long unsigned;
    using EdgeIndexType = uint64_t;
    using AtomIndexType = uint64_t; // 65k max is sufficient


    /* derived types */
    using EdgeType = std::tuple<
        AtomIndexType,
        AtomIndexType,
        BondType
    >;

    /*!
     * Contains static functions for common operations with EdgeTypes
     */
    struct EdgeTypeHelpers {
        /*!
         * Lexicographical less-than operator for a tuple and two atom indices
         * \param tuple The tuple 
         * \param a The first atom index
         * \param a The second atom index
         * \returns Whether the tuple is smaller than the tuple defined by the
         *  two atom indices.
         */
        static bool _tupleLexicographicallySmallerIdxs(
            const EdgeType& tuple,
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept {
            return (
                std::get<0>(tuple) < a
                || (
                    std::get<0>(tuple) == a
                    && std::get<1>(tuple) < b 
                )
            );
        }

        /*!
         * Lexicographical greater than operator for a tuple and two atom
         * indices
         * \param tuple The tuple 
         * \param a The first atom index
         * \param a The second atom index
         * \returns Whether the tuple is greater than the tuple defined by the
         *  two atom indices.
         */
        static bool _tupleLexicographicallyGreaterIdxs(
            const EdgeType& tuple,
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept {
            return (
                std::get<0>(tuple) > a 
                || (
                    std::get<0>(tuple) == a 
                    && std::get<1>(tuple) > b 
                )
            );
        }

        /*!
         * Lexicographical less-than operator for two EdgeTypes
         * \param a The first tuple 
         * \param b The second tuple
         * \returns Whether the first tuple is smaller than the second
         *  lexicographically
         */
        static bool tupleSmaller(
            const EdgeType& a,
            const EdgeType& b
        ) noexcept {
            return (
                std::get<0>(a) < std::get<0>(b)
                || (
                    std::get<0>(a) == std::get<0>(b)
                    && std::get<1>(a) < std::get<1>(b)
                )
            );
        }
    };
}

#endif
