#include <vector>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <cassert>

#include "common_typedefs.hpp"
#include "AdjacencyList.hpp"
#include "EdgeList.hpp"

/* TODO
 * - Nothing ensures that the molecular graph stays connected. 
 *   Either store chords or detect disconnects in remove_* functions
 *
 *   -> If the state of the ConnectivityManager is to be connected at ALL
 *   times, then the API of add_atom must change and remove_atom and
 *   remove_bond need added checks to ensure there is only one component in the
 *   graph. If add_bond is to exist only as::
 *
 *   AtomIndexType add_atom(
 *       const AtomIndexType& bonded_to,
 *       const BondType& bond_type
 *   )
 *
 *   then the ConnectivityManager must initialize with at least one atom. 
 *   This however awkardly forces dependencies of the class to initialize
 *   similarly :(
 *
 * - How to represent aromatic structures internally? With alternating bonds or
 *   with an additional bond type? Clear preference for an additional bond type
 *   for easier parsing and deduction of min-max distances for distance
 *   geometry
 *
 * - Graph algorithms -> ring detection, component detection
 */

using namespace std;

template<typename T, typename U>
ostream& operator << (
    ostream& os,
    const std::pair<T, U>& pair
) {
    os << "(" << pair.first << ", " << pair.second << ")";
    return os;
}

namespace MoleculeManip {

    class ConnectivityManager {
    private:
    /* Private members */
        EdgeList _edges;
        AdjacencyList _adjacencies;
        std::vector<bool> _atom_exists;

        /*!
         * Checks whether the supplied atom index is valid, i.e. in range.
         * \param[in] a The atom index to check
         * \returns whether the supplied atom index is valid, i.e. in range.
         */
        inline bool _valid_atom_index(const AtomIndexType& a) const noexcept {
            return (
                a < _adjacencies.size()
                && _atom_exists[a]
            );
        }

        /*!
         * Checks whether the supplied atom indices are valid, i.e. in range 
         * and not identical.
         * \param[in] a The first atom index to check
         * \param[in] b The second atom index to check
         * \returns whether the supplied atom indices are valid.
         */
        inline bool _valid_atom_indices(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept {
            return (
                _valid_atom_index(a)
                && _valid_atom_index(b)
                && a != b
            );
        }

    public:
        /*!
         * Add an atom to the ConnectivityManager.
         * \returns The index of the new atom.
         * \note is O(1)
         */
        AtomIndexType add_atom() noexcept {
            // add an empty vector to _adjacencies at the end.
            _adjacencies.add_slot();
            // add an entry to _atom_exists
            _atom_exists.emplace_back(true);

            return _atom_exists.size() - 1;
        };

        /*!
         * Adds a bond to the ConnectivityManager. 
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to register
         * \pre a, b must fulfill a != b
         * \note is O( log(E) ), where
         * - E is the number of stored edges
         */
        void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept {
            assert(_valid_atom_indices(a, b));

            _edges.add(
                make_tuple(
                    a,
                    b,
                    bond_type
                )
            );

            _adjacencies.add_adjacency(a, b);
        }

        /*!
         * Checks whether two atoms are bonded.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns Whether two atoms are bonded or not.
         * \pre a, b must fulfill a != b
         * \note is O( B ), where
         * - B is the number of stored bonds for atom a
         */
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept {
            assert(_valid_atom_indices(a, b));

            /* could binary search the edge list immediately, but linear search
             * in the adjacency list is probably faster if we don't need to 
             * know the bond type
             */
            return _adjacencies.is_adjacent(a, b);
        }

        /*!
         * Checks whether two atoms are bonded by the specified BondType.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to check
         * \returns Whether the two atoms are bonded as specified.
         * \pre a, b must fulfill: a != b.
         * \note is O( log(E) ), where
         * - E is the number of stored edges
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept {
            assert(_valid_atom_indices(a, b));

            /* for this case, checking the adjacency list is not good enough
             * the bond information is only in the edge list, binary search it
             */
            auto found_and_pos_pair = _edges.search(
                a,
                b
            );
            return (
                found_and_pos_pair.first
                && (
                    bond_type == std::get<2>(
                        _edges.get(found_and_pos_pair.second)
                    )
                )
            );
        }

        /*!
         * Gets the BondType between two atoms.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns The bondtype between the atoms.
         * \pre a, b must fulfill: a != b.
         * \note is O( log(E) ), where
         * - E is the number of stored edges
         */
        BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept {
            assert(_valid_atom_indices(a, b));

            auto found_and_pos_pair = _edges.search(
                a,
                b
            );
            /* found_and_pos_pair.second is the index where the found tuple is
             *  stored
             * _edges[found_and_pos_pair.second] is the found tuple
             * std::get<2>(_edges[found_and_pos_pair.second]) is the 3rd
             *  element in the tuple -> the bondtype
             */
            return std::get<2>(
                _edges.get(found_and_pos_pair.second)
            );
        }

        /*!
         * Get a vector of all bonds outwards from an atom.
         * \param[in] a The atom index
         * \returns A vector of all bonds outwards from an atom as std::pairs,
         *  first is the BondType, second the atom index.
         * \note is O( B log (E) ), where
         * - B is the number of bonds from atom a
         * - E is the number of stored edges
         * -> O( log (E) ) with small constants
         */
        vector<
            pair<
                AtomIndexType,
                BondType
            >
        > get_bond_pairs(const AtomIndexType& a) const noexcept {
            assert(_valid_atom_index(a));

            /* We can get a full list of connections and how they are
             * ordered in the edge list from the adjacencies.
             */
            auto bonded_to = _adjacencies.get_adjacencies(a);
            vector<
                pair<
                    AtomIndexType,
                    BondType
                >
            > bond_pairs;

            /* Get the bond type for every bonded atom index by binary 
             * searching the edges
             */
            for(const auto& bonded_atom_index : bonded_to) {
                auto found_and_index_pair = _edges.search(
                    std::min(a, bonded_atom_index),
                    std::max(a, bonded_atom_index)
                );

                if(found_and_index_pair.first) {
                    bond_pairs.emplace_back(
                        bonded_atom_index,
                        std::get<2>(
                            _edges.get(found_and_index_pair.second)
                        )
                    );
                }
            }

            return bond_pairs;
        }

        /*!
         * Get a vector of all bonded atom indices.
         * \param[in] a The atom index
         * \returns A vector all bonded atom indices.
         * \note is O(1)
         */
        vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const noexcept {
            assert(_valid_atom_index(a));

            return _adjacencies.get_adjacencies(a);
        }

        /*!
         * Removes an atom from the ConnectivityManager. 
         * \param[in] a The atom index
         * \note is O( B log (E) + sum_all(e_i) ), where 
         * - B is the number of bonds the atom being removed has
         * - E is the number of stored edges
         * - e_i are the number of edges of the atoms bonded to the atom being
         *   removed
         * -> O( log (E) ) with some constants involved
         */
        void remove_atom(
            const AtomIndexType& a
        ) noexcept {
            assert(_valid_atom_index(a));

            // must update edges and adjacencies 
            auto bonded_to = _adjacencies.get_adjacencies(a);

            // erase all edges to and from this atom
            for(const auto& bonded_atom_index : bonded_to) {
                _edges.remove(
                    std::min(a, bonded_atom_index),
                    std::max(a, bonded_atom_index)
                );
            }

            // remove all other mentions in _adjacencies
            for(const auto& bonded_atom_index : bonded_to) {
                _adjacencies.remove_adjacency(a, bonded_atom_index);
            }

            // set _atom_exists to false for this index
            _atom_exists[a] = false;
        }

        /*!
         * Removes a bond from the ConnectivityManager. If the bond doesn't exist,
         * the function does not throw.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \pre a, b must fulfill: a != b.
         * \note is O( e_a + e_b + log (E) ), where
         * - e_i is the number of edges involving atom i
         * - E is the number of stored edges
         * -> O( log(E) ) with negligible constants
         */
        void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept {
            assert(_valid_atom_indices(a, b));

            // must update _edges and _adjacencies
            _adjacencies.remove_adjacency(a, b);

            _edges.remove(
                a,
                b
            );
        };

        void reset() {
            _edges = EdgeList();
            _adjacencies = AdjacencyList();
            _atom_exists = std::vector<bool>();
        }

        /*!
         * A debug helper function. Checks that the internal state is valid.
         * \note is O(E)
         */
        std::pair<bool, std::string> validate() const noexcept {
            // basic checks
            if(_adjacencies.size() != _atom_exists.size()) {
                return make_pair(
                    false,
                    "The size of _adjacencies and _atom_exists is mismatched."
                );
            }

            // is the edge list ordered
            if(!_edges.is_ordered()) {
                return make_pair(
                    false,
                    "The list of edges is unordered."
                );
            }

            /* all bonds specified in _adjacencies are present in _edges and
             *  all specified atom indices are valid
             */
            for(AtomIndexType from = 0; from < _adjacencies.size(); from++) {
                for(const auto& to : _adjacencies.get_adjacencies(from)) {
                    if(!_valid_atom_indices(from, to)) {
                        return std::make_pair(
                            false,
                            "One or both atom indices ("s
                                + std::to_string(from)
                                + ", "
                                + std::to_string(to)
                                +") are invalid."
                        );
                    }
                    auto found_and_pos_pair = _edges.search(
                        std::min(from, to),
                        std::max(from, to)
                    );
                    if(!found_and_pos_pair.first) {
                        return std::make_pair(
                            false,
                            "The edge ("s 
                                + std::to_string(from)
                                + ", "
                                + std::to_string(to)
                                +") in _adjacencies is not in _edges."
                        );
                    }
                }
            }

            /* all edges in _edges are also in _adjacencies and involve only 
             * valid atom indices
             */
            for(const auto& edge : _edges) {
                if(!_valid_atom_indices(
                    std::get<0>(edge), 
                    std::get<1>(edge)
                )) {
                    return std::make_pair(
                        false,
                        "One or both atom indices ("s
                            + std::to_string(std::get<0>(edge))
                            + ", "
                            + std::to_string(std::get<1>(edge))
                            +") are invalid."
                    );
                }
                auto check_edge_in_adjacencies = [](
                    const vector<AtomIndexType>& haystack,
                    const AtomIndexType& needle
                ) -> bool {
                    return std::find(
                        haystack.begin(),
                        haystack.end(),
                        needle
                    ) != haystack.end();
                };

                // check forwards
                if(!check_edge_in_adjacencies(
                    _adjacencies.get_adjacencies(std::get<0>(edge)),
                    std::get<1>(edge)
                )) {
                    return std::make_pair(
                        false,
                        "The edge("s
                            + std::to_string(std::get<0>(edge))
                            + ", "
                            + std::to_string(std::get<1>(edge))
                            + ") in _edges is not in _adjacencies (forward)"
                    );
                }

                // check backwards
                if(!check_edge_in_adjacencies(
                    _adjacencies.get_adjacencies(std::get<1>(edge)),
                    std::get<0>(edge)
                )) {
                    return std::make_pair(
                        false,
                        "The edge("s
                            + std::to_string(std::get<0>(edge))
                            + ", "
                            + std::to_string(std::get<1>(edge))
                            + ") in _edges is not in _adjacencies (backward)"
                    );
                }
            }
            
            // passed all
            return std::make_pair(
                true,
                "Passed all tests"
            );
        }
    };
}
