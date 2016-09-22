#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <functional>
#include <iostream>

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

    /* Common typedefs */
    enum BondType : uint8_t {
        Single,
        Double,
        Triple,
        Quadruple,
        Quintuple,
        Sextuple
    };
    using unsigned_type = unsigned;
    using EdgeIndexType = uint32_t;
    using AtomIndexType = uint16_t; // 65k max is sufficient

    /* derived types */
    using TupleType = tuple<
        AtomIndexType,
        AtomIndexType,
        BondType
    >;

    struct TupleTypeHelpers {
        static bool _tuple_lex_lt_than_idxs(
            const TupleType& tuple,
            const AtomIndexType& a,
            const AtomIndexType& b
        ) {
            return (
                std::get<0>(tuple) < a
                || (
                    std::get<0>(tuple) == a
                    && std::get<1>(tuple) < b 
                )
            );
        }

        static bool _tuple_lex_gt_than_idxs(
            const TupleType& tuple,
            const AtomIndexType& a,
            const AtomIndexType& b
        ) {
            return (
                std::get<0>(tuple) > a 
                || (
                    std::get<0>(tuple) == a 
                    && std::get<1>(tuple) > b 
                )
            );
        }

        static bool tuple_lt(
            const TupleType& a,
            const TupleType& b
        ) {
            return (
                std::get<0>(a) < std::get<0>(b)
                || (
                    std::get<0>(a) == std::get<0>(b)
                    && std::get<1>(a) < std::get<1>(b)
                )
            );
        }

        static std::pair<bool, EdgeIndexType> binary_search(
            const std::vector<TupleType>& container,
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept {
            if(container.size() == 0) return make_pair(false, 0);

            // bind the lexicographical comparison functions to the indices
            auto tuple_is_smaller = std::bind(
                _tuple_lex_lt_than_idxs,
                std::placeholders::_1,
                a,
                b
            );
            auto tuple_is_greater = std::bind(
                _tuple_lex_gt_than_idxs,
                std::placeholders::_1,
                a,
                b
            );

            unsigned_type left = 0;
            unsigned_type right = container.size() - 1;
            unsigned_type middle;

            while(left <= right) {
                middle = std::floor(
                    ((double) left + (double) right) / 2.0
                );

                if(tuple_is_smaller(container[middle])) {
                    left = middle + 1;
                    continue;
                } else if(tuple_is_greater(container[middle])) {
                    // Underflow guard
                    if(middle == 0) right = middle;
                    else right = middle - 1;
                    continue;
                } else { // implicit if tuple_equals(container[middle])
                    return make_pair(
                        true, 
                        middle
                    );
                }
            }

            return make_pair(false, 0);
        }

        static EdgeIndexType binary_search_insert_position(
            const vector<TupleType>& container,
            const AtomIndexType& a, 
            const AtomIndexType& b
        ) noexcept {
            if(container.size() == 0) {
                return 0;
            }

            auto tuple_is_smaller = std::bind(
                _tuple_lex_lt_than_idxs,
                std::placeholders::_1,
                a,
                b
            );
            auto tuple_is_greater = std::bind(
                _tuple_lex_gt_than_idxs,
                std::placeholders::_1,
                a,
                b
            );

            unsigned_type left = 0;
            unsigned_type right = container.size() - 1;
            unsigned_type middle;

            /* some debug help functions */
            /* auto print_state = [&](const unsigned_type& counter) -> void {
                std::cout << counter << ": " << left << " " << middle << " " << right << std::endl;
            }; */
            /* auto print_tuple = [](const TupleType& tuple) -> void {
                std::cout << "(" << std::get<0>(tuple) << ", ";
                std::cout << std::get<1>(tuple) << ", ";
                std::cout << std::get<2>(tuple) << ")" << std::endl;
            };
            */

            while(left <= right) {

                middle = std::floor(
                    ((double) left + (double) right) / 2.0
                );

                if(left == right) {
                    if(tuple_is_smaller(container[middle])) {
                        return middle + 1;
                    } else {
                        return middle;
                    }
                }

                if(tuple_is_smaller(container[middle])) {
                    left = middle + 1;
                    continue;
                }

                if(tuple_is_greater(container[middle])) {
                    // Must guard against underflow 
                    if(middle == 0) right = middle;
                    else right = middle - 1; 
                    continue;
                }
            }

            // "majority vote"
            if(middle == left) {
                return left;
            } else if(middle == right) {
                return right;
            } else {
                return 0;
            }
        }
    };

    class AdjacencyList {
    public:
        virtual ~AdjacencyList() {};

    /* Public member functions (in alphabetical order) */

        /*!
         * Add an atom to the AdjacencyList.
         */
        virtual void add_atom() noexcept = 0;

        /*!
         * Adds a bond to the AdjacencyList. No checks are made, and no
         * exceptions are thrown.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to register
         */
        virtual void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept = 0;

        /*!
         * Checks whether two atoms are bonded.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns Whether two atoms are bonded or not.
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept = 0;

        /*!
         * Checks whether two atoms are bonded by the specified BondType.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to check
         * \returns Whether the two atoms are bonded as specified.
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept = 0;

        /*!
         * Gets the BondType between two atoms.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns The bondtype between the atoms.
         */
        virtual BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept = 0;

        /*!
         * Get a vector of all bonds outwards from an atom.
         * \param[in] a The atom index
         * \returns A vector of all bonds outwards from an atom as std::pairs,
         * first is the BondType, second the atom index.
         */
        virtual vector<
            pair<
                AtomIndexType,
                BondType
            >
        > get_bond_pairs(const AtomIndexType& a) const noexcept = 0;

        /*!
         * Get a vector of all bonded atom indices.
         * \param[in] a The atom index
         * \returns A vector all bonded atom indices.
         */
        virtual vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const noexcept = 0;

        /*!
         * Checks whether the list is ordered.
         * \returns whether the list is ordered.
         */
        virtual bool is_ordered() const noexcept = 0;

        /*!
         * Removes an atom from the AdjacencyList. 
         * \param[in] a The atom index
         */
        virtual void remove_atom(
            const AtomIndexType& a
        ) noexcept = 0;

        /*!
         * Removes a bond from the AdjacencyList. If the bond doesn't exist,
         * the function call is ignored.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         */
        virtual void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept = 0;

    };


    /* Implementation complexity comparison
     * ====================================
     *
     *                             Minimal   Fast
     *                             -------   -------
     * add_atom                    1         1
     * add_bond                    log E     log E
     * bond_exists (1)             log E     log b 
     * bond_exists (2)             log E     log E
     * get_bond_type               log E     log E
     * get_bond_pairs              E         log E
     * get_bonded_atom_indices     E         1
     * remove_atom                 E         log E
     * remove_bond                 log E     log E
     *
     * storage                     E         ~3E
     *
     *
     * E is the number of stored edges
     * b is the number of stored edges for atom a, b << E
     *
     *
     * Summary of implementation differences
     * =====================================
     *
     * Minimal uses only an ordered edge list, containing tuples (i, j, bty),
     * where i and j are atom indices and bty is a bondtype. Fast uses the same
     * edge list and additionally a ragged two-dimensional vector containing 
     * atom indices (vector< vector<AtomIndexType> >). This helps speed up some
     * lookup algorithms considerably.
     *
     * All implementations make use of binary searches in the ordered edge
     * list.
     */

    class MinimalAdjacencyList : public AdjacencyList {
    private:
    /* Private members */
        std::vector<TupleType> _edges;

    public:
        /* Do nothing.
         * In this implementation, atoms exist only when bonds to them exist.
         * Is O(1)
         */
        void add_atom() noexcept {}

        /* a, b must fulfill: a != b.
         * Is O( log(E) ), where
         * - E is the number of stored edges
         */
        void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept override final {
            auto pos = TupleTypeHelpers::binary_search_insert_position(
                _edges,
                std::min(a, b), 
                std::max(a, b)
            );
            _edges.insert(
                _edges.begin() + pos,
                make_tuple(
                    std::min(a, b),
                    std::max(a, b),
                    bond_type
                )
            );
        }

        /*
         * Is O( log(E) ), where
         * - E is the number of stored edges
         */
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final {
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            return found_and_pos_pair.first;
        }

        /* a, b must fulfill a != b
         * Is O( log(E) ), where
         * - E is the number of stored edges
         */
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept override final {
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            return (
                found_and_pos_pair.first
                && (
                    bond_type == std::get<2>(
                        _edges[found_and_pos_pair.second]
                    )
                )
            );
        }

        /* A bond must exist for a, b for this to work correctly.
         * a, b must fulfill a != b
         * Is O( log(E) ), where
         * - E is the number of stored edges
         */
        BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final {
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            /* found_and_pos_pair.second is the index where the found tuple is
             *  stored
             * _edges[found_and_pos_pair.second] is the found tuple
             * std::get<2>(_edges[found_and_pos_pair.second]) is the 3rd
             *  element in the tuple -> the bondtype
             */
            return std::get<2>(
                _edges[found_and_pos_pair.second]
            );
        }

        /* Is O(E), where
         * - E is the number of stored edges
         */
        vector<
            pair<
                AtomIndexType,
                BondType
            >
        > get_bond_pairs(const AtomIndexType& a) const noexcept override final {
            vector<
                pair<
                    AtomIndexType,
                    BondType
                >
            > bond_list;

            auto it = _edges.begin();
            auto end = _edges.end();
            while(
                it != end 
                && std::get<0>(*it) <= a
            ) {
                if(std::get<0>(*it) == a) {
                    bond_list.emplace_back(
                        std::get<1>(*it),
                        std::get<2>(*it)
                    );
                } else if(std::get<1>(*it) == a) {
                    bond_list.emplace_back(
                        std::get<0>(*it),
                        std::get<2>(*it)
                    );
                }

                it++;
            }

            return bond_list;
        }

        /* is O(E), where
         * - E is the number of stored edges
         */
        vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const noexcept override final {
            vector<AtomIndexType> bonded_atom_indices;

            auto it = _edges.begin();
            auto end = _edges.end();
            while(
                it != end 
                && std::get<0>(*it) <= a
            ) {
                if(std::get<0>(*it) == a) {
                    bonded_atom_indices.emplace_back(
                        std::get<1>(*it)
                    );
                } else if(std::get<1>(*it) == a) {
                    bonded_atom_indices.emplace_back(
                        std::get<0>(*it)
                    );
                }

                it++;
            }

            return bonded_atom_indices;
        }

        bool is_ordered() const noexcept override final {
            auto lt = [](
                const TupleType& a,
                const TupleType& b
            ) {
                return (
                    std::get<0>(a) < std::get<0>(b)
                    || (
                        std::get<0>(a) == std::get<0>(b)
                        && std::get<1>(a) < std::get<1>(b)
                    )
                );
            };

            /* - Can't think of a way to put this in functional paradigm
             * - Not allowed to use indices if you want to be able to use list 
             *   in template instantiation
             * - List will be terrible anyway because random access is O(N)
             */
            for(EdgeIndexType i = 0; i < _edges.size() - 1; i++) {
                if(!lt(
                    _edges[i],
                    _edges[i+1]
                )) {
                    return false;
                }
            }
            return true;
        }

        /* is O(E), where
         * - E is the number of stored edges
         */
        void remove_atom(
            const AtomIndexType& a
        ) noexcept override final {
            /* check if the AtomIndexType is present in the edges, remove if so
             */

            _edges.erase(
                std::remove_if(
                    _edges.begin(),
                    _edges.end(),
                    [&a](const auto& tuple) {
                        return (
                            std::get<0>(tuple) == a
                            || std::get<1>(tuple) == a
                        );
                    }
                )
            );

        }

        /* a, b must fulfill a != b 
         * is O( log(E) ), where
         * - E is the number of stored edges
         */
        void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept override final {
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            if(found_and_pos_pair.first) {
                _edges.erase(
                    _edges.begin() + found_and_pos_pair.second
                );
            }
        }
    };

    class FastAdjacencyList : public AdjacencyList {
    private:
    /* Private members */
        std::vector<TupleType> _edges;
        std::vector<
            std::vector<
                AtomIndexType
            >
        > _adjacencies;

    public:
        /* is O(1) */
        void add_atom() noexcept override final {
            // add an empty vector to _adjacencies at the end.
            _adjacencies.emplace_back();
        };

        /* is O( log(E) ), where
         * - E is the number of stored edges
         */
        void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept override final {
            auto insert_position = TupleTypeHelpers::binary_search_insert_position(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );

            _edges.insert(
                _edges.begin() + insert_position,
                make_tuple(
                    a,
                    b,
                    bond_type
                )
            );

            _adjacencies[a].emplace_back(b);
            _adjacencies[b].emplace_back(a);
        }

        /* is O( B ), where
         * - B is the number of stored bonds for atom a
         */
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final {
            /* could binary search the edge list immediately, but linear search
             * in the adjacency list is probably faster if we don't need to 
             * know the bond type
             */
            return std::find(
                _adjacencies[a].begin(),
                _adjacencies[a].end(),
                b
            ) != _adjacencies[a].end();
        }

        /* is O( log(E) ), where
         * - E is the number of stored edges
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept {
            /* for this case, checking the adjacency list is not good enough
             * the bond information is only in the edge list, binary search it
             */
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            return (
                found_and_pos_pair.first
                && (
                    bond_type == std::get<2>(
                        _edges[found_and_pos_pair.second]
                    )
                )
            );
        }

        /* is O( log(E) ), where
         * - E is the number of stored edges
         */
        BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final {
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            /* found_and_pos_pair.second is the index where the found tuple is
             *  stored
             * _edges[found_and_pos_pair.second] is the found tuple
             * std::get<2>(_edges[found_and_pos_pair.second]) is the 3rd
             *  element in the tuple -> the bondtype
             */
            return std::get<2>(
                _edges[found_and_pos_pair.second]
            );
        }

        /* is O( B log (E) ), where
         * - B is the number of bonds from atom a
         * - E is the number of stored edges
         *
         * => O( log (E) ) with small constants
         */
        vector<
            pair<
                AtomIndexType,
                BondType
            >
        > get_bond_pairs(const AtomIndexType& a) const noexcept override final {
            /* We can get a full list of connections and how they are
             * ordered in the edge list from the adjacencies.
             */
            auto bonded_to = _adjacencies[a];
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
                auto found_and_index_pair = TupleTypeHelpers::binary_search(
                    _edges,
                    std::min(a, bonded_atom_index),
                    std::max(a, bonded_atom_index)
                );

                if(found_and_index_pair.first) {
                    bond_pairs.emplace_back(
                        bonded_atom_index,
                        std::get<2>(
                            _edges[found_and_index_pair.second]
                        )
                    );
                }
            }

            return bond_pairs;
        }

        /* is O(1) */
        vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const noexcept override final {
            return _adjacencies[a];
        }

        /* is O(E) */
        bool is_ordered() const noexcept override final {
            for(EdgeIndexType i = 0; i < _edges.size() - 1; i++) {
                if(!TupleTypeHelpers::tuple_lt(
                    _edges[i],
                    _edges[i+1]
                )) {
                    return false;
                }
            }
            return true;
        }

        /* is O( B log (E) + sum_all(e_i) ), where 
         * - B is the number of bonds the atom being removed has
         * - E is the number of stored edges
         * - e_i are the number of edges of the atoms bonded to the atom being
         *   removed
         *
         * => O( log (E) ) with some constants involved
         */
        void remove_atom(
            const AtomIndexType& a
        ) noexcept override final {
            // must update edges and adjacencies 
            auto bonded_to = _adjacencies[a];

            // erase all edges to and from this atom
            for(const auto& bonded_atom_index : bonded_to) {
                auto found_and_index_pair = TupleTypeHelpers::binary_search(
                    _edges,
                    std::min(a, bonded_atom_index),
                    std::max(a, bonded_atom_index)
                );

                if(found_and_index_pair.first) {
                    _edges.erase(
                        _edges.begin() + found_and_index_pair.second
                    );
                }
            }

            // remove all other mentions in _adjacencies
            for(const auto& bonded_atom_index : bonded_to) {
                _adjacencies[bonded_atom_index].erase(
                    std::remove(
                        _adjacencies[bonded_atom_index].begin(),
                        _adjacencies[bonded_atom_index].end(),
                        a
                    )
                );
            }

            // remove the atom from adjacencies
            _adjacencies.erase(
                _adjacencies.begin() + a
            );

        }

        /* is O( e_a + e_b + log (E) ), where
         * - e_i is the number of edges involving atom i
         * - E is the number of stored edges
         *
         * => O( log(E) ) with negligible constants
         */
        void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept {
            // must update _edges and _adjacencies
            _adjacencies[a].erase(
                std::remove(
                    _adjacencies[a].begin(),
                    _adjacencies[a].end(),
                    b
                )
            );
            _adjacencies[b].erase(
                std::remove(
                    _adjacencies[b].begin(),
                    _adjacencies[b].end(),
                    a
                )
            );
            auto found_and_pos_pair = TupleTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            if(found_and_pos_pair.first) {
                _edges.erase(
                    _edges.begin() + found_and_pos_pair.second
                );
            }
        };
    };
}
