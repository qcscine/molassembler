#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <functional>
#include <iostream>
#include <cassert>

/* TODO
 * - Nothing ensures that the molecular graph stays connected. 
 *   Either store chords or detect disconnects in remove_* functions
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

    /* Common typedefs */
    enum BondType : uint8_t {
        Single,
        Double,
        Triple,
        Quadruple,
        Quintuple,
        Sextuple
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
    using EdgeType = tuple<
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
        static bool _tuple_lex_lt_than_idxs(
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
        static bool _tuple_lex_gt_than_idxs(
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
        static bool tuple_lt(
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

        /*!
         * Checks whether a list of EdgeTypes is ordered.
         * \param 
         * \returns whether the list is ordered.
         * \note is O(E) 
         */
        static bool is_ordered(
            const std::vector<EdgeType>& container
        ) noexcept {
            /* this is a necessary underflow guard:
             * in the loop condition all types are unsigned:
             *  i < container.size() - 1
             * if container.size() is 0, the expression is an underflow
             */
            if(container.size() == 0) return true;

            for(EdgeIndexType i = 0; i < container.size() - 1; i++) {
                if(!tuple_lt(
                    container[i],
                    container[i+1]
                )) {
                    return false;
                }
            }

            return true;
        }

        /*!
         * Binary searches an ordered vector of edge tuples for a specified 
         * edge, returning whether it was found and the position in the 
         * container it was found at.
         * \param[in] container The vector containing the edge tuples.
         * \param[in] a The first atom index of the edge
         * \param[in] b The second atom index of the edge
         * \returns A pair containing whether it was found, and the position in
         *  the container it was found at. 
         * \pre a, b must fulfill a < b
         * \note Is O( log(N) ), where
         * - N is the number of edges stored in the container
         */
        static std::pair<bool, EdgeIndexType> binary_search(
            const std::vector<EdgeType>& container,
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept {
            /* important safeguard, otherwise
             *  unsigned_type right = container.size() - 1
             * is an underflow
             */
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

            // initialize L, M, R
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
                    // underflow safeguard
                    if(middle == 0) right = middle;
                    else right = middle - 1;
                    continue;
                } else { // not smaller and not greater -> equals
                    // found it!
                    return make_pair(
                        true, 
                        middle
                    );
                }
            }

            // search fails
            return make_pair(false, 0);
        }

        /*!
         * Binary searches an ordered vector of edge tuples and returns an
         * index before which the specified edge can be inserted.
         * \param[in] container The vector containing the edge tuples.
         * \param[in] a The first atom index of the edge
         * \param[in] b The second atom index of the edge
         * \returns The position within the container before which the
         *  specified edge can be inserted.
         * \pre a, b must fulfill a < b
         * \note Is O( log(N) ), where
         * - N is the number of edges stored in the container
         */
        static EdgeIndexType binary_search_insert_position(
            const vector<EdgeType>& container,
            const AtomIndexType& a, 
            const AtomIndexType& b
        ) noexcept {
            /* important safeguard, otherwise
             *  unsigned_type right = container.size() - 1;
             * causes an underflow
             */
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

            // initializations
            unsigned_type left = 0;
            unsigned_type right = container.size() - 1;
            unsigned_type middle;

            /* some debug help functions */
            /* auto print_state = [&](const unsigned_type& counter) -> void {
                std::cout << counter << ": " << left << " " << middle << " " << right << std::endl;
            }; */
            /* auto print_tuple = [](const EdgeType& tuple) -> void {
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
                    // underflow safeguard
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

    /*!
     * Stores information about interconnectivity between atoms.
     */
    class AbstractConnectivityManager {
    public:
        virtual ~AbstractConnectivityManager() {};

    /* Public member functions (in alphabetical order) */

        /*!
         * Add an atom to the ConnectivityManager.
         * \returns The index of the new atom.
         */
        virtual AtomIndexType add_atom() = 0;

        /*!
         * Adds a bond to the ConnectivityManager. 
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to register
         */
        virtual void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) = 0;

        /*!
         * Checks whether two atoms are bonded.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns Whether two atoms are bonded or not.
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const = 0;

        /* Checks whether a specified edge exists.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to check
         * \returns Whether the two atoms are bonded as specified.
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const = 0;

        /*!
         * Gets the BondType between two atoms.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns The bondtype between the atoms.
         */
        virtual BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const = 0;

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
        > get_bond_pairs(const AtomIndexType& a) const = 0;

        /*!
         * Get a vector of all bonded atom indices.
         * \param[in] a The atom index
         * \returns A vector all bonded atom indices.
         */
        virtual vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const = 0;

        /*!
         * Removes an atom from the ConnectivityManager. 
         * \param[in] a The atom index
         */
        virtual void remove_atom(
            const AtomIndexType& a
        ) = 0;

        /*!
         * Removes a bond from the ConnectivityManager. If the bond doesn't exist,
         * the function call is ignored.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         */
        virtual void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) = 0;

        virtual void reset() = 0;

        /*!
         * Validate the current state of the instance.
         * \returns A pair containing whether the state is valid, and an 
         *  explanation.
         */
        virtual std::pair<bool, std::string> validate() const noexcept = 0;

    };


    /* Implementation complexity comparison
     * ====================================
     *
     * in big-O notation
     * 
     * function name               Minimal   Fast
     * ------------------------    -------   -------
     * add_atom                    1         1
     * add_bond                    log E     log E
     * bond_exists (1)             log E     b_i
     * bond_exists (2)             log E     log E
     * get_bond_type               log E     log E
     * get_bond_pairs              E         log E
     * get_bonded_atom_indices     E         1
     * remove_atom                 E         log E
     * remove_bond                 log E     log E
     * validate                    E         E
     *
     * storage                     E         E
     *
     *
     * E is the number of stored edges
     * b_i is the number of stored edges for atom i, b_i << E,
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
     *
     *
     * Performance
     * ===========
     *
     * Benchmarking reveals that execution times cross for get_bond_pairs at 
     * roughly N = 300 atoms with E = 4N edges. Above 300 atoms, Fast is faster.
     * In general, graph modification operations are more expensive for Fast,
     * but access to the data is cheaper.
     *
     * Detailed graphs are in graphs/speed.pdf
     */

    class MinimalConnectivityManager : public AbstractConnectivityManager {
    private:
    /* Private members */
        std::vector<EdgeType> _edges;
        std::vector<bool> _atom_exists;

        /*!
         * Checks whether the supplied atom index is valid, i.e. in range.
         * \param[in] a The atom index to check
         * \returns whether the supplied atom index is valid, i.e. in range.
         */
        inline bool _valid_atom_index(const AtomIndexType& a) const noexcept {
            bool retval = (
                a < _atom_exists.size()
                && _atom_exists[a]
            );
            if(!retval) {
                std::cout << "Called _valid_atom_index with (" << a << "). " 
                    << "_atom_exists.size() = " << _atom_exists.size()
                    << ", _atom_exists[a] = " << (_atom_exists[a] ? "T" : "F" )
                    << std::endl;
            }
            return retval;
        }
        /*inline bool _valid_atom_index(const AtomIndexType& a) const noexcept {
            return (
                a < _atom_exists.size()
                && _atom_exists[a]
            );
        }*/

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
            bool retval = (
                _valid_atom_index(a)
                && _valid_atom_index(b)
                && a != b
            );
            if(!retval) {
                std::cout << "Called _valid_atom_indices with (" << a << ", " 
                    << b << ")." << std::endl;
            }
            return retval;
        }
        /*inline bool _valid_atom_indices(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept {
            return (
                _valid_atom_index(a)
                && _valid_atom_index(b)
                && a != b
            );
        }*/

    public:

        /*!
         * Adds an atom to the ConnectivityManager. In this implementation,
         * atoms exist only when bonds to them exist.
         * \returns The index of the new atom.
         * \note Is O(1)
         */
        AtomIndexType add_atom() noexcept override final {
            _atom_exists.emplace_back(true);
            return _atom_exists.size() - 1;
        }

        /*!
         * Adds a bond to the ConnectivityManager. 
         * \pre a, b must fulfill: a != b.
         * \note Is O( log(E) ), where
         *  - E is the number of stored edges
         */
        void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept override final {
            assert(_valid_atom_indices(a, b));

            auto pos = EdgeTypeHelpers::binary_search_insert_position(
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

        /*!
         * Checks whether two atoms are bonded.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns Whether two atoms are bonded or not.
         * \pre a, b must fulfill: a != b.
         * \note Is O( log(E) ), where
         * - E is the number of stored edges
         */
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final {
            assert(_valid_atom_indices(a, b));

            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
                _edges,
                std::min(a, b),
                std::max(a, b)
            );
            return found_and_pos_pair.first;
        }

        /*!
         * Checks whether two atoms are bonded by the specified BondType.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to check
         * \returns Whether the two atoms are bonded as specified.
         * \pre a, b must fulfill a != b
         * \note Is O( log(E) ), where
         * - E is the number of stored edges
         */
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept override final {
            assert(_valid_atom_indices(a, b));
            
            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
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

        /*!
         * Gets the BondType between two atoms.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns The bondtype between the atoms.
         * \pre A bond must exist for a, b for this to work correctly. In case 
         *  no bond exists for a,b this will erroneously return
         *  BondType::Single
         * \pre a, b must fulfill a != b
         * \note Is O( log(E) ), where
         * - E is the number of stored edges
         */
        BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final {
            assert(_valid_atom_indices(a, b));

            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
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

        /*!
         * Get a vector of all bonds outwards from an atom.
         * \param[in] a The atom index
         * \returns A vector of all bonds outwards from an atom as std::pairs,
         *  first is the BondType, second the atom index.
         * \note Is O(E), where
         * - E is the number of stored edges
         */
        vector<
            pair<
                AtomIndexType,
                BondType
            >
        > get_bond_pairs(const AtomIndexType& a) const noexcept override final {
            assert(_valid_atom_index(a));

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

        /*!
         * Get a vector of all bonded atom indices.
         * \param[in] a The atom index
         * \returns A vector all bonded atom indices.
         * \note is O(E), where
         * - E is the number of stored edges
         */
        vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const noexcept override final {
            assert(_valid_atom_index(a));

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

        /*!
         * Removes an atom from the ConnectivityManager. 
         * \param[in] a The atom index
         * \note is O(E), where
         * - E is the number of stored edges
         */
        void remove_atom(
            const AtomIndexType& a
        ) noexcept override final {
            assert(_valid_atom_index(a));

            // set the index to false
            _atom_exists[a] = false;

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
                ),
                _edges.end()
            );

        }

        /*!
         * Removes a bond from the ConnectivityManager. If the bond doesn't exist,
         * the function call is ignored.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \pre a, b must fulfill a != b 
         * \note is O( log(E) ), where
         * - E is the number of stored edges
         */
        void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept override final {
            assert(_valid_atom_indices(a, b));

            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
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

        void reset() override final {
            _edges = std::vector<EdgeType>();
            _atom_exists = std::vector<bool>();
        }

        /*!
         * A debug helper function. Checks that the internal state is valid.
         * \note is O(E)
         */
        std::pair<bool, std::string> validate() const noexcept override final {
            // is the edge list ordered
            if(!EdgeTypeHelpers::is_ordered(_edges)) {
                return make_pair(
                    false,
                    "The list of edges is unordered."
                );
            }

            // every edge contains valid atom indices
            for(const auto& edge: _edges) {
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
            }

            return make_pair(
                true,
                "Passed all checks"
            );
        }
    };

    class FastConnectivityManager : public AbstractConnectivityManager {
    private:
    /* Private members */
        std::vector<EdgeType> _edges;
        std::vector<
            std::vector<
                AtomIndexType
            >
        > _adjacencies;
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
        AtomIndexType add_atom() noexcept override final {
            // add an empty vector to _adjacencies at the end.
            _adjacencies.emplace_back();
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
        ) noexcept override final {
            assert(_valid_atom_indices(a, b));

            auto insert_position = EdgeTypeHelpers::binary_search_insert_position(
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
        ) const noexcept override final {
            assert(_valid_atom_indices(a, b));

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
        ) const noexcept override final {
            assert(_valid_atom_indices(a, b));

            /* for this case, checking the adjacency list is not good enough
             * the bond information is only in the edge list, binary search it
             */
            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
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
        ) const noexcept override final {
            assert(_valid_atom_indices(a, b));

            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
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
        > get_bond_pairs(const AtomIndexType& a) const noexcept override final {
            assert(_valid_atom_index(a));

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
                auto found_and_index_pair = EdgeTypeHelpers::binary_search(
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

        /*!
         * Get a vector of all bonded atom indices.
         * \param[in] a The atom index
         * \returns A vector all bonded atom indices.
         * \note is O(1)
         */
        vector<AtomIndexType> get_bonded_atom_indices(
            const AtomIndexType& a
        ) const noexcept override final {
            assert(_valid_atom_index(a));

            return _adjacencies[a];
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
        ) noexcept override final {
            assert(_valid_atom_index(a));

            // must update edges and adjacencies 
            auto bonded_to = _adjacencies[a];

            // erase all edges to and from this atom
            for(const auto& bonded_atom_index : bonded_to) {
                auto found_and_index_pair = EdgeTypeHelpers::binary_search(
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

            // remove all entries from this atom in adjacencies
            _adjacencies[a] = std::vector<AtomIndexType>();

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
        ) noexcept override final {
            assert(_valid_atom_indices(a, b));

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
            auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
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

        void reset() override final {
            _edges = std::vector<EdgeType>();
            _adjacencies = std::vector<
                std::vector<
                    AtomIndexType
                >
            >();
            _atom_exists = std::vector<bool>();
        }

        /*!
         * A debug helper function. Checks that the internal state is valid.
         * \note is O(E)
         */
        std::pair<bool, std::string> validate() const noexcept override final {
            // basic checks
            if(_adjacencies.size() != _atom_exists.size()) {
                return make_pair(
                    false,
                    "The size of _adjacencies and _atom_exists is mismatched."
                );
            }

            // is the edge list ordered
            if(!EdgeTypeHelpers::is_ordered(_edges)) {
                return make_pair(
                    false,
                    "The list of edges is unordered."
                );
            }

            /* all bonds specified in _adjacencies are present in _edges and
             *  all specified atom indices are valid
             */
            for(AtomIndexType from = 0; from < _adjacencies.size(); from++) {
                for(const auto& to : _adjacencies[from]) {
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
                    auto found_and_pos_pair = EdgeTypeHelpers::binary_search(
                        _edges,
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
                    _adjacencies[std::get<0>(edge)],
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
                    _adjacencies[std::get<1>(edge)],
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
